from __future__ import annotations

import argparse
import urllib.request

import numpy as np
import pandas as pd

from config import DOC_TABLE_DIR, PROCESSED_DIR, RAW_DIR, SOURCE_WORKBOOK, TABLE_DIR, ZENODO_SOURCE_URL
from utils import ensure_dirs, short_taxon_name


def download_source() -> None:
    ensure_dirs(RAW_DIR)
    if SOURCE_WORKBOOK.exists():
        print(f"Found {SOURCE_WORKBOOK}")
        return
    print(f"Downloading {ZENODO_SOURCE_URL}")
    urllib.request.urlretrieve(ZENODO_SOURCE_URL, SOURCE_WORKBOOK)


def load_species_matrix(top_n: int, prevalence: float) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metadata = pd.read_excel(SOURCE_WORKBOOK, sheet_name="subject_metadata")
    metadata = metadata[["sample_name", "Case_status", "Age_at_collection", "Sex", "BMI", "total_sequences"]].copy()
    metadata = metadata.rename(columns={"Case_status": "diagnosis"})
    metadata = metadata[metadata["diagnosis"].isin(["PD", "Control"])].copy()
    metadata["label"] = (metadata["diagnosis"] == "PD").astype(int)

    taxa = pd.read_excel(SOURCE_WORKBOOK, sheet_name="metaphlan_rel_ab")
    taxa = taxa[taxa["clade_name"].astype(str).str.contains(r"\|s__|^s__", regex=True)].copy()
    taxa = taxa.set_index("clade_name")
    abundance = taxa.transpose()
    abundance.index.name = "sample_name"
    abundance = abundance.apply(pd.to_numeric, errors="coerce").fillna(0)
    if abundance.to_numpy().max() > 1:
        abundance = abundance / 100.0

    common = metadata["sample_name"][metadata["sample_name"].isin(abundance.index)]
    metadata = metadata.set_index("sample_name").loc[common].reset_index()
    abundance = abundance.loc[common]

    prevalence_mask = (abundance > 0).mean(axis=0) >= prevalence
    abundance = abundance.loc[:, prevalence_mask]
    variable_features = abundance.var(axis=0).sort_values(ascending=False).head(top_n).index
    abundance = abundance.loc[:, variable_features]

    transformed = np.log1p(abundance * 1_000_000)
    transformed.columns = [short_taxon_name(c) for c in transformed.columns]
    abundance.columns = transformed.columns

    feature_dict = pd.DataFrame(
        {
            "feature": transformed.columns,
            "mean_relative_abundance": abundance.mean(axis=0).to_numpy(),
            "prevalence": (abundance > 0).mean(axis=0).to_numpy(),
        }
    ).sort_values("mean_relative_abundance", ascending=False)
    return metadata, abundance, transformed, feature_dict


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n", type=int, default=300)
    parser.add_argument("--prevalence", type=float, default=0.05)
    args = parser.parse_args()

    ensure_dirs(PROCESSED_DIR, TABLE_DIR, DOC_TABLE_DIR)
    download_source()
    metadata, abundance, features, feature_dict = load_species_matrix(args.top_n, args.prevalence)

    metadata.to_csv(PROCESSED_DIR / "metadata.csv", index=False)
    abundance.to_csv(PROCESSED_DIR / "species_relative_abundance.csv")
    features.to_csv(PROCESSED_DIR / "species_ml_features.csv")
    feature_dict.to_csv(TABLE_DIR / "feature_dictionary.csv", index=False)
    feature_dict.head(50).to_csv(DOC_TABLE_DIR / "top_species_features.csv", index=False)

    manifest = pd.DataFrame(
        [
            {
                "dataset": "Wallen_2022_PD_gut_metagenomics",
                "source": "Zenodo 10.5281/zenodo.7246185",
                "n_samples": metadata.shape[0],
                "n_pd": int((metadata["diagnosis"] == "PD").sum()),
                "n_control": int((metadata["diagnosis"] == "Control").sum()),
                "n_species_features": features.shape[1],
            }
        ]
    )
    manifest.to_csv(TABLE_DIR / "dataset_manifest.csv", index=False)
    manifest.to_csv(DOC_TABLE_DIR / "dataset_manifest.csv", index=False)
    print(manifest)


if __name__ == "__main__":
    main()
