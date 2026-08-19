from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from config import DOC_FIGURE_DIR, DOC_TABLE_DIR, FIGURE_DIR, PROCESSED_DIR, TABLE_DIR
from utils import ensure_dirs, save_fig

import seaborn as sns


def bh_fdr(pvalues: list[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def main() -> None:
    ensure_dirs(TABLE_DIR, FIGURE_DIR, DOC_TABLE_DIR, DOC_FIGURE_DIR)
    meta = pd.read_csv(PROCESSED_DIR / "metadata.csv")
    X_log = pd.read_csv(PROCESSED_DIR / "species_ml_features.csv", index_col=0)
    X_raw = pd.read_csv(PROCESSED_DIR / "species_relative_abundance.csv", index_col=0)
    meta = meta.set_index("sample_name").loc[X_log.index]
    pd_mask = meta["diagnosis"] == "PD"

    rows = []
    for feature in X_log.columns:
        pd_vals = X_log.loc[pd_mask, feature]
        control_vals = X_log.loc[~pd_mask, feature]
        stat = mannwhitneyu(pd_vals, control_vals, alternative="two-sided")
        pd_mean = X_raw.loc[pd_mask, feature].mean()
        control_mean = X_raw.loc[~pd_mask, feature].mean()
        log2_fc = np.log2((pd_mean + 1e-8) / (control_mean + 1e-8))
        rows.append(
            {
                "taxon": feature,
                "pd_mean_relative_abundance": pd_mean,
                "control_mean_relative_abundance": control_mean,
                "log2_fold_change_PD_vs_control": log2_fc,
                "p_value": stat.pvalue,
            }
        )
    out = pd.DataFrame(rows)
    out["fdr"] = bh_fdr(out["p_value"].tolist())
    out = out.sort_values("p_value")
    out.to_csv(TABLE_DIR / "differential_species.csv", index=False)
    out.head(50).to_csv(DOC_TABLE_DIR / "top_differential_species.csv", index=False)

    top = out.reindex(out["log2_fold_change_PD_vs_control"].abs().sort_values(ascending=False).head(20).index).copy()
    top = top.sort_values("log2_fold_change_PD_vs_control")
    ax = sns.barplot(data=top, x="log2_fold_change_PD_vs_control", y="taxon", hue=top["log2_fold_change_PD_vs_control"] > 0, dodge=False)
    ax.legend_.remove()
    ax.set_xlabel("log2 fold change: PD vs control")
    ax.set_ylabel("Species")
    ax.set_title("Gut microbial species shifted in Parkinson's disease")
    save_fig(FIGURE_DIR / "differential_species_log2fc.png", close=False)
    save_fig(DOC_FIGURE_DIR / "differential_species_log2fc.png")
    print(out.head(10))


if __name__ == "__main__":
    main()
