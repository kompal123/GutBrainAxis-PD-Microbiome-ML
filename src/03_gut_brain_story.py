from __future__ import annotations

import pandas as pd

from config import DOC_FIGURE_DIR, DOC_TABLE_DIR, FIGURE_DIR, TABLE_DIR
from utils import ensure_dirs, save_fig

import seaborn as sns


MECHANISM_RULES = {
    "SCFA / anti-inflammatory": ["Faecalibacterium", "Roseburia", "Eubacterium", "Butyrivibrio", "Agathobacter", "Anaerostipes"],
    "Mucin / gut barrier": ["Akkermansia", "Bacteroides", "Bifidobacterium"],
    "Opportunistic / inflammatory": ["Escherichia", "Klebsiella", "Enterococcus", "Streptococcus", "Bilophila", "Desulfovibrio"],
    "Neuroactive metabolism hypothesis": ["Bacteroides", "Parabacteroides", "Alistipes", "Lactobacillus", "Clostridium"],
    "Oral-gut translocation": ["Streptococcus", "Veillonella", "Rothia", "Actinomyces"],
}


def annotate_taxon(taxon: str) -> str:
    matches = [label for label, keys in MECHANISM_RULES.items() if any(key.lower() in taxon.lower() for key in keys)]
    return "; ".join(matches) if matches else "Unmapped / literature follow-up"


def main() -> None:
    ensure_dirs(TABLE_DIR, FIGURE_DIR, DOC_TABLE_DIR, DOC_FIGURE_DIR)
    diff = pd.read_csv(TABLE_DIR / "differential_species.csv")
    imp = pd.read_csv(TABLE_DIR / "ml_taxa_importance.csv")
    story = diff.merge(imp, on="taxon", how="left")
    story["direction"] = story["log2_fold_change_PD_vs_control"].map(lambda x: "PD-enriched" if x > 0 else "Control-enriched")
    story["gut_brain_axis_hypothesis"] = story["taxon"].map(annotate_taxon)
    story = story.sort_values(["importance", "p_value"], ascending=[False, True])
    story.to_csv(TABLE_DIR / "gut_brain_axis_taxa_story.csv", index=False)
    story.head(50).to_csv(DOC_TABLE_DIR / "gut_brain_axis_taxa_story.csv", index=False)

    plot_df = story.head(50).copy()
    exploded = plot_df.assign(gut_brain_axis_hypothesis=plot_df["gut_brain_axis_hypothesis"].str.split("; ")).explode("gut_brain_axis_hypothesis")
    counts = exploded["gut_brain_axis_hypothesis"].value_counts().reset_index()
    counts.columns = ["mechanism", "n_top_taxa"]
    fig = sns.barplot(data=counts, x="n_top_taxa", y="mechanism", color="#8d6e63")
    fig.set_xlabel("Number of top taxa")
    fig.set_ylabel("Gut-brain axis mechanism")
    fig.set_title("Biological themes among top PD microbiome signals")
    save_fig(FIGURE_DIR / "gut_brain_mechanism_summary.png", close=False)
    save_fig(DOC_FIGURE_DIR / "gut_brain_mechanism_summary.png")
    print(story.head(15)[["taxon", "direction", "importance", "gut_brain_axis_hypothesis"]])


if __name__ == "__main__":
    main()
