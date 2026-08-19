from __future__ import annotations

import joblib
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import ConfusionMatrixDisplay, accuracy_score, average_precision_score, balanced_accuracy_score, classification_report, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_validate, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from config import DOC_FIGURE_DIR, DOC_TABLE_DIR, FIGURE_DIR, MODEL_DIR, PROCESSED_DIR, RANDOM_STATE, TABLE_DIR
from utils import ensure_dirs, save_fig

import matplotlib.pyplot as plt
import seaborn as sns


def models() -> dict:
    return {
        "logistic_regression": Pipeline(
            [("scaler", StandardScaler()), ("model", LogisticRegression(max_iter=3000, class_weight="balanced"))]
        ),
        "random_forest": RandomForestClassifier(n_estimators=400, min_samples_leaf=3, class_weight="balanced", random_state=RANDOM_STATE, n_jobs=1),
        "extra_trees": ExtraTreesClassifier(n_estimators=500, min_samples_leaf=2, class_weight="balanced", random_state=RANDOM_STATE, n_jobs=1),
    }


def main() -> None:
    ensure_dirs(TABLE_DIR, FIGURE_DIR, MODEL_DIR, DOC_TABLE_DIR, DOC_FIGURE_DIR)
    meta = pd.read_csv(PROCESSED_DIR / "metadata.csv").set_index("sample_name")
    X = pd.read_csv(PROCESSED_DIR / "species_ml_features.csv", index_col=0)
    y = meta.loc[X.index, "label"].astype(int)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    rows = []
    for name, model in models().items():
        scores = cross_validate(
            model,
            X,
            y,
            cv=cv,
            scoring={"roc_auc": "roc_auc", "average_precision": "average_precision", "balanced_accuracy": "balanced_accuracy", "f1": "f1"},
            n_jobs=1,
        )
        rows.append(
            {
                "model": name,
                "roc_auc_mean": scores["test_roc_auc"].mean(),
                "roc_auc_sd": scores["test_roc_auc"].std(),
                "pr_auc_mean": scores["test_average_precision"].mean(),
                "balanced_accuracy_mean": scores["test_balanced_accuracy"].mean(),
                "f1_mean": scores["test_f1"].mean(),
            }
        )
    summary = pd.DataFrame(rows).sort_values("roc_auc_mean", ascending=False)
    summary.to_csv(TABLE_DIR / "pd_classifier_benchmark.csv", index=False)
    summary.to_csv(DOC_TABLE_DIR / "pd_classifier_benchmark.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 4.2))
    sns.barplot(data=summary, x="roc_auc_mean", y="model", color="#3a7ca5", ax=ax)
    ax.set_xlim(0, 1)
    ax.set_xlabel("5-fold ROC-AUC")
    ax.set_ylabel("Model")
    ax.set_title("Parkinson's disease classifier from gut species")
    save_fig(FIGURE_DIR / "pd_classifier_benchmark.png")
    save_fig(DOC_FIGURE_DIR / "pd_classifier_benchmark.png")

    train_idx, test_idx = train_test_split(X.index, test_size=0.2, stratify=y, random_state=RANDOM_STATE)
    best_name = summary.iloc[0]["model"]
    best = models()[best_name]
    best.fit(X.loc[train_idx], y.loc[train_idx])
    pred = best.predict(X.loc[test_idx])
    prob = best.predict_proba(X.loc[test_idx])[:, 1]
    heldout = pd.DataFrame(
        [
            {
                "best_model": best_name,
                "accuracy": accuracy_score(y.loc[test_idx], pred),
                "balanced_accuracy": balanced_accuracy_score(y.loc[test_idx], pred),
                "f1": f1_score(y.loc[test_idx], pred),
                "roc_auc": roc_auc_score(y.loc[test_idx], prob),
                "pr_auc": average_precision_score(y.loc[test_idx], prob),
                "n_train": len(train_idx),
                "n_test": len(test_idx),
            }
        ]
    )
    heldout.to_csv(TABLE_DIR / "pd_heldout_metrics.csv", index=False)
    heldout.to_csv(DOC_TABLE_DIR / "pd_heldout_metrics.csv", index=False)
    pd.DataFrame(classification_report(y.loc[test_idx], pred, output_dict=True, zero_division=0)).transpose().to_csv(TABLE_DIR / "pd_heldout_report.csv")

    fig, ax = plt.subplots(figsize=(5.2, 4.7))
    ConfusionMatrixDisplay.from_predictions(y.loc[test_idx], pred, display_labels=["Control", "PD"], cmap="Blues", colorbar=False, ax=ax)
    ax.set_title("Held-out PD classification")
    save_fig(FIGURE_DIR / "pd_confusion_matrix.png")
    save_fig(DOC_FIGURE_DIR / "pd_confusion_matrix.png")

    perm = permutation_importance(best, X.loc[test_idx], y.loc[test_idx], n_repeats=20, random_state=RANDOM_STATE, n_jobs=1)
    importance = pd.DataFrame({"taxon": X.columns, "importance": perm.importances_mean}).sort_values("importance", ascending=False)
    importance.to_csv(TABLE_DIR / "ml_taxa_importance.csv", index=False)
    importance.head(50).to_csv(DOC_TABLE_DIR / "top_ml_taxa_importance.csv", index=False)
    top = importance.head(20).sort_values("importance")
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.barh(top["taxon"], top["importance"], color="#6a994e")
    ax.set_xlabel("Permutation importance")
    ax.set_ylabel("Species")
    ax.set_title("Microbial features driving PD prediction")
    save_fig(FIGURE_DIR / "ml_taxa_importance.png")
    save_fig(DOC_FIGURE_DIR / "ml_taxa_importance.png")

    joblib.dump(best, MODEL_DIR / "pd_microbiome_classifier.joblib")
    print(summary)
    print(heldout)


if __name__ == "__main__":
    main()
