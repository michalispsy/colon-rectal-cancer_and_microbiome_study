#!/usr/bin/env python3
"""
STEP 5.1 — SHAP ANALYSIS PER MODEL
===================================

Implements:
  5.1a All Universal SHAP
       - Top features per LODO fold
       - Consensus ranking / stability score
       - PA vs rCLR contribution split

  5.1b OnlyMale vs OnlyFemale SHAP
       - Male-only vs Female-only top features
       - Overlap analysis: shared, male-specific, female-specific

Expected Phase 4 input structure:
  OUTPUTS/phase_4/4.1/4.1_prep/Model_A/<Study>/X_train.csv
  OUTPUTS/phase_4/4.1/4.1_prep/Model_A/<Study>/y_train.csv
  OUTPUTS/phase_4/4.1/4.1_prep/Model_A/<Study>/X_test.csv
  OUTPUTS/phase_4/4.1/4.1_prep/Model_A/<Study>/y_test.csv

The script retrains the Phase 4 classifiers fold-by-fold using the same
prepared train/test matrices and then computes SHAP values on the held-out
LODO test fold.

Run from repository root, e.g.:
  python step5_1_shap_analysis.py \
    --input-dir OUTPUTS/phase_4/4.1/4.1_prep \
    --output-dir OUTPUTS/phase_5/5.1_shap \
    --top-k 20 \
    --background-size 200 \
    --explain-size 300

Notes:
  - Requires: pandas, numpy, scikit-learn, matplotlib, shap
  - If shap is missing: pip install shap
"""

from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

try:
    import shap
except ImportError as exc:
    raise ImportError(
        "The package 'shap' is required for Step 5.1. Install it with: pip install shap"
    ) from exc

import matplotlib
matplotlib.use("Agg")  # non-interactive backend; avoids Tkinter errors on headless/terminal runs
import matplotlib.pyplot as plt

RANDOM_STATE = 42

# Use the best/most useful Phase 4 algorithm choices by default.
# All can also be changed to RandomForest if you prefer a purely non-linear SHAP analysis.
DEFAULT_MODEL_ALGORITHMS = {
    "All": "LASSO",
    "OnlyMale": "RandomForest",
    "OnlyFemale": "RandomForest",
}

# Phase 4 prepared matrices may still be stored under the original directory names.
# The keys below are the names used in Phase 5 outputs/reports.
# The values are the source directories created by Phase 4.
MODEL_SOURCE_DIRS = {
    "All": "Model_A",
    "OnlyMale": "Model_B",
    "OnlyFemale": "Model_C",
    "Adenoma": "Model_D",
}


def make_classifier(algorithm: str):
    """Return a classifier matching Phase 4 settings."""
    if algorithm == "RandomForest":
        return RandomForestClassifier(
            n_estimators=500,
            random_state=RANDOM_STATE,
            class_weight="balanced",
            n_jobs=-1,
        )
    if algorithm == "LASSO":
        return LogisticRegression(
            penalty="l1",
            solver="liblinear",
            random_state=RANDOM_STATE,
            class_weight="balanced",
            max_iter=1000,
        )
    raise ValueError(f"Unknown algorithm: {algorithm}")


def infer_feature_type(feature: str) -> str:
    if feature.startswith("PA_"):
        return "PA"
    if feature.startswith("rCLR_"):
        return "rCLR"
    return "covariate"


def infer_species(feature: str) -> str:
    if feature.startswith("PA_"):
        return feature.replace("PA_", "", 1)
    if feature.startswith("rCLR_"):
        return feature.replace("rCLR_", "", 1)
    return feature


def safe_auc(y_true: pd.Series, y_score: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return np.nan
    return float(roc_auc_score(y_true, y_score))


def sample_rows(df: pd.DataFrame, max_n: int, random_state: int = RANDOM_STATE) -> pd.DataFrame:
    """Subsample rows for faster SHAP while keeping all rows if n <= max_n."""
    if max_n <= 0 or len(df) <= max_n:
        return df
    return df.sample(n=max_n, random_state=random_state)


def extract_binary_shap_values(raw_shap_values) -> np.ndarray:
    """
    Convert SHAP output to a 2D array: samples x features for the positive class.

    SHAP return formats vary by explainer/model/version:
      - list of arrays [class0, class1]
      - array with shape (samples, features)
      - array with shape (samples, features, classes)
      - Explanation object with .values
    """
    values = getattr(raw_shap_values, "values", raw_shap_values)

    if isinstance(values, list):
        if len(values) == 2:
            return np.asarray(values[1])
        return np.asarray(values[0])

    values = np.asarray(values)

    if values.ndim == 3:
        # Common new format: samples x features x classes
        if values.shape[-1] == 2:
            return values[:, :, 1]
        # Less common: classes x samples x features
        if values.shape[0] == 2:
            return values[1, :, :]

    if values.ndim == 2:
        return values

    raise ValueError(f"Unsupported SHAP values shape: {values.shape}")


def compute_shap_matrix(
    clf,
    algorithm: str,
    X_train: pd.DataFrame,
    X_explain: pd.DataFrame,
    background_size: int,
) -> np.ndarray:
    """Compute SHAP values for the positive class."""
    if algorithm == "RandomForest":
        explainer = shap.TreeExplainer(clf)
        raw_values = explainer.shap_values(X_explain, check_additivity=False)
        return extract_binary_shap_values(raw_values)

    if algorithm == "LASSO":
        # For linear models, use a small background dataset from train.
        background = sample_rows(X_train, background_size)
        explainer = shap.LinearExplainer(clf, background)
        raw_values = explainer.shap_values(X_explain)
        return extract_binary_shap_values(raw_values)

    raise ValueError(f"Unsupported algorithm for SHAP: {algorithm}")


def fold_feature_importance(
    shap_values: np.ndarray,
    X_explain: pd.DataFrame,
    model_name: str,
    algorithm: str,
    study: str,
    auc: float,
) -> pd.DataFrame:
    """Create per-fold feature importance table."""
    if shap_values.shape[1] != X_explain.shape[1]:
        raise ValueError(
            f"SHAP feature dimension mismatch for {model_name}/{study}: "
            f"{shap_values.shape[1]} vs {X_explain.shape[1]}"
        )

    mean_abs = np.abs(shap_values).mean(axis=0)
    mean_signed = shap_values.mean(axis=0)

    out = pd.DataFrame(
        {
            "Model": model_name,
            "Algorithm": algorithm,
            "Study": study,
            "Fold_AUC": auc,
            "Feature": X_explain.columns,
            "Feature_Type": [infer_feature_type(c) for c in X_explain.columns],
            "Species": [infer_species(c) for c in X_explain.columns],
            "MeanAbsSHAP": mean_abs,
            "MeanSignedSHAP": mean_signed,
        }
    )
    out = out.sort_values("MeanAbsSHAP", ascending=False).reset_index(drop=True)
    out["Rank"] = np.arange(1, len(out) + 1)
    return out


def run_model_shap(
    input_dir: Path,
    output_dir: Path,
    model_name: str,
    algorithm: str,
    top_k: int,
    background_size: int,
    explain_size: int,
) -> pd.DataFrame:
    """Run fold-by-fold SHAP for one model."""
    source_model_name = MODEL_SOURCE_DIRS.get(model_name, model_name)
    model_dir = input_dir / source_model_name
    if not model_dir.exists():
        raise FileNotFoundError(
            f"Missing model directory: {model_dir}. "
            f"For output model name '{model_name}', expected Phase 4 source directory '{source_model_name}'."
        )

    model_out = output_dir / model_name
    model_out.mkdir(parents=True, exist_ok=True)

    fold_tables: List[pd.DataFrame] = []
    studies = sorted([p.name for p in model_dir.iterdir() if p.is_dir()])

    for study in studies:
        fold_dir = model_dir / study
        required = ["X_train.csv", "y_train.csv", "X_test.csv", "y_test.csv"]
        if not all((fold_dir / f).exists() for f in required):
            print(f"[SKIP] {model_name}/{study}: missing required files")
            continue

        X_train = pd.read_csv(fold_dir / "X_train.csv", index_col=0)
        y_train = pd.read_csv(fold_dir / "y_train.csv", index_col=0)["Label"]
        X_test = pd.read_csv(fold_dir / "X_test.csv", index_col=0)
        y_test = pd.read_csv(fold_dir / "y_test.csv", index_col=0)["Label"]

        # Keep column order identical and numeric.
        X_test = X_test[X_train.columns]
        X_train = X_train.apply(pd.to_numeric)
        X_test = X_test.apply(pd.to_numeric)

        clf = make_classifier(algorithm)
        clf.fit(X_train, y_train)

        y_score = clf.predict_proba(X_test)[:, 1]
        auc = safe_auc(y_test, y_score)

        # Explain held-out samples. Subsample only for computational speed.
        X_explain = sample_rows(X_test, explain_size, random_state=RANDOM_STATE)
        shap_values = compute_shap_matrix(
            clf=clf,
            algorithm=algorithm,
            X_train=X_train,
            X_explain=X_explain,
            background_size=background_size,
        )

        fold_imp = fold_feature_importance(
            shap_values=shap_values,
            X_explain=X_explain,
            model_name=model_name,
            algorithm=algorithm,
            study=study,
            auc=auc,
        )
        fold_tables.append(fold_imp)

        fold_imp.to_csv(model_out / f"{model_name}_{algorithm}_{study}_all_features_shap.csv", index=False)
        fold_imp.head(top_k).to_csv(model_out / f"{model_name}_{algorithm}_{study}_top{top_k}_shap.csv", index=False)
        print(f"[OK] {model_name}/{study}/{algorithm}: AUC={auc:.3f}, explained n={len(X_explain)}")

    if not fold_tables:
        raise RuntimeError(f"No SHAP fold tables were created for {model_name}")

    all_folds = pd.concat(fold_tables, ignore_index=True)
    all_folds.to_csv(model_out / f"{model_name}_{algorithm}_all_folds_shap.csv", index=False)
    return all_folds


def consensus_ranking(all_folds: pd.DataFrame, top_k: int) -> pd.DataFrame:
    """
    Build consensus feature ranking across LODO folds.

    StabilityScore = fraction of folds where feature appears in top_k.
    ConsensusScore rewards high mean SHAP and stable repeated selection.
    """
    n_folds = all_folds["Study"].nunique()
    top = all_folds[all_folds["Rank"] <= top_k].copy()

    summary = (
        all_folds.groupby(["Model", "Algorithm", "Feature", "Feature_Type", "Species"], as_index=False)
        .agg(
            MeanAbsSHAP_Mean=("MeanAbsSHAP", "mean"),
            MeanAbsSHAP_Median=("MeanAbsSHAP", "median"),
            MeanSignedSHAP_Mean=("MeanSignedSHAP", "mean"),
            MeanRank=("Rank", "mean"),
            MedianRank=("Rank", "median"),
            BestRank=("Rank", "min"),
            N_Folds_Evaluated=("Study", "nunique"),
        )
    )

    stable = (
        top.groupby("Feature", as_index=False)
        .agg(N_Folds_TopK=("Study", "nunique"))
    )
    summary = summary.merge(stable, on="Feature", how="left")
    summary["N_Folds_TopK"] = summary["N_Folds_TopK"].fillna(0).astype(int)
    summary["Total_Folds"] = n_folds
    summary["StabilityScore"] = summary["N_Folds_TopK"] / n_folds

    # Scale SHAP to 0-1 for a combined score.
    max_shap = summary["MeanAbsSHAP_Mean"].max()
    summary["MeanAbsSHAP_Scaled"] = summary["MeanAbsSHAP_Mean"] / max_shap if max_shap > 0 else 0
    summary["ConsensusScore"] = 0.65 * summary["StabilityScore"] + 0.35 * summary["MeanAbsSHAP_Scaled"]

    summary = summary.sort_values(
        ["ConsensusScore", "StabilityScore", "MeanAbsSHAP_Mean"],
        ascending=False,
    ).reset_index(drop=True)
    summary["ConsensusRank"] = np.arange(1, len(summary) + 1)
    return summary


def summarize_feature_types(consensus: pd.DataFrame, top_n: int) -> pd.DataFrame:
    top = consensus.head(top_n).copy()
    out = (
        top.groupby("Feature_Type", as_index=False)
        .agg(
            N_Features=("Feature", "count"),
            MeanAbsSHAP_Total=("MeanAbsSHAP_Mean", "sum"),
            MeanAbsSHAP_Mean=("MeanAbsSHAP_Mean", "mean"),
        )
        .sort_values("MeanAbsSHAP_Total", ascending=False)
    )
    out["Percent_TopN"] = 100 * out["N_Features"] / len(top)
    out["Percent_SHAP_TopN"] = 100 * out["MeanAbsSHAP_Total"] / out["MeanAbsSHAP_Total"].sum()
    return out


def species_level_consensus(consensus: pd.DataFrame) -> pd.DataFrame:
    """Merge PA and rCLR entries into species-level importance."""
    microbe = consensus[consensus["Feature_Type"].isin(["PA", "rCLR"])].copy()
    if microbe.empty:
        return pd.DataFrame()

    out = (
        microbe.groupby(["Model", "Algorithm", "Species"], as_index=False)
        .agg(
            Species_MeanAbsSHAP=("MeanAbsSHAP_Mean", "sum"),
            Species_MaxStability=("StabilityScore", "max"),
            Species_BestRank=("ConsensusRank", "min"),
            PA_SHAP=("MeanAbsSHAP_Mean", lambda x: 0.0),
        )
    )

    # Add explicit PA/rCLR split.
    pivot = (
        microbe.pivot_table(
            index=["Model", "Algorithm", "Species"],
            columns="Feature_Type",
            values="MeanAbsSHAP_Mean",
            aggfunc="sum",
            fill_value=0.0,
        )
        .reset_index()
    )
    if "PA" not in pivot.columns:
        pivot["PA"] = 0.0
    if "rCLR" not in pivot.columns:
        pivot["rCLR"] = 0.0

    pivot = pivot.rename(columns={"PA": "PA_SHAP", "rCLR": "rCLR_SHAP"})
    out = out.drop(columns=["PA_SHAP"]).merge(pivot, on=["Model", "Algorithm", "Species"], how="left")
    out["Dominant_Signal"] = np.where(out["PA_SHAP"] >= out["rCLR_SHAP"], "PA", "rCLR")
    out = out.sort_values("Species_MeanAbsSHAP", ascending=False).reset_index(drop=True)
    out["SpeciesRank"] = np.arange(1, len(out) + 1)
    return out


def overlap_analysis(
    consensus_b: pd.DataFrame,
    consensus_c: pd.DataFrame,
    top_n: int,
    output_dir: Path,
) -> Dict[str, pd.DataFrame]:
    """Compare OnlyMale and OnlyFemale consensus features."""
    b_top = consensus_b.head(top_n).copy()
    c_top = consensus_c.head(top_n).copy()

    b_set = set(b_top["Feature"])
    c_set = set(c_top["Feature"])

    shared = sorted(b_set & c_set)
    male_specific = sorted(b_set - c_set)
    female_specific = sorted(c_set - b_set)

    def attach_info(features: Iterable[str], label: str) -> pd.DataFrame:
        rows = []
        for f in features:
            b_row = b_top[b_top["Feature"] == f]
            c_row = c_top[c_top["Feature"] == f]
            row = {
                "Category": label,
                "Feature": f,
                "Feature_Type": infer_feature_type(f),
                "Species": infer_species(f),
                "Male_Rank": int(b_row["ConsensusRank"].iloc[0]) if not b_row.empty else np.nan,
                "Female_Rank": int(c_row["ConsensusRank"].iloc[0]) if not c_row.empty else np.nan,
                "Male_ConsensusScore": float(b_row["ConsensusScore"].iloc[0]) if not b_row.empty else np.nan,
                "Female_ConsensusScore": float(c_row["ConsensusScore"].iloc[0]) if not c_row.empty else np.nan,
                "Male_Stability": float(b_row["StabilityScore"].iloc[0]) if not b_row.empty else np.nan,
                "Female_Stability": float(c_row["StabilityScore"].iloc[0]) if not c_row.empty else np.nan,
            }
            rows.append(row)
        return pd.DataFrame(rows)

    shared_df = attach_info(shared, "shared")
    male_df = attach_info(male_specific, "male_specific")
    female_df = attach_info(female_specific, "female_specific")

    combined = pd.concat([shared_df, male_df, female_df], ignore_index=True)
    combined.to_csv(output_dir / f"OnlyMale_vs_OnlyFemale_top{top_n}_feature_overlap.csv", index=False)

    stats = pd.DataFrame(
        [{
            "Top_N": top_n,
            "N_Male_Top": len(b_set),
            "N_Female_Top": len(c_set),
            "N_Shared": len(shared),
            "N_Male_Specific": len(male_specific),
            "N_Female_Specific": len(female_specific),
            "Jaccard_Index": len(shared) / len(b_set | c_set) if len(b_set | c_set) else np.nan,
            "Overlap_Percent_of_Male": 100 * len(shared) / len(b_set) if b_set else np.nan,
            "Overlap_Percent_of_Female": 100 * len(shared) / len(c_set) if c_set else np.nan,
        }]
    )
    stats.to_csv(output_dir / f"OnlyMale_vs_OnlyFemale_top{top_n}_overlap_stats.csv", index=False)

    return {
        "combined": combined,
        "stats": stats,
        "shared": shared_df,
        "male_specific": male_df,
        "female_specific": female_df,
    }


def plot_top_features(consensus: pd.DataFrame, output_path: Path, title: str, top_n: int = 20) -> None:
    df = consensus.head(top_n).iloc[::-1].copy()
    plt.figure(figsize=(10, max(5, 0.34 * len(df))))
    plt.barh(df["Feature"], df["ConsensusScore"])
    plt.xlabel("Consensus score")
    plt.ylabel("Feature")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_feature_type_split(type_summary: pd.DataFrame, output_path: Path, title: str) -> None:
    df = type_summary.copy()
    plt.figure(figsize=(7, 4.5))
    plt.bar(df["Feature_Type"], df["Percent_SHAP_TopN"])
    plt.ylabel("% of top-feature SHAP contribution")
    plt.xlabel("Feature type")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_overlap_counts(stats: pd.DataFrame, output_path: Path) -> None:
    row = stats.iloc[0]
    labels = ["Shared", "Male-specific", "Female-specific"]
    values = [row["N_Shared"], row["N_Male_Specific"], row["N_Female_Specific"]]
    plt.figure(figsize=(7, 4.5))
    plt.bar(labels, values)
    plt.ylabel("Number of features")
    plt.title(f"OnlyMale vs OnlyFemale top-{int(row['Top_N'])} SHAP overlap")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def write_report(
    output_dir: Path,
    model_consensus: Dict[str, pd.DataFrame],
    type_summaries: Dict[str, pd.DataFrame],
    overlap: Dict[str, pd.DataFrame],
    top_n: int,
) -> None:
    lines = []
    lines.append("# Step 5.1 — SHAP Analysis per Model\n")
    lines.append("## What was computed\n")
    lines.append("- Fold-level SHAP values on each held-out LODO test set.\n")
    lines.append("- Top features per fold.\n")
    lines.append("- Consensus ranking across folds using stability and mean absolute SHAP.\n")
    lines.append("- PA vs rCLR contribution split.\n")
    lines.append("- OnlyMale vs OnlyFemale overlap analysis.\n\n")

    for model_name, cons in model_consensus.items():
        algo = cons["Algorithm"].iloc[0]
        lines.append(f"## {model_name} ({algo})\n")
        lines.append(f"Top {top_n} consensus features:\n\n")
        cols = ["ConsensusRank", "Feature", "Feature_Type", "Species", "StabilityScore", "MeanAbsSHAP_Mean", "ConsensusScore"]
        lines.append(cons.head(top_n)[cols].to_markdown(index=False))
        lines.append("\n\nFeature-type split among top features:\n\n")
        lines.append(type_summaries[model_name].to_markdown(index=False))
        lines.append("\n\n")

    lines.append("## OnlyMale vs OnlyFemale overlap\n")
    lines.append(overlap["stats"].to_markdown(index=False))
    lines.append("\n\n")
    if not overlap["combined"].empty:
        lines.append(overlap["combined"].to_markdown(index=False))
        lines.append("\n")

    (output_dir / "step5_1_shap_report.md").write_text("".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 5.1 SHAP analysis for Phase 4 LODO models")
    parser.add_argument("--input-dir", default="OUTPUTS/phase_4/4.1/4.1_prep", help="Phase 4 LODO prep directory")
    parser.add_argument("--output-dir", default="OUTPUTS/phase_5/5.1_shap", help="Output directory")
    parser.add_argument("--top-k", type=int, default=20, help="Top K features per fold for stability score")
    parser.add_argument("--top-n-overlap", type=int, default=20, help="Top N consensus features for B vs C overlap")
    parser.add_argument("--background-size", type=int, default=200, help="Background samples for LinearExplainer")
    parser.add_argument("--explain-size", type=int, default=300, help="Max held-out samples explained per fold; <=0 means all")
    parser.add_argument(
        "--model-a-algorithm",
        choices=["LASSO", "RandomForest"],
        default=DEFAULT_MODEL_ALGORITHMS["All"],
        help="Algorithm to use for All SHAP",
    )
    parser.add_argument(
        "--model-b-algorithm",
        choices=["LASSO", "RandomForest"],
        default=DEFAULT_MODEL_ALGORITHMS["OnlyMale"],
        help="Algorithm to use for OnlyMale SHAP",
    )
    parser.add_argument(
        "--model-c-algorithm",
        choices=["LASSO", "RandomForest"],
        default=DEFAULT_MODEL_ALGORITHMS["OnlyFemale"],
        help="Algorithm to use for OnlyFemale SHAP",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    figures_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    selected = {
        "All": args.model_a_algorithm,
        "OnlyMale": args.model_b_algorithm,
        "OnlyFemale": args.model_c_algorithm,
    }

    model_consensus: Dict[str, pd.DataFrame] = {}
    type_summaries: Dict[str, pd.DataFrame] = {}

    for model_name, algorithm in selected.items():
        print(f"\n{'='*70}\nRunning SHAP: {model_name} / {algorithm}\n{'='*70}")
        all_folds = run_model_shap(
            input_dir=input_dir,
            output_dir=output_dir,
            model_name=model_name,
            algorithm=algorithm,
            top_k=args.top_k,
            background_size=args.background_size,
            explain_size=args.explain_size,
        )

        cons = consensus_ranking(all_folds, top_k=args.top_k)
        model_consensus[model_name] = cons
        cons.to_csv(output_dir / model_name / f"{model_name}_{algorithm}_consensus_ranking.csv", index=False)

        species_cons = species_level_consensus(cons)
        species_cons.to_csv(output_dir / model_name / f"{model_name}_{algorithm}_species_level_consensus.csv", index=False)

        type_summary = summarize_feature_types(cons, top_n=args.top_k)
        type_summaries[model_name] = type_summary
        type_summary.to_csv(output_dir / model_name / f"{model_name}_{algorithm}_PA_vs_rCLR_summary_top{args.top_k}.csv", index=False)

        plot_top_features(
            cons,
            figures_dir / f"{model_name}_{algorithm}_top{args.top_k}_consensus_features.png",
            title=f"{model_name} {algorithm}: top {args.top_k} SHAP consensus features",
            top_n=args.top_k,
        )
        plot_feature_type_split(
            type_summary,
            figures_dir / f"{model_name}_{algorithm}_PA_vs_rCLR_top{args.top_k}.png",
            title=f"{model_name} {algorithm}: PA vs rCLR contribution",
        )

    overlap = overlap_analysis(
        consensus_b=model_consensus["OnlyMale"],
        consensus_c=model_consensus["OnlyFemale"],
        top_n=args.top_n_overlap,
        output_dir=output_dir,
    )
    plot_overlap_counts(overlap["stats"], figures_dir / f"OnlyMale_vs_OnlyFemale_top{args.top_n_overlap}_overlap_counts.png")

    write_report(
        output_dir=output_dir,
        model_consensus=model_consensus,
        type_summaries=type_summaries,
        overlap=overlap,
        top_n=args.top_k,
    )

    print(f"\nDone. Results saved to: {output_dir}")
    print("Main files:")
    print(f"  - {output_dir / 'step5_1_shap_report.md'}")
    print(f"  - {output_dir / 'OnlyMale_vs_OnlyFemale_top20_feature_overlap.csv'}")
    print(f"  - {figures_dir}")


if __name__ == "__main__":
    main()
