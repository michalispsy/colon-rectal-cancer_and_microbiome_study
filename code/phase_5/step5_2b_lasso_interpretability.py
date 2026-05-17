#!/usr/bin/env python3
"""
STEP 5.3 — LASSO LOGISTIC REGRESSION INTERPRETABILITY
=====================================================

Purpose
-------
Explore the interpretability of sparse LASSO Logistic Regression across the
same Phase 4 LODO folds used for the Random Forest SHAP analysis.

Important terminology:
  - Random Forest / LASSO Logistic Regression = models / algorithms
  - All / OnlyMale / OnlyFemale = arrangements

Implemented arrangements:
  - All        = previous Model_A, universal CRC vs Control
  - OnlyMale   = previous Model_B, male-only CRC vs Control
  - OnlyFemale = previous Model_C, female-only CRC vs Control

This script intentionally ignores Adenoma / Models D,E,F unless you explicitly
extend it later.

What this script does
---------------------
For each arrangement and each LODO fold:
  1. Load the already prepared Phase 4 train/test matrices.
  2. Fit an L1-penalized logistic regression model.
  3. Optionally tune C inside the training set using Stratified CV.
  4. Evaluate the fold on the held-out study.
  5. Extract all coefficients and the non-zero selected features.

Then, across folds:
  6. Build a consensus coefficient ranking.
  7. Summarize selection stability, coefficient direction, PA/rCLR split.
  8. Build species-level consensus by merging PA and rCLR features.
  9. Compare LASSO-selected species with Random Forest SHAP species from 5.1.
 10. Optionally compare LASSO-selected species with Phase 3 formal differential
     abundance results.

Expected inputs
---------------
Phase 4 prepared matrices:
  OUTPUTS/phase_4/4.1/4.1_prep/Model_A/<Study>/X_train.csv
  OUTPUTS/phase_4/4.1/4.1_prep/Model_A/<Study>/y_train.csv
  ...

Optional Phase 5.1 Random Forest SHAP directory:
  OUTPUTS/phase_5/5.1_shap

Optional Phase 3 directory:
  OUTPUTS/phase_3

Example
-------
Run from repository root:

  python code/phase_5/step5_3_lasso_interpretability.py \
    --input-dir OUTPUTS/phase_4/4.1/4.1_prep \
    --rf-shap-dir OUTPUTS/phase_5/5.1_shap \
    --phase3-dir OUTPUTS/phase_3 \
    --output-dir OUTPUTS/phase_5/5.3_lasso_interpretability \
    --top-k 20 \
    --tune-c

Outputs
-------
  OUTPUTS/phase_5/5.3_lasso_interpretability/
    step5_3_lasso_interpretability_report.md
    lasso_model_performance_summary.csv
    lasso_arrangement_top20_overlap.csv
    lasso_arrangement_top20_overlap_stats.csv
    All/
      All_LASSO_all_folds_coefficients.csv
      All_LASSO_nonzero_coefficients_per_fold.csv
      All_LASSO_consensus_coefficients.csv
      All_LASSO_species_level_consensus.csv
      All_LASSO_feature_type_summary_top20.csv
      All_LASSO_top20_vs_RF_SHAP.csv              [if --rf-shap-dir exists]
      All_LASSO_top20_vs_Phase3.csv               [if --phase3-dir exists]
    OnlyMale/
      ...
    OnlyFemale/
      ...
    figures/
      All_LASSO_top20_consensus_coefficients.png
      All_LASSO_feature_type_top20.png
      OnlyMale_LASSO_top20_consensus_coefficients.png
      OnlyFemale_LASSO_top20_consensus_coefficients.png
      lasso_arrangement_top20_overlap.png
      lasso_performance_auc.png
"""

from __future__ import annotations

import argparse
import re
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score
from sklearn.model_selection import StratifiedKFold

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

RANDOM_STATE = 42

# Output/report names use arrangements.
# Source directories use Phase 4's historical names.
ARRANGEMENT_SOURCE_DIRS = {
    "All": "Model_A",
    "OnlyMale": "Model_B",
    "OnlyFemale": "Model_C",
}

PHASE3_FILTERS = {
    "All": {
        "3.1": {"comparison": {"CRC_vs_Control_all"}},
        "3.2": {"comparison": {"CRC_vs_Control__all"}},
        "3.3": {"comparison": {"CRC_vs_Control"}, "stratum": {"All"}},
    },
    "OnlyMale": {
        "3.1": {"comparison": {"CRC_vs_Control_Male"}},
        "3.2": {"comparison": {"CRC_vs_Control__Male"}},
        "3.3": {"comparison": {"CRC_vs_Control"}, "stratum": {"Male"}},
    },
    "OnlyFemale": {
        "3.1": {"comparison": {"CRC_vs_Control_Female"}},
        "3.2": {"comparison": {"CRC_vs_Control__Female"}},
        "3.3": {"comparison": {"CRC_vs_Control"}, "stratum": {"Female"}},
    },
}


def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def infer_feature_type(feature: str) -> str:
    if feature.startswith("PA_") or feature.startswith("PA__"):
        return "PA"
    if feature.startswith("rCLR_") or feature.startswith("rCLR__"):
        return "rCLR"
    return "covariate"


def infer_species(feature: str) -> str:
    s = str(feature)
    s = re.sub(r"^(PA|rCLR|CLR|clr|pa)__", "", s)
    s = re.sub(r"^(PA|rCLR|CLR|clr|pa)_", "", s)
    return s


def clean_species_name(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    s = re.sub(r"^(PA|rCLR|CLR|clr|pa)__", "", s)
    s = re.sub(r"^(PA|rCLR|CLR|clr|pa)_", "", s)
    s = s.replace("s__", "")
    s = s.replace("_", " ")
    s = re.sub(r"\s+", " ", s).strip()
    return s


def canonical_key(x: object) -> str:
    return clean_species_name(x).lower()


def safe_auc(y_true: pd.Series, y_score: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return np.nan
    return float(roc_auc_score(y_true, y_score))


def safe_auprc(y_true: pd.Series, y_score: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return np.nan
    return float(average_precision_score(y_true, y_score))


def safe_f1(y_true: pd.Series, y_pred: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return np.nan
    return float(f1_score(y_true, y_pred))


def build_lasso_model(
    tune_c: bool,
    c_value: float,
    cv_folds: int,
    c_grid: np.ndarray,
    random_state: int,
):
    """
    Create either:
      - LogisticRegressionCV with L1 penalty to tune C inside train folds, or
      - LogisticRegression with fixed C.

    The L1 penalty gives sparse coefficients, i.e. many exact zeros.
    """
    if tune_c:
        return LogisticRegressionCV(
            Cs=c_grid,
            cv=StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=random_state),
            penalty="l1",
            solver="liblinear",
            scoring="roc_auc",
            class_weight="balanced",
            max_iter=3000,
            random_state=random_state,
            n_jobs=1,
            refit=True,
        )

    return LogisticRegression(
        penalty="l1",
        solver="liblinear",
        C=c_value,
        class_weight="balanced",
        max_iter=3000,
        random_state=random_state,
    )


def get_model_c(clf) -> float:
    if hasattr(clf, "C_"):
        c = getattr(clf, "C_")
        try:
            return float(np.ravel(c)[0])
        except Exception:
            return np.nan
    if hasattr(clf, "C"):
        return float(getattr(clf, "C"))
    return np.nan


def load_fold_data(fold_dir: Path) -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame, pd.Series]:
    X_train = pd.read_csv(fold_dir / "X_train.csv", index_col=0)
    y_train = pd.read_csv(fold_dir / "y_train.csv", index_col=0)["Label"]
    X_test = pd.read_csv(fold_dir / "X_test.csv", index_col=0)
    y_test = pd.read_csv(fold_dir / "y_test.csv", index_col=0)["Label"]

    X_test = X_test[X_train.columns]
    X_train = X_train.apply(pd.to_numeric)
    X_test = X_test.apply(pd.to_numeric)
    return X_train, y_train, X_test, y_test


def run_arrangement_lasso(
    input_dir: Path,
    output_dir: Path,
    arrangement: str,
    top_k: int,
    coef_threshold: float,
    tune_c: bool,
    c_value: float,
    cv_folds: int,
    c_grid: np.ndarray,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Fit LASSO Logistic Regression for one arrangement across LODO folds."""
    source_dir = input_dir / ARRANGEMENT_SOURCE_DIRS[arrangement]
    if not source_dir.exists():
        raise FileNotFoundError(
            f"Missing Phase 4 source directory for {arrangement}: {source_dir}"
        )

    arrangement_out = output_dir / arrangement
    safe_mkdir(arrangement_out)

    all_coef_tables: List[pd.DataFrame] = []
    performance_rows: List[dict] = []

    studies = sorted([p.name for p in source_dir.iterdir() if p.is_dir()])
    for study in studies:
        fold_dir = source_dir / study
        required = ["X_train.csv", "y_train.csv", "X_test.csv", "y_test.csv"]
        if not all((fold_dir / f).exists() for f in required):
            print(f"[SKIP] {arrangement}/{study}: missing required files")
            continue

        X_train, y_train, X_test, y_test = load_fold_data(fold_dir)

        # CV cannot have more splits than the smallest class count.
        min_class_count = int(y_train.value_counts().min())
        effective_cv = max(2, min(cv_folds, min_class_count))
        clf = build_lasso_model(
            tune_c=tune_c,
            c_value=c_value,
            cv_folds=effective_cv,
            c_grid=c_grid,
            random_state=RANDOM_STATE,
        )
        clf.fit(X_train, y_train)

        y_score = clf.predict_proba(X_test)[:, 1]
        y_pred = (y_score >= 0.5).astype(int)

        auc = safe_auc(y_test, y_score)
        auprc = safe_auprc(y_test, y_score)
        f1 = safe_f1(y_test, y_pred)
        selected_c = get_model_c(clf)

        coef = np.ravel(clf.coef_)
        nonzero = np.abs(coef) > coef_threshold
        n_nonzero = int(nonzero.sum())

        fold_table = pd.DataFrame({
            "Arrangement": arrangement,
            "Model": "LASSO_LogisticRegression",
            "Study": study,
            "Selected_C": selected_c,
            "Fold_AUC": auc,
            "Fold_AUPRC": auprc,
            "Fold_F1": f1,
            "Feature": X_train.columns,
            "Feature_Type": [infer_feature_type(c) for c in X_train.columns],
            "Species": [infer_species(c) for c in X_train.columns],
            "Species_clean": [clean_species_name(infer_species(c)) for c in X_train.columns],
            "Species_key": [canonical_key(infer_species(c)) for c in X_train.columns],
            "Coefficient": coef,
            "AbsCoefficient": np.abs(coef),
            "IsNonZero": nonzero,
            "Direction": np.where(coef > coef_threshold, "CRC-associated",
                          np.where(coef < -coef_threshold, "Control-associated", "zero")),
        }).sort_values("AbsCoefficient", ascending=False).reset_index(drop=True)

        fold_table["AbsCoefRank"] = np.arange(1, len(fold_table) + 1)
        all_coef_tables.append(fold_table)

        fold_table.to_csv(
            arrangement_out / f"{arrangement}_LASSO_{study}_all_coefficients.csv",
            index=False,
        )
        fold_table[fold_table["IsNonZero"]].head(top_k).to_csv(
            arrangement_out / f"{arrangement}_LASSO_{study}_top{top_k}_nonzero_coefficients.csv",
            index=False,
        )

        performance_rows.append({
            "Arrangement": arrangement,
            "Study": study,
            "N_train": len(X_train),
            "N_test": len(X_test),
            "Selected_C": selected_c,
            "AUC": auc,
            "AUPRC": auprc,
            "F1": f1,
            "N_nonzero_features": n_nonzero,
            "Percent_nonzero_features": 100 * n_nonzero / X_train.shape[1],
            "N_features_total": X_train.shape[1],
            "CV_folds_used_for_C": effective_cv if tune_c else 0,
        })

        print(
            f"[OK] {arrangement}/{study}/LASSO: "
            f"AUC={auc:.3f}, C={selected_c:.4g}, nonzero={n_nonzero}/{X_train.shape[1]}"
        )

    if not all_coef_tables:
        raise RuntimeError(f"No LASSO coefficient tables were created for {arrangement}")

    all_coefs = pd.concat(all_coef_tables, ignore_index=True)
    perf = pd.DataFrame(performance_rows)

    all_coefs.to_csv(arrangement_out / f"{arrangement}_LASSO_all_folds_coefficients.csv", index=False)
    all_coefs[all_coefs["IsNonZero"]].to_csv(
        arrangement_out / f"{arrangement}_LASSO_nonzero_coefficients_per_fold.csv",
        index=False,
    )
    perf.to_csv(arrangement_out / f"{arrangement}_LASSO_fold_performance.csv", index=False)

    return all_coefs, perf


def consensus_coefficients(all_coefs: pd.DataFrame, top_k: int, coef_threshold: float) -> pd.DataFrame:
    """Build consensus feature ranking across LODO folds."""
    n_folds = all_coefs["Study"].nunique()
    nz = all_coefs[all_coefs["IsNonZero"]].copy()

    summary = (
        all_coefs.groupby(
            ["Arrangement", "Model", "Feature", "Feature_Type", "Species", "Species_clean", "Species_key"],
            as_index=False,
        )
        .agg(
            MeanCoefficient=("Coefficient", "mean"),
            MedianCoefficient=("Coefficient", "median"),
            MeanAbsCoefficient=("AbsCoefficient", "mean"),
            MedianAbsCoefficient=("AbsCoefficient", "median"),
            MaxAbsCoefficient=("AbsCoefficient", "max"),
            MeanAbsCoefRank=("AbsCoefRank", "mean"),
            BestAbsCoefRank=("AbsCoefRank", "min"),
            N_Folds_Evaluated=("Study", "nunique"),
        )
    )

    nz_summary = (
        nz.groupby("Feature", as_index=False)
        .agg(
            N_Folds_NonZero=("Study", "nunique"),
            N_Positive=("Coefficient", lambda x: int((x > coef_threshold).sum())),
            N_Negative=("Coefficient", lambda x: int((x < -coef_threshold).sum())),
            MeanNonZeroCoefficient=("Coefficient", "mean"),
            MeanAbsNonZeroCoefficient=("AbsCoefficient", "mean"),
        )
    )

    top_nz = (
        nz.sort_values(["Study", "AbsCoefficient"], ascending=[True, False])
        .groupby("Study")
        .head(top_k)
    )
    top_summary = (
        top_nz.groupby("Feature", as_index=False)
        .agg(N_Folds_TopK_NonZero=("Study", "nunique"))
    )

    summary = summary.merge(nz_summary, on="Feature", how="left")
    summary = summary.merge(top_summary, on="Feature", how="left")

    for col in ["N_Folds_NonZero", "N_Positive", "N_Negative", "N_Folds_TopK_NonZero"]:
        summary[col] = summary[col].fillna(0).astype(int)
    for col in ["MeanNonZeroCoefficient", "MeanAbsNonZeroCoefficient"]:
        summary[col] = summary[col].fillna(0.0)

    summary["Total_Folds"] = n_folds
    summary["SelectionFrequency"] = summary["N_Folds_NonZero"] / n_folds
    summary["TopKSelectionFrequency"] = summary["N_Folds_TopK_NonZero"] / n_folds
    summary["SignedConsistency"] = summary[["N_Positive", "N_Negative"]].max(axis=1) / summary["N_Folds_NonZero"].replace(0, np.nan)
    summary["SignedConsistency"] = summary["SignedConsistency"].fillna(0.0)

    max_abs_nonzero = summary["MeanAbsNonZeroCoefficient"].max()
    summary["MeanAbsNonZeroCoef_Scaled"] = (
        summary["MeanAbsNonZeroCoefficient"] / max_abs_nonzero if max_abs_nonzero > 0 else 0.0
    )

    # ConsensusScore favors repeated selection, strong coefficient magnitude among selected folds,
    # and stable coefficient direction.
    summary["ConsensusScore"] = (
        0.55 * summary["SelectionFrequency"]
        + 0.25 * summary["MeanAbsNonZeroCoef_Scaled"]
        + 0.20 * summary["SignedConsistency"]
    )

    summary["ConsensusDirection"] = np.where(
        summary["N_Positive"] > summary["N_Negative"],
        "CRC-associated",
        np.where(summary["N_Negative"] > summary["N_Positive"], "Control-associated", "mixed/zero"),
    )

    summary = summary.sort_values(
        ["ConsensusScore", "SelectionFrequency", "MeanAbsNonZeroCoefficient"],
        ascending=False,
    ).reset_index(drop=True)
    summary["ConsensusRank"] = np.arange(1, len(summary) + 1)
    return summary


def species_level_consensus(consensus: pd.DataFrame) -> pd.DataFrame:
    """Aggregate PA/rCLR features to species-level LASSO importance."""
    microbe = consensus[consensus["Feature_Type"].isin(["PA", "rCLR"])].copy()
    if microbe.empty:
        return pd.DataFrame()

    base = (
        microbe.groupby(["Arrangement", "Model", "Species_clean", "Species_key"], as_index=False)
        .agg(
            Species_MaxSelectionFrequency=("SelectionFrequency", "max"),
            Species_MaxTopKSelectionFrequency=("TopKSelectionFrequency", "max"),
            Species_SumMeanAbsNonZeroCoef=("MeanAbsNonZeroCoefficient", "sum"),
            Species_BestConsensusRank=("ConsensusRank", "min"),
            N_Selected_FeatureTypes=("Feature_Type", lambda x: int(x.nunique())),
            Any_CRC_Associated=("ConsensusDirection", lambda x: bool((x == "CRC-associated").any())),
            Any_Control_Associated=("ConsensusDirection", lambda x: bool((x == "Control-associated").any())),
        )
    )

    pivot_coef = microbe.pivot_table(
        index=["Arrangement", "Model", "Species_clean", "Species_key"],
        columns="Feature_Type",
        values="MeanAbsNonZeroCoefficient",
        aggfunc="sum",
        fill_value=0.0,
    ).reset_index()
    if "PA" not in pivot_coef.columns:
        pivot_coef["PA"] = 0.0
    if "rCLR" not in pivot_coef.columns:
        pivot_coef["rCLR"] = 0.0
    pivot_coef = pivot_coef.rename(columns={"PA": "PA_MeanAbsNonZeroCoef", "rCLR": "rCLR_MeanAbsNonZeroCoef"})

    pivot_freq = microbe.pivot_table(
        index=["Arrangement", "Model", "Species_clean", "Species_key"],
        columns="Feature_Type",
        values="SelectionFrequency",
        aggfunc="max",
        fill_value=0.0,
    ).reset_index()
    if "PA" not in pivot_freq.columns:
        pivot_freq["PA"] = 0.0
    if "rCLR" not in pivot_freq.columns:
        pivot_freq["rCLR"] = 0.0
    pivot_freq = pivot_freq.rename(columns={"PA": "PA_SelectionFrequency", "rCLR": "rCLR_SelectionFrequency"})

    out = base.merge(pivot_coef, on=["Arrangement", "Model", "Species_clean", "Species_key"], how="left")
    out = out.merge(pivot_freq, on=["Arrangement", "Model", "Species_clean", "Species_key"], how="left")
    out["Dominant_Signal"] = np.where(
        out["PA_MeanAbsNonZeroCoef"] >= out["rCLR_MeanAbsNonZeroCoef"],
        "PA",
        "rCLR",
    )
    out = out.sort_values(
        ["Species_MaxSelectionFrequency", "Species_SumMeanAbsNonZeroCoef"],
        ascending=False,
    ).reset_index(drop=True)
    out["SpeciesRank"] = np.arange(1, len(out) + 1)
    return out


def feature_type_summary(consensus: pd.DataFrame, top_n: int) -> pd.DataFrame:
    top = consensus.head(top_n).copy()
    out = (
        top.groupby("Feature_Type", as_index=False)
        .agg(
            N_Features=("Feature", "count"),
            MeanAbsNonZeroCoef_Total=("MeanAbsNonZeroCoefficient", "sum"),
            Mean_SelectionFrequency=("SelectionFrequency", "mean"),
        )
        .sort_values("MeanAbsNonZeroCoef_Total", ascending=False)
    )
    out["Percent_TopN"] = 100 * out["N_Features"] / len(top)
    total = out["MeanAbsNonZeroCoef_Total"].sum()
    out["Percent_Coefficient_TopN"] = 100 * out["MeanAbsNonZeroCoef_Total"] / total if total > 0 else 0.0
    return out


def arrangement_overlap(species_cons: Dict[str, pd.DataFrame], top_n: int, output_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sets = {}
    for arr, df in species_cons.items():
        sets[arr] = set(df.head(top_n)["Species_key"])

    rows = []
    all_arrangements = list(sets)
    for arr in all_arrangements:
        others = set().union(*[sets[o] for o in all_arrangements if o != arr])
        for key in sorted(sets[arr]):
            cats = [a for a in all_arrangements if key in sets[a]]
            rows.append({
                "Species_key": key,
                "Arrangements": ";".join(cats),
                "N_Arrangements": len(cats),
                "In_All": key in sets.get("All", set()),
                "In_OnlyMale": key in sets.get("OnlyMale", set()),
                "In_OnlyFemale": key in sets.get("OnlyFemale", set()),
                "Specificity": "shared" if len(cats) > 1 else f"{arr}_specific",
            })
    combined = pd.DataFrame(rows).drop_duplicates("Species_key").sort_values(["N_Arrangements", "Species_key"], ascending=[False, True])

    shared_all_three = sets["All"] & sets["OnlyMale"] & sets["OnlyFemale"]
    shared_all_male = sets["All"] & sets["OnlyMale"]
    shared_all_female = sets["All"] & sets["OnlyFemale"]
    shared_male_female = sets["OnlyMale"] & sets["OnlyFemale"]

    stats = pd.DataFrame([{
        "Top_N": top_n,
        "N_All_Top": len(sets["All"]),
        "N_OnlyMale_Top": len(sets["OnlyMale"]),
        "N_OnlyFemale_Top": len(sets["OnlyFemale"]),
        "N_Shared_All_Three": len(shared_all_three),
        "N_Shared_All_OnlyMale": len(shared_all_male),
        "N_Shared_All_OnlyFemale": len(shared_all_female),
        "N_Shared_OnlyMale_OnlyFemale": len(shared_male_female),
        "Jaccard_All_OnlyMale": len(shared_all_male) / len(sets["All"] | sets["OnlyMale"]) if sets["All"] | sets["OnlyMale"] else np.nan,
        "Jaccard_All_OnlyFemale": len(shared_all_female) / len(sets["All"] | sets["OnlyFemale"]) if sets["All"] | sets["OnlyFemale"] else np.nan,
        "Jaccard_OnlyMale_OnlyFemale": len(shared_male_female) / len(sets["OnlyMale"] | sets["OnlyFemale"]) if sets["OnlyMale"] | sets["OnlyFemale"] else np.nan,
    }])

    combined.to_csv(output_dir / f"lasso_arrangement_top{top_n}_overlap.csv", index=False)
    stats.to_csv(output_dir / f"lasso_arrangement_top{top_n}_overlap_stats.csv", index=False)
    return combined, stats


def find_species_consensus_file(shap_dir: Path, arrangement: str) -> Optional[Path]:
    """Find Random Forest SHAP species-level consensus from Step 5.1."""
    candidates = []
    # Expected new names.
    candidates.extend(sorted((shap_dir / arrangement).glob(f"{arrangement}_RandomForest_species_level_consensus.csv")))
    candidates.extend(sorted((shap_dir / arrangement).glob("*RandomForest*species_level_consensus.csv")))
    # Fallback old names if needed.
    old = {"All": "Model_A", "OnlyMale": "Model_B", "OnlyFemale": "Model_C"}[arrangement]
    candidates.extend(sorted((shap_dir / old).glob("*RandomForest*species_level_consensus.csv")))
    candidates.extend(sorted(shap_dir.glob(f"**/{arrangement}*RandomForest*species_level_consensus.csv")))
    candidates.extend(sorted(shap_dir.glob(f"**/{old}*RandomForest*species_level_consensus.csv")))
    for p in candidates:
        if p.exists():
            return p
    return None


def compare_with_rf_shap(
    lasso_species: pd.DataFrame,
    rf_shap_dir: Path,
    arrangement: str,
    top_k: int,
    output_dir: Path,
) -> Optional[pd.DataFrame]:
    path = find_species_consensus_file(rf_shap_dir, arrangement)
    if path is None:
        print(f"[WARN] RF SHAP species consensus not found for {arrangement}; skipping RF comparison.")
        return None

    rf = pd.read_csv(path)
    if "Species_key" not in rf.columns:
        if "Species" in rf.columns:
            rf["Species_key"] = rf["Species"].map(canonical_key)
        elif "Species_clean" in rf.columns:
            rf["Species_key"] = rf["Species_clean"].map(canonical_key)
        else:
            print(f"[WARN] RF file {path} has no Species/Species_key column; skipping.")
            return None

    if "SpeciesRank" not in rf.columns:
        rf["SpeciesRank"] = np.arange(1, len(rf) + 1)

    rf_top = set(rf.head(top_k)["Species_key"])
    rf_any = set(rf["Species_key"])

    out = lasso_species.copy()
    out["In_RF_SHAP_TopK"] = out["Species_key"].isin(rf_top)
    out["In_RF_SHAP_Any"] = out["Species_key"].isin(rf_any)

    # Add RF rank if available.
    rank_map = rf.drop_duplicates("Species_key").set_index("Species_key")["SpeciesRank"].to_dict()
    out["RF_SHAP_SpeciesRank"] = out["Species_key"].map(rank_map)

    out.head(top_k).to_csv(output_dir / arrangement / f"{arrangement}_LASSO_top{top_k}_vs_RF_SHAP.csv", index=False)
    return out


def find_existing_file(candidates: Iterable[Path], required: bool) -> Optional[Path]:
    for p in candidates:
        if p.exists():
            return p
    if required:
        raise FileNotFoundError("Could not find any of:\n" + "\n".join(map(str, candidates)))
    return None


def phase3_file_candidates(phase3_dir: Path, subdir: str, filename: str) -> List[Path]:
    return [
        phase3_dir / subdir / filename,
        phase3_dir / filename,
    ] + sorted(phase3_dir.glob(f"**/{filename}"))


def load_phase3_tables(phase3_dir: Path) -> Dict[str, pd.DataFrame]:
    files = {
        "3.1": ("3.1", "step3_1_blocked_wilcoxon_significant_results.csv"),
        "3.2": ("3.2", "step3_2_cmh_significant_results.csv"),
        "3.3": ("3.3", "step3_3_meta_analysis_significant_results.csv"),
        "3.5": ("3.5", "step3_5_biomarker_categories.csv"),
    }
    out = {}
    for step, (subdir, filename) in files.items():
        required = step in {"3.1", "3.2", "3.3"}
        path = find_existing_file(phase3_file_candidates(phase3_dir, subdir, filename), required=required)
        if path is None:
            continue
        df = pd.read_csv(path)
        if "species" not in df.columns and "Species" in df.columns:
            df = df.rename(columns={"Species": "species"})
        if "species" not in df.columns:
            raise ValueError(f"{path} does not contain a species/Species column")
        df = df.copy()
        df["Species_clean"] = df["species"].map(clean_species_name)
        df["Species_key"] = df["species"].map(canonical_key)
        df["Phase3Step"] = step
        out[step] = df
    return out


def apply_filters(df: pd.DataFrame, filters: Dict[str, Set[str]]) -> pd.DataFrame:
    out = df.copy()
    for col, allowed in filters.items():
        if col not in out.columns:
            return out.iloc[0:0].copy()
        out = out[out[col].astype(str).isin(allowed)]
    return out.copy()


def build_phase3_sets(phase3: Dict[str, pd.DataFrame]) -> Dict[str, Dict[str, Set[str]]]:
    sets: Dict[str, Dict[str, Set[str]]] = {}
    for arr in ["All", "OnlyMale", "OnlyFemale"]:
        sets[arr] = {}
        for step in ["3.1", "3.2", "3.3"]:
            df = phase3.get(step, pd.DataFrame())
            filtered = apply_filters(df, PHASE3_FILTERS[arr][step]) if not df.empty else df
            sets[arr][step] = set(filtered["Species_key"].dropna()) if not filtered.empty else set()
        sets[arr]["formal_union"] = sets[arr]["3.1"] | sets[arr]["3.2"] | sets[arr]["3.3"]
        sets[arr]["supported_2plus"] = {
            s for s in sets[arr]["formal_union"]
            if sum(s in sets[arr][step] for step in ["3.1", "3.2", "3.3"]) >= 2
        }
    return sets


def load_phase35_annotations(phase3: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    df = phase3.get("3.5", pd.DataFrame()).copy()
    if df.empty:
        return df
    keep = [
        "Species_key", "Species_clean", "categories", "is_early", "is_late",
        "is_progression", "is_ordered_increasing", "is_ordered_decreasing",
        "best_trend_direction", "best_trend_q",
    ]
    keep = [c for c in keep if c in df.columns]
    return df[keep].drop_duplicates("Species_key")


def compare_with_phase3(
    lasso_species: pd.DataFrame,
    phase3_sets: Dict[str, Dict[str, Set[str]]],
    phase35: pd.DataFrame,
    arrangement: str,
    top_k: int,
    output_dir: Path,
) -> pd.DataFrame:
    out = lasso_species.copy()
    out["In_Phase3_1_BlockedWilcoxon_rCLR"] = out["Species_key"].isin(phase3_sets[arrangement]["3.1"])
    out["In_Phase3_2_CMH_PA"] = out["Species_key"].isin(phase3_sets[arrangement]["3.2"])
    out["In_Phase3_3_MetaAnalysis_rCLR"] = out["Species_key"].isin(phase3_sets[arrangement]["3.3"])
    out["In_Any_Formal_Phase3"] = out["Species_key"].isin(phase3_sets[arrangement]["formal_union"])
    out["In_2plus_Formal_Phase3"] = out["Species_key"].isin(phase3_sets[arrangement]["supported_2plus"])
    out["N_Formal_Phase3_Supports"] = (
        out["In_Phase3_1_BlockedWilcoxon_rCLR"].astype(int)
        + out["In_Phase3_2_CMH_PA"].astype(int)
        + out["In_Phase3_3_MetaAnalysis_rCLR"].astype(int)
    )
    if not phase35.empty:
        out = out.merge(phase35, on="Species_key", how="left", suffixes=("", "_phase35"))

    out.head(top_k).to_csv(output_dir / arrangement / f"{arrangement}_LASSO_top{top_k}_vs_Phase3.csv", index=False)
    out.to_csv(output_dir / arrangement / f"{arrangement}_LASSO_all_species_vs_Phase3.csv", index=False)
    return out


def plot_top_coefficients(consensus: pd.DataFrame, output_path: Path, arrangement: str, top_n: int) -> None:
    df = consensus.head(top_n).iloc[::-1].copy()
    values = df["MeanNonZeroCoefficient"].where(df["MeanNonZeroCoefficient"] != 0, df["MeanCoefficient"])
    plt.figure(figsize=(10, max(5, 0.36 * len(df))))
    plt.barh(df["Feature"], values)
    plt.axvline(0, linewidth=1)
    plt.xlabel("Mean non-zero coefficient")
    plt.ylabel("Feature")
    plt.title(f"{arrangement} LASSO: top {top_n} consensus coefficients")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_feature_type_summary(summary: pd.DataFrame, output_path: Path, arrangement: str, top_n: int) -> None:
    plt.figure(figsize=(7, 4.5))
    plt.bar(summary["Feature_Type"], summary["Percent_Coefficient_TopN"])
    plt.ylabel("% of top-feature coefficient magnitude")
    plt.xlabel("Feature type")
    plt.title(f"{arrangement} LASSO: PA/rCLR split in top {top_n}")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_overlap_stats(stats: pd.DataFrame, output_path: Path) -> None:
    row = stats.iloc[0]
    labels = ["All∩OnlyMale", "All∩OnlyFemale", "OnlyMale∩OnlyFemale", "All three"]
    values = [
        row["N_Shared_All_OnlyMale"],
        row["N_Shared_All_OnlyFemale"],
        row["N_Shared_OnlyMale_OnlyFemale"],
        row["N_Shared_All_Three"],
    ]
    plt.figure(figsize=(8, 4.5))
    plt.bar(labels, values)
    plt.ylabel("Shared top species")
    plt.title(f"LASSO top-{int(row['Top_N'])} species overlap across arrangements")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_performance(perf: pd.DataFrame, output_path: Path) -> None:
    order = ["All", "OnlyMale", "OnlyFemale"]
    data = [perf.loc[perf["Arrangement"] == arr, "AUC"].dropna().values for arr in order]
    plt.figure(figsize=(8, 4.8))
    plt.boxplot(data, labels=order, showmeans=True)
    plt.ylabel("LODO AUC")
    plt.title("LASSO Logistic Regression performance across LODO folds")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def df_to_md(df: pd.DataFrame, max_rows: Optional[int] = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    try:
        return df.to_markdown(index=False)
    except Exception:
        return df.to_string(index=False)


def write_report(
    output_dir: Path,
    performance_summary: pd.DataFrame,
    arrangement_consensus: Dict[str, pd.DataFrame],
    species_consensus: Dict[str, pd.DataFrame],
    type_summaries: Dict[str, pd.DataFrame],
    overlap_stats: pd.DataFrame,
    rf_comparisons: Dict[str, pd.DataFrame],
    phase3_comparisons: Dict[str, pd.DataFrame],
    top_k: int,
    tune_c: bool,
) -> None:
    lines: List[str] = []
    lines.append("# Step 5.3 — LASSO Logistic Regression Interpretability\n")
    lines.append("\n## Purpose\n")
    lines.append(
        "This analysis explores interpretability using an L1-penalized Logistic Regression model. "
        "The L1 penalty induces sparsity by shrinking some feature coefficients exactly to zero. "
        "Therefore, selected non-zero features form a compact and directly interpretable microbial signature.\n"
    )
    lines.append("\n## Arrangements and model\n")
    lines.append("- Arrangements: **All**, **OnlyMale**, **OnlyFemale**.\n")
    lines.append("- Model: **LASSO Logistic Regression**.\n")
    lines.append("- Validation: same LODO folds as Phase 4.\n")
    lines.append(f"- C tuning inside train folds: **{'yes' if tune_c else 'no'}**.\n")
    lines.append("\n## Performance and sparsity summary\n")
    lines.append(df_to_md(performance_summary))
    lines.append("\n")

    for arr in ["All", "OnlyMale", "OnlyFemale"]:
        cons = arrangement_consensus[arr]
        sp = species_consensus[arr]
        lines.append(f"\n## {arr}\n")
        lines.append("Top consensus coefficients:\n")
        cols = [
            "ConsensusRank", "Feature", "Feature_Type", "Species_clean", "ConsensusDirection",
            "SelectionFrequency", "N_Folds_NonZero", "MeanNonZeroCoefficient",
            "MeanAbsNonZeroCoefficient", "ConsensusScore",
        ]
        lines.append(df_to_md(cons.head(top_k)[cols]))
        lines.append("\n\nFeature-type split among top coefficients:\n")
        lines.append(df_to_md(type_summaries[arr]))
        lines.append("\n\nTop species-level LASSO consensus:\n")
        sp_cols = [
            "SpeciesRank", "Species_clean", "Dominant_Signal", "Species_MaxSelectionFrequency",
            "Species_SumMeanAbsNonZeroCoef", "Any_CRC_Associated", "Any_Control_Associated",
        ]
        lines.append(df_to_md(sp.head(top_k)[sp_cols]))
        lines.append("\n")

        if arr in rf_comparisons:
            top = rf_comparisons[arr].head(top_k)
            n_rf = int(top["In_RF_SHAP_TopK"].sum())
            lines.append(f"\nLASSO vs Random Forest SHAP: {n_rf}/{len(top)} top LASSO species are also top RF SHAP species.\n")

        if arr in phase3_comparisons:
            top = phase3_comparisons[arr].head(top_k)
            n_any = int(top["In_Any_Formal_Phase3"].sum())
            n_meta = int(top["In_Phase3_3_MetaAnalysis_rCLR"].sum())
            lines.append(
                f"LASSO vs Phase 3: {n_any}/{len(top)} top LASSO species overlap with formal Phase 3 results; "
                f"{n_meta}/{len(top)} overlap with random-effects meta-analysis.\n"
            )

    lines.append("\n## Arrangement overlap\n")
    lines.append(df_to_md(overlap_stats))
    lines.append("\n")

    (output_dir / "step5_3_lasso_interpretability_report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step 5.3 LASSO Logistic Regression interpretability")
    p.add_argument("--input-dir", default="OUTPUTS/phase_4/4.1/4.1_prep", help="Phase 4 LODO prep directory")
    p.add_argument("--output-dir", default="OUTPUTS/phase_5/5.3_lasso_interpretability", help="Output directory")
    p.add_argument("--rf-shap-dir", default=None, help="Optional Step 5.1 Random Forest SHAP directory")
    p.add_argument("--phase3-dir", default=None, help="Optional Phase 3 output directory")
    p.add_argument("--top-k", type=int, default=20, help="Top K features/species for reports and overlap")
    p.add_argument("--coef-threshold", type=float, default=1e-8, help="Absolute coefficient threshold for non-zero selection")
    p.add_argument("--tune-c", action="store_true", help="Use LogisticRegressionCV to tune C inside each train fold")
    p.add_argument("--c-value", type=float, default=1.0, help="Fixed C if --tune-c is not used")
    p.add_argument("--cv-folds", type=int, default=5, help="Inner CV folds for C tuning")
    p.add_argument(
        "--c-grid",
        default="1e-3,3e-3,1e-2,3e-2,1e-1,3e-1,1,3,10,30,100",
        help="Comma-separated C grid for --tune-c",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    figures_dir = output_dir / "figures"
    safe_mkdir(output_dir)
    safe_mkdir(figures_dir)

    c_grid = np.array([float(x) for x in args.c_grid.split(",") if x.strip()])

    arrangement_consensus: Dict[str, pd.DataFrame] = {}
    species_cons: Dict[str, pd.DataFrame] = {}
    type_summaries: Dict[str, pd.DataFrame] = {}
    performance_tables: List[pd.DataFrame] = []

    for arr in ["All", "OnlyMale", "OnlyFemale"]:
        print(f"\n{'='*72}\nRunning LASSO interpretability: {arr}\n{'='*72}")
        all_coefs, perf = run_arrangement_lasso(
            input_dir=input_dir,
            output_dir=output_dir,
            arrangement=arr,
            top_k=args.top_k,
            coef_threshold=args.coef_threshold,
            tune_c=args.tune_c,
            c_value=args.c_value,
            cv_folds=args.cv_folds,
            c_grid=c_grid,
        )

        cons = consensus_coefficients(all_coefs, top_k=args.top_k, coef_threshold=args.coef_threshold)
        sp = species_level_consensus(cons)
        ftype = feature_type_summary(cons, top_n=args.top_k)

        arrangement_consensus[arr] = cons
        species_cons[arr] = sp
        type_summaries[arr] = ftype
        performance_tables.append(perf)

        cons.to_csv(output_dir / arr / f"{arr}_LASSO_consensus_coefficients.csv", index=False)
        sp.to_csv(output_dir / arr / f"{arr}_LASSO_species_level_consensus.csv", index=False)
        ftype.to_csv(output_dir / arr / f"{arr}_LASSO_feature_type_summary_top{args.top_k}.csv", index=False)

        plot_top_coefficients(
            cons,
            figures_dir / f"{arr}_LASSO_top{args.top_k}_consensus_coefficients.png",
            arrangement=arr,
            top_n=args.top_k,
        )
        plot_feature_type_summary(
            ftype,
            figures_dir / f"{arr}_LASSO_feature_type_top{args.top_k}.png",
            arrangement=arr,
            top_n=args.top_k,
        )

    perf_all = pd.concat(performance_tables, ignore_index=True)
    perf_all.to_csv(output_dir / "lasso_fold_performance_all_arrangements.csv", index=False)
    performance_summary = (
        perf_all.groupby("Arrangement", as_index=False)
        .agg(
            N_folds=("Study", "nunique"),
            Median_AUC=("AUC", "median"),
            IQR_AUC_low=("AUC", lambda x: x.quantile(0.25)),
            IQR_AUC_high=("AUC", lambda x: x.quantile(0.75)),
            Median_AUPRC=("AUPRC", "median"),
            Median_F1=("F1", "median"),
            Median_N_nonzero_features=("N_nonzero_features", "median"),
            IQR_N_nonzero_low=("N_nonzero_features", lambda x: x.quantile(0.25)),
            IQR_N_nonzero_high=("N_nonzero_features", lambda x: x.quantile(0.75)),
            Median_Percent_nonzero_features=("Percent_nonzero_features", "median"),
            Median_Selected_C=("Selected_C", "median"),
        )
    )
    performance_summary.to_csv(output_dir / "lasso_model_performance_summary.csv", index=False)
    plot_performance(perf_all, figures_dir / "lasso_performance_auc.png")

    overlap, overlap_stats = arrangement_overlap(species_cons, top_n=args.top_k, output_dir=output_dir)
    plot_overlap_stats(overlap_stats, figures_dir / f"lasso_arrangement_top{args.top_k}_overlap.png")

    rf_comparisons: Dict[str, pd.DataFrame] = {}
    if args.rf_shap_dir:
        rf_dir = Path(args.rf_shap_dir)
        if rf_dir.exists():
            for arr in ["All", "OnlyMale", "OnlyFemale"]:
                comp = compare_with_rf_shap(
                    lasso_species=species_cons[arr],
                    rf_shap_dir=rf_dir,
                    arrangement=arr,
                    top_k=args.top_k,
                    output_dir=output_dir,
                )
                if comp is not None:
                    rf_comparisons[arr] = comp
        else:
            print(f"[WARN] --rf-shap-dir does not exist: {rf_dir}")

    phase3_comparisons: Dict[str, pd.DataFrame] = {}
    if args.phase3_dir:
        phase3_dir = Path(args.phase3_dir)
        if phase3_dir.exists():
            phase3 = load_phase3_tables(phase3_dir)
            phase3_sets = build_phase3_sets(phase3)
            phase35 = load_phase35_annotations(phase3)
            for arr in ["All", "OnlyMale", "OnlyFemale"]:
                comp = compare_with_phase3(
                    lasso_species=species_cons[arr],
                    phase3_sets=phase3_sets,
                    phase35=phase35,
                    arrangement=arr,
                    top_k=args.top_k,
                    output_dir=output_dir,
                )
                phase3_comparisons[arr] = comp
        else:
            print(f"[WARN] --phase3-dir does not exist: {phase3_dir}")

    write_report(
        output_dir=output_dir,
        performance_summary=performance_summary,
        arrangement_consensus=arrangement_consensus,
        species_consensus=species_cons,
        type_summaries=type_summaries,
        overlap_stats=overlap_stats,
        rf_comparisons=rf_comparisons,
        phase3_comparisons=phase3_comparisons,
        top_k=args.top_k,
        tune_c=args.tune_c,
    )

    print("\nStep 5.3 complete.")
    print(f"Outputs written to: {output_dir}")
    print("Main report:", output_dir / "step5_3_lasso_interpretability_report.md")
    print("Performance summary:", output_dir / "lasso_model_performance_summary.csv")
    print("Figures:", figures_dir)


if __name__ == "__main__":
    main()
