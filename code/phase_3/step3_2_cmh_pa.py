#!/usr/bin/env python3
"""
STEP 3.2 | PA features differential prevalence analysis
Binary PA features — Cochran-Mantel-Haenszel tests blocked by study

Run from:
    project/code/

Expected project structure:
    project/
    ├── code/
    │   └── step32_cmh_pa.py
    ├── data/
    │   └── crc_study_final_data/
    │       └── species_level/
    │           ├── X_explore.csv
    │           └── metadata.csv
    └── OUTPUTS/
        └── step_3_2/

Outputs:
    OUTPUTS/step_3_2/step32_cmh_all_results.csv
    OUTPUTS/step_3_2/step32_cmh_significant_results.csv
    OUTPUTS/step_3_2/step32_cmh_summary.png
    OUTPUTS/step_3_2/step32_final_report.txt
"""

from __future__ import annotations

import math
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=RuntimeWarning)

# =============================================================================
# Paths: run this script from the code/ directory
# =============================================================================
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

DATA_DIR = (PROJECT_ROOT / "data" / "crc_study_final_data" / "species_level").resolve()
OUTPUT_DIR = (PROJECT_ROOT / "OUTPUTS" / "phase_3" / "3.2").resolve()

X_PATH = (DATA_DIR / "X_explore.csv").resolve()
METADATA_PATH = (DATA_DIR / "metadata.csv").resolve()

# =============================================================================
# Settings
# =============================================================================
SAMPLE_COL = "Sample"
STUDY_COL = "Study"
CONDITION_COL = "Condition"
SEX_COL = "Gender"
ALPHA = 0.05
MIN_STRATA = 2

# Comparisons.
# Each comparison is defined as label, grouping column, positive/first group, negative/second group,
# and optional filters.
COMPARISONS = [
    {
        "label": "CRC_vs_Control__all",
        "group_col": CONDITION_COL,
        "group1": "CRC",
        "group2": "Control",
        "filters": {},
        "description": "CRC vs Control, all samples",
    },
    {
        "label": "CRC_vs_Control__Male",
        "group_col": CONDITION_COL,
        "group1": "CRC",
        "group2": "Control",
        "filters": {SEX_COL: "Male"},
        "description": "CRC vs Control, males only",
    },
    {
        "label": "CRC_vs_Control__Female",
        "group_col": CONDITION_COL,
        "group1": "CRC",
        "group2": "Control",
        "filters": {SEX_COL: "Female"},
        "description": "CRC vs Control, females only",
    },
    {
        "label": "Adenoma_vs_Control__all",
        "group_col": CONDITION_COL,
        "group1": "Adenoma",
        "group2": "Control",
        "filters": {},
        "description": "Adenoma vs Control, all samples",
    },
    {
        "label": "CRC_vs_Adenoma__all",
        "group_col": CONDITION_COL,
        "group1": "CRC",
        "group2": "Adenoma",
        "filters": {},
        "description": "CRC vs Adenoma, all samples",
    },
    # Sex comparisons within each condition. These correspond to the Step 3.1
    # Male vs Female comparison stratified by condition, but using PA features.
    {
        "label": "Male_vs_Female__Control",
        "group_col": SEX_COL,
        "group1": "Male",
        "group2": "Female",
        "filters": {CONDITION_COL: "Control"},
        "description": "Male vs Female within Control samples",
    },
    {
        "label": "Male_vs_Female__Adenoma",
        "group_col": SEX_COL,
        "group1": "Male",
        "group2": "Female",
        "filters": {CONDITION_COL: "Adenoma"},
        "description": "Male vs Female within Adenoma samples",
    },
    {
        "label": "Male_vs_Female__CRC",
        "group_col": SEX_COL,
        "group1": "Male",
        "group2": "Female",
        "filters": {CONDITION_COL: "CRC"},
        "description": "Male vs Female within CRC samples",
    },
]


def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_inputs() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[str]]:
    if not X_PATH.exists():
        raise FileNotFoundError(f"X_explore.csv not found at: {X_PATH}")
    if not METADATA_PATH.exists():
        raise FileNotFoundError(f"metadata.csv not found at: {METADATA_PATH}")

    x = pd.read_csv(X_PATH)
    metadata = pd.read_csv(METADATA_PATH)

    required_x = {SAMPLE_COL}
    required_meta = {SAMPLE_COL, STUDY_COL, CONDITION_COL, SEX_COL}
    missing_x = required_x - set(x.columns)
    missing_meta = required_meta - set(metadata.columns)
    if missing_x:
        raise ValueError(f"X_explore.csv is missing required columns: {sorted(missing_x)}")
    if missing_meta:
        raise ValueError(f"metadata.csv is missing required columns: {sorted(missing_meta)}")

    pa_cols = [c for c in x.columns if c.startswith("PA_")]
    if len(pa_cols) == 0:
        raise ValueError("No PA_ columns found in X_explore.csv. Step 3.2 requires binary PA features.")

    # Keep only PA columns and metadata columns needed for this analysis.
    x_small = x[[SAMPLE_COL] + pa_cols].copy()
    df = metadata[[SAMPLE_COL, STUDY_COL, CONDITION_COL, SEX_COL]].merge(
        x_small, on=SAMPLE_COL, how="inner"
    )

    # Drop samples with missing essential metadata.
    df = df.dropna(subset=[STUDY_COL, CONDITION_COL, SEX_COL]).copy()

    # Ensure PA features are binary integer 0/1.
    for col in pa_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
        df[col] = (df[col] > 0).astype(int)

    return x, metadata, df, pa_cols


def species_name_from_pa_col(col: str) -> str:
    return col.replace("PA_", "", 1)


def apply_filters(df: pd.DataFrame, filters: Dict[str, str]) -> pd.DataFrame:
    out = df
    for col, value in filters.items():
        out = out[out[col] == value]
    return out.copy()


def build_2x2_tables_by_study(
    df: pd.DataFrame,
    feature_col: str,
    group_col: str,
    group1: str,
    group2: str,
) -> Tuple[List[np.ndarray], Dict[str, int]]:
    """
    Build one 2x2 table per study:

        [[group1 present, group1 absent],
         [group2 present, group2 absent]]

    Only studies containing both groups are used.
    """
    tables: List[np.ndarray] = []
    stats = {
        "n_strata_total": int(df[STUDY_COL].nunique()),
        "n_strata_used": 0,
        "n_strata_skipped_missing_group": 0,
        "n_strata_skipped_empty": 0,
    }

    for study, sub in df.groupby(STUDY_COL):
        sub = sub[sub[group_col].isin([group1, group2])]
        if sub[group_col].nunique() < 2:
            stats["n_strata_skipped_missing_group"] += 1
            continue

        g1 = sub[sub[group_col] == group1]
        g2 = sub[sub[group_col] == group2]
        if len(g1) == 0 or len(g2) == 0:
            stats["n_strata_skipped_empty"] += 1
            continue

        a = int(g1[feature_col].sum())
        b = int(len(g1) - a)
        c = int(g2[feature_col].sum())
        d = int(len(g2) - c)
        table = np.array([[a, b], [c, d]], dtype=float)
        tables.append(table)
        stats["n_strata_used"] += 1

    return tables, stats


def cmh_test(
    tables: List[np.ndarray],
) -> Tuple[float, float, float, float, float]:
    """
    Return CMH statistic, p-value, pooled common odds ratio,
    pooled common odds ratio lower CI, pooled common odds ratio upper CI.

    shift_zeros=True applies a small continuity correction to strata with zero cells,
    which stabilizes common odds ratio estimation in sparse binary microbiome data.
    """
    if len(tables) < MIN_STRATA:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    st = StratifiedTable(tables, shift_zeros=True)
    test = st.test_null_odds()
    stat = float(test.statistic)
    pval = float(test.pvalue)

    try:
        common_or = float(st.oddsratio_pooled)
    except Exception:
        common_or = np.nan

    try:
        ci_low, ci_high = st.oddsratio_pooled_confint()
        ci_low = float(ci_low)
        ci_high = float(ci_high)
    except Exception:
        ci_low, ci_high = np.nan, np.nan

    return stat, pval, common_or, ci_low, ci_high


def run_one_comparison(df: pd.DataFrame, pa_cols: List[str], comp: Dict[str, object]) -> pd.DataFrame:
    label = str(comp["label"])
    group_col = str(comp["group_col"])
    group1 = str(comp["group1"])
    group2 = str(comp["group2"])
    filters = dict(comp.get("filters", {}))
    description = str(comp.get("description", label))

    sub = apply_filters(df, filters)
    sub = sub[sub[group_col].isin([group1, group2])].copy()

    rows = []
    n_group1 = int((sub[group_col] == group1).sum())
    n_group2 = int((sub[group_col] == group2).sum())

    for feature in pa_cols:
        tables, table_stats = build_2x2_tables_by_study(sub, feature, group_col, group1, group2)

        if n_group1 > 0:
            prev1 = float(sub.loc[sub[group_col] == group1, feature].mean())
        else:
            prev1 = np.nan
        if n_group2 > 0:
            prev2 = float(sub.loc[sub[group_col] == group2, feature].mean())
        else:
            prev2 = np.nan

        stat, pval, common_or, ci_low, ci_high = cmh_test(tables)

        rows.append(
            {
                "comparison": label,
                "description": description,
                "feature": feature,
                "species": species_name_from_pa_col(feature),
                "group_col": group_col,
                "group1": group1,
                "group2": group2,
                "filters": "; ".join([f"{k}={v}" for k, v in filters.items()]) if filters else "none",
                "n_group1": n_group1,
                "n_group2": n_group2,
                "prevalence_group1": prev1,
                "prevalence_group2": prev2,
                "prevalence_diff_group1_minus_group2": prev1 - prev2 if pd.notna(prev1) and pd.notna(prev2) else np.nan,
                "n_strata_total": table_stats["n_strata_total"],
                "n_strata_used": table_stats["n_strata_used"],
                "cmh_statistic": stat,
                "p_value": pval,
                "common_odds_ratio": common_or,
                "common_odds_ratio_ci_low": ci_low,
                "common_odds_ratio_ci_high": ci_high,
            }
        )

    res = pd.DataFrame(rows)

    # FDR correction is done within each comparison across all PA species.
    valid = res["p_value"].notna()
    res["q_value_bh"] = np.nan
    res["significant_fdr_0_05"] = False
    if valid.sum() > 0:
        _, qvals, _, _ = multipletests(res.loc[valid, "p_value"].values, alpha=ALPHA, method="fdr_bh")
        res.loc[valid, "q_value_bh"] = qvals
        res.loc[valid, "significant_fdr_0_05"] = res.loc[valid, "q_value_bh"] < ALPHA

    # Useful sort: significant first, then q-value, then absolute prevalence difference.
    res = res.sort_values(
        by=["comparison", "q_value_bh", "p_value", "prevalence_diff_group1_minus_group2"],
        ascending=[True, True, True, False],
        na_position="last",
    ).reset_index(drop=True)
    return res


def run_all_tests(df: pd.DataFrame, pa_cols: List[str]) -> pd.DataFrame:
    all_results = []
    for comp in COMPARISONS:
        all_results.append(run_one_comparison(df, pa_cols, comp))
    return pd.concat(all_results, ignore_index=True)


def make_summary_figure(results: pd.DataFrame, output_path: Path) -> None:
    """Create one compact figure summarizing significant counts and effect directions."""
    summary = (
        results.groupby("comparison")
        .agg(
            n_tested=("p_value", lambda s: int(s.notna().sum())),
            n_significant=("significant_fdr_0_05", "sum"),
            median_abs_prev_diff=("prevalence_diff_group1_minus_group2", lambda s: float(np.nanmedian(np.abs(s)))),
        )
        .reset_index()
    )
    summary["n_significant"] = summary["n_significant"].astype(int)
    summary = summary.sort_values("n_significant", ascending=True)

    fig, axes = plt.subplots(1, 2, figsize=(15, 7))

    axes[0].barh(summary["comparison"], summary["n_significant"])
    axes[0].set_xlabel("Number of significant PA features, FDR < 0.05")
    axes[0].set_ylabel("")
    axes[0].set_title("CMH significant species by comparison")

    axes[1].barh(summary["comparison"], summary["median_abs_prev_diff"])
    axes[1].set_xlabel("Median absolute prevalence difference")
    axes[1].set_ylabel("")
    axes[1].set_title("Typical absolute prevalence shift")

    fig.suptitle("Step 3.2 — PA features: study-blocked CMH tests", fontsize=14)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def format_float(x: object, digits: int = 4) -> str:
    try:
        if pd.isna(x):
            return "NA"
        return f"{float(x):.{digits}g}"
    except Exception:
        return str(x)


def generate_report(
    df: pd.DataFrame,
    pa_cols: List[str],
    results: pd.DataFrame,
    report_path: Path,
) -> None:
    lines: List[str] = []
    lines.append("STEP 3.2 | PA FEATURES — COCHRAN-MANTEL-HAENSZEL TESTS")
    lines.append("=" * 78)
    lines.append("")

    lines.append("INPUTS")
    lines.append("-" * 78)
    lines.append(f"X_explore path: {X_PATH}")
    lines.append(f"Metadata path:  {METADATA_PATH}")
    lines.append(f"Merged samples used: {len(df)}")
    lines.append(f"PA features tested: {len(pa_cols)}")
    lines.append(f"Studies: {df[STUDY_COL].nunique()}")
    lines.append("")

    lines.append("Condition counts:")
    for k, v in df[CONDITION_COL].value_counts().items():
        lines.append(f"  {k}: {v}")
    lines.append("")

    lines.append("Sex counts:")
    for k, v in df[SEX_COL].value_counts().items():
        lines.append(f"  {k}: {v}")
    lines.append("")

    lines.append("METHOD")
    lines.append("-" * 78)
    lines.append("For each PA species feature, a Cochran-Mantel-Haenszel test was run")
    lines.append("using Study as the blocking/stratification variable. Each study contributes")
    lines.append("a 2x2 table of presence/absence by group, when both groups are present.")
    lines.append("FDR correction was performed separately within each comparison using")
    lines.append("Benjamini-Hochberg correction.")
    lines.append("")
    lines.append("For sparse 2x2 tables, a continuity correction is used for stable common")
    lines.append("odds-ratio estimation. The hypothesis test remains the stratified association")
    lines.append("test analogous to R's mantelhaen.test(presence, condition, study).")
    lines.append("")

    lines.append("COMPARISON SUMMARY")
    lines.append("-" * 78)
    summary = (
        results.groupby("comparison")
        .agg(
            description=("description", "first"),
            n_tested=("p_value", lambda s: int(s.notna().sum())),
            n_significant=("significant_fdr_0_05", "sum"),
            min_q=("q_value_bh", "min"),
            median_abs_prev_diff=("prevalence_diff_group1_minus_group2", lambda s: float(np.nanmedian(np.abs(s)))),
            median_strata_used=("n_strata_used", "median"),
        )
        .reset_index()
        .sort_values("comparison")
    )

    for _, row in summary.iterrows():
        lines.append(f"{row['comparison']}")
        lines.append(f"  Description: {row['description']}")
        lines.append(f"  Tested species: {int(row['n_tested'])}")
        lines.append(f"  Significant species, FDR < 0.05: {int(row['n_significant'])}")
        lines.append(f"  Minimum q-value: {format_float(row['min_q'], 4)}")
        lines.append(f"  Median |prevalence difference|: {format_float(row['median_abs_prev_diff'], 4)}")
        lines.append(f"  Median number of study strata used: {format_float(row['median_strata_used'], 4)}")
        lines.append("")

    lines.append("TOP SPECIES PER COMPARISON")
    lines.append("-" * 78)
    for comp in sorted(results["comparison"].unique()):
        sub = results[results["comparison"] == comp].copy()
        sub = sub[sub["q_value_bh"].notna()].sort_values(["q_value_bh", "p_value"]).head(10)
        lines.append(comp)
        if sub.empty:
            lines.append("  No valid tests.")
        else:
            for _, r in sub.iterrows():
                direction = "higher in " + str(r["group1"]) if r["prevalence_diff_group1_minus_group2"] > 0 else "higher in " + str(r["group2"])
                lines.append(
                    "  - "
                    + str(r["species"])
                    + f" | q={format_float(r['q_value_bh'], 4)}"
                    + f", p={format_float(r['p_value'], 4)}"
                    + f", common OR={format_float(r['common_odds_ratio'], 4)}"
                    + f", prev_diff={format_float(r['prevalence_diff_group1_minus_group2'], 4)}"
                    + f" ({direction})"
                    + f", strata={int(r['n_strata_used'])}"
                )
        lines.append("")

    lines.append("INTERPRETATION GUIDE")
    lines.append("-" * 78)
    lines.append("- common OR > 1: the species is more likely to be present in group1 than group2, after blocking by study.")
    lines.append("- common OR < 1: the species is more likely to be present in group2 than group1, after blocking by study.")
    lines.append("- q_value_bh < 0.05: significant after FDR correction within that comparison.")
    lines.append("- These tests examine presence/absence only. Abundance shifts among present samples are tested separately in Step 3.1 using rCLR features.")
    lines.append("")

    lines.append("FILES WRITTEN")
    lines.append("-" * 78)
    lines.append(f"All results CSV:         {OUTPUT_DIR / 'step3_2_cmh_all_results.csv'}")
    lines.append(f"Significant results CSV: {OUTPUT_DIR / 'step3_2_cmh_significant_results.csv'}")
    lines.append(f"Summary figure:          {OUTPUT_DIR / 'step3_2_cmh_summary.png'}")
    lines.append(f"Final report:            {report_path}")
    lines.append("")

    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    safe_mkdir(OUTPUT_DIR)
    _, _, df, pa_cols = load_inputs()

    results = run_all_tests(df, pa_cols)
    sig = results[results["significant_fdr_0_05"]].copy()

    all_path = OUTPUT_DIR / "step3_2_cmh_all_results.csv"
    sig_path = OUTPUT_DIR / "step3_2_cmh_significant_results.csv"
    fig_path = OUTPUT_DIR / "step3_2_cmh_summary.png"
    report_path = OUTPUT_DIR / "step3_2_final_report.txt"

    results.to_csv(all_path, index=False)
    sig.to_csv(sig_path, index=False)
    make_summary_figure(results, fig_path)
    generate_report(df, pa_cols, results, report_path)

    print("Step 3.2 complete.")
    print(f"All results:         {all_path}")
    print(f"Significant results: {sig_path}")
    print(f"Summary figure:      {fig_path}")
    print(f"Final report:        {report_path}")


if __name__ == "__main__":
    main()
