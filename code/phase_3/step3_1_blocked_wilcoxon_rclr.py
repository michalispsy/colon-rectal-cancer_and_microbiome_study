#!/usr/bin/env python3
"""
STEP 3.1 — Differential Abundance Analysis on rCLR features
Blocked Wilcoxon / van Elteren test in Python

Run from the project/code directory:
    python step31_blocked_wilcoxon_rclr.py

Expected project structure:
    project/
    ├── code/
    │   └── step31_blocked_wilcoxon_rclr.py
    ├── data/
    │   └── crc_study_final_data/
    │       └── species_level/
    │           ├── X_explore.csv
    │           └── metadata.csv
    └── OUTPUTS/
        └── step_3_1/

Outputs:
    step31_blocked_wilcoxon_all_results.csv
    step31_blocked_wilcoxon_significant_results.csv
    step31_blocked_wilcoxon_summary.png
    step31_final_report.txt

Notes:
- Uses only rCLR_ features from X_explore.csv.
- Blocks/stratifies by Study.
- Implements a van Elteren test, the stratified Wilcoxon rank-sum analogue.
- FDR correction is Benjamini-Hochberg, applied separately per comparison family.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import rankdata, norm
from statsmodels.stats.multitest import multipletests


# =============================================================================
# Paths — run this script from code/, but paths are resolved from the script file
# =============================================================================
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

DATA_DIR = (PROJECT_ROOT / "data" / "crc_study_final_data" / "species_level").resolve()
OUTPUT_DIR = (PROJECT_ROOT / "OUTPUTS" / "phase_3" / "3.1" ).resolve()

X_PATH = (DATA_DIR / "X_explore.csv").resolve()
METADATA_PATH = (DATA_DIR / "metadata.csv").resolve()


# =============================================================================
# User-adjustable parameters
# =============================================================================
SAMPLE_COL = "Sample"
STUDY_COL = "Study"
CONDITION_COL = "Condition"
SEX_COL = "Gender"

RCLR_PREFIX = "rCLR_"
FDR_ALPHA = 0.05
MIN_PER_GROUP_PER_STUDY = 2     # study contributes only if both groups have >= this many samples
MIN_CONTRIBUTING_STUDIES = 2    # feature/comparison tested only if >= this many valid study blocks

RANDOM_STATE = 42


# =============================================================================
# Helpers
# =============================================================================
def clean_species_name(feature_name: str) -> str:
    """Remove rCLR_ prefix for cleaner reporting."""
    return feature_name.replace(RCLR_PREFIX, "", 1)


def safe_median(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else np.nan


def mann_whitney_u_from_ranks(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float]:
    """
    Return U_x, expected U_x, variance of U_x with tie correction.

    This is used inside each study block. The block-level U statistics are then
    summed across studies for the van Elteren / blocked Wilcoxon test.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]

    n1 = len(x)
    n2 = len(y)
    n = n1 + n2

    if n1 == 0 or n2 == 0:
        return np.nan, np.nan, np.nan

    combined = np.concatenate([x, y])
    ranks = rankdata(combined, method="average")
    rank_sum_x = ranks[:n1].sum()

    u_x = rank_sum_x - n1 * (n1 + 1) / 2.0
    expected = n1 * n2 / 2.0

    # Tie-corrected variance for Mann-Whitney U.
    # var(U) = n1*n2/12 * [(N+1) - sum(t^3-t)/(N*(N-1))]
    _, counts = np.unique(combined, return_counts=True)
    tie_sum = np.sum(counts**3 - counts)

    if n <= 1:
        variance = np.nan
    else:
        variance = (n1 * n2 / 12.0) * ((n + 1.0) - tie_sum / (n * (n - 1.0)))

    if variance <= 0 or not np.isfinite(variance):
        variance = np.nan

    return float(u_x), float(expected), float(variance)


def van_elteren_test(
    df: pd.DataFrame,
    feature: str,
    group_col: str,
    group_a: str,
    group_b: str,
    block_col: str,
    min_per_group_per_block: int = 2,
) -> Dict[str, object]:
    """
    Stratified Wilcoxon rank-sum test, also called the van Elteren test.

    For each block/study:
        - compute Mann-Whitney U for group_a vs group_b
        - compute its null expectation and tie-corrected variance
    Across blocks:
        Z = sum(U - E[U]) / sqrt(sum(Var[U]))

    Positive Z / effect means group_a tends to have higher feature values.
    """
    total_u_minus_e = 0.0
    total_var = 0.0
    contributing_blocks: List[str] = []
    n_a_total = 0
    n_b_total = 0

    block_details = []

    for block, sub in df.groupby(block_col, observed=True):
        vals_a = sub.loc[sub[group_col] == group_a, feature].to_numpy(dtype=float)
        vals_b = sub.loc[sub[group_col] == group_b, feature].to_numpy(dtype=float)

        vals_a = vals_a[np.isfinite(vals_a)]
        vals_b = vals_b[np.isfinite(vals_b)]

        if len(vals_a) < min_per_group_per_block or len(vals_b) < min_per_group_per_block:
            continue

        u, expected, variance = mann_whitney_u_from_ranks(vals_a, vals_b)
        if not np.isfinite(u) or not np.isfinite(expected) or not np.isfinite(variance):
            continue

        total_u_minus_e += (u - expected)
        total_var += variance
        contributing_blocks.append(str(block))
        n_a_total += len(vals_a)
        n_b_total += len(vals_b)
        block_details.append((str(block), len(vals_a), len(vals_b)))

    if total_var <= 0 or len(contributing_blocks) == 0:
        z = np.nan
        p_value = np.nan
    else:
        z = total_u_minus_e / np.sqrt(total_var)
        p_value = 2.0 * norm.sf(abs(z))

    all_a = df.loc[df[group_col] == group_a, feature].to_numpy(dtype=float)
    all_b = df.loc[df[group_col] == group_b, feature].to_numpy(dtype=float)
    median_a = safe_median(all_a)
    median_b = safe_median(all_b)
    median_diff = median_a - median_b

    return {
        "feature": feature,
        "species": clean_species_name(feature),
        "group_a": group_a,
        "group_b": group_b,
        "n_group_a": n_a_total,
        "n_group_b": n_b_total,
        "n_blocks": len(contributing_blocks),
        "blocks": ";".join(contributing_blocks),
        "z": z,
        "p_value": p_value,
        "median_group_a": median_a,
        "median_group_b": median_b,
        "median_difference_a_minus_b": median_diff,
        "direction": (
            f"higher_in_{group_a}" if median_diff > 0 else
            f"higher_in_{group_b}" if median_diff < 0 else
            "no_median_difference"
        ),
    }


def build_analysis_dataset(x: pd.DataFrame, metadata: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
    """Merge X_explore with metadata and return rCLR feature list."""
    required_x = {SAMPLE_COL}
    required_meta = {SAMPLE_COL, STUDY_COL, CONDITION_COL, SEX_COL}

    missing_x = required_x - set(x.columns)
    missing_meta = required_meta - set(metadata.columns)

    if missing_x:
        raise ValueError(f"X_explore.csv is missing required columns: {sorted(missing_x)}")
    if missing_meta:
        raise ValueError(f"metadata.csv is missing required columns: {sorted(missing_meta)}")

    rclr_features = [c for c in x.columns if c.startswith(RCLR_PREFIX)]
    if not rclr_features:
        raise ValueError(f"No rCLR features found. Expected columns starting with '{RCLR_PREFIX}'.")

    merged = metadata.merge(x[[SAMPLE_COL] + rclr_features], on=SAMPLE_COL, how="inner")

    if merged.empty:
        raise ValueError("Merged metadata + X_explore is empty. Check Sample IDs.")

    return merged, rclr_features


def apply_fdr(results: pd.DataFrame) -> pd.DataFrame:
    """Apply Benjamini-Hochberg FDR separately within each comparison."""
    out = results.copy()
    out["p_adj_BH"] = np.nan
    out["significant_FDR_0_05"] = False

    for comparison, idx in out.groupby("comparison", observed=True).groups.items():
        idx = list(idx)
        pvals = out.loc[idx, "p_value"].astype(float)
        valid = pvals.notna() & np.isfinite(pvals)

        if valid.sum() == 0:
            continue

        valid_idx = pvals.index[valid]
        reject, qvals, _, _ = multipletests(pvals.loc[valid_idx], alpha=FDR_ALPHA, method="fdr_bh")
        out.loc[valid_idx, "p_adj_BH"] = qvals
        out.loc[valid_idx, "significant_FDR_0_05"] = reject

    return out


def run_comparison(
    df: pd.DataFrame,
    rclr_features: List[str],
    comparison_name: str,
    group_col: str,
    group_a: str,
    group_b: str,
    filter_condition: Optional[str] = None,
    filter_sex: Optional[str] = None,
) -> pd.DataFrame:
    """Run blocked Wilcoxon for one comparison across all rCLR features."""
    sub = df.copy()

    if filter_condition is not None:
        sub = sub[sub[CONDITION_COL] == filter_condition].copy()
    if filter_sex is not None:
        sub = sub[sub[SEX_COL] == filter_sex].copy()

    sub = sub[sub[group_col].isin([group_a, group_b])].copy()

    rows = []
    for feature in rclr_features:
        result = van_elteren_test(
            df=sub,
            feature=feature,
            group_col=group_col,
            group_a=group_a,
            group_b=group_b,
            block_col=STUDY_COL,
            min_per_group_per_block=MIN_PER_GROUP_PER_STUDY,
        )
        result["comparison"] = comparison_name
        result["test"] = "van_elteren_blocked_wilcoxon"
        result["blocking_variable"] = STUDY_COL
        if filter_condition is not None:
            result["filter_condition"] = filter_condition
        if filter_sex is not None:
            result["filter_sex"] = filter_sex
        rows.append(result)

    return pd.DataFrame(rows)


def make_summary_plot(results: pd.DataFrame, output_path: Path) -> None:
    """Create one compact figure summarizing Step 3.1 results."""
    summary = (
        results.groupby("comparison", observed=True)
        .agg(
            n_tested=("p_value", lambda s: int(np.isfinite(s).sum())),
            n_significant=("significant_FDR_0_05", "sum"),
            min_q=("p_adj_BH", "min"),
        )
        .reset_index()
    )
    summary["min_minus_log10_q"] = -np.log10(summary["min_q"].clip(lower=1e-300))

    # Keep labels readable.
    labels = summary["comparison"].tolist()
    x = np.arange(len(labels))

    fig, axes = plt.subplots(2, 1, figsize=(max(10, len(labels) * 1.2), 8), constrained_layout=True)

    axes[0].bar(x, summary["n_significant"])
    axes[0].set_ylabel("Significant species\n(FDR < 0.05)")
    axes[0].set_title("Step 3.1 rCLR differential abundance: significant species per comparison")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=35, ha="right")

    axes[1].bar(x, summary["min_minus_log10_q"])
    axes[1].set_ylabel("Best -log10(q)")
    axes[1].set_title("Strongest FDR-adjusted signal per comparison")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=35, ha="right")

    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def write_report(
    report_path: Path,
    df: pd.DataFrame,
    rclr_features: List[str],
    results: pd.DataFrame,
) -> None:
    """Write human-readable report with the main findings."""
    lines = []
    lines.append("STEP 3.1 — rCLR FEATURES: BLOCKED WILCOXON / VAN ELTEREN TEST")
    lines.append("=" * 78)
    lines.append("")

    lines.append("INPUTS")
    lines.append("-" * 78)
    lines.append(f"X_explore path: {X_PATH}")
    lines.append(f"Metadata path:  {METADATA_PATH}")
    lines.append(f"Merged samples: {df.shape[0]}")
    lines.append(f"rCLR features tested: {len(rclr_features)}")
    lines.append(f"Blocking variable: {STUDY_COL}")
    lines.append(f"FDR method: Benjamini-Hochberg, separately within each comparison")
    lines.append(f"FDR threshold: {FDR_ALPHA}")
    lines.append(f"Minimum per group per study block: {MIN_PER_GROUP_PER_STUDY}")
    lines.append(f"Minimum contributing study blocks required: {MIN_CONTRIBUTING_STUDIES}")
    lines.append("")

    lines.append("SAMPLE COUNTS")
    lines.append("-" * 78)
    lines.append("Condition counts:")
    for k, v in df[CONDITION_COL].value_counts(dropna=False).items():
        lines.append(f"  {k}: {v}")
    lines.append("Sex counts:")
    for k, v in df[SEX_COL].value_counts(dropna=False).items():
        lines.append(f"  {k}: {v}")
    lines.append(f"Number of studies: {df[STUDY_COL].nunique()}")
    lines.append("")

    lines.append("INTERPRETATION OF EFFECT DIRECTION")
    lines.append("-" * 78)
    lines.append("median_difference_a_minus_b = median(group_a) - median(group_b)")
    lines.append("Positive value means the rCLR feature is higher in group_a.")
    lines.append("Negative value means the rCLR feature is higher in group_b.")
    lines.append("For condition comparisons, group_a is the disease/advanced group where applicable.")
    lines.append("")

    lines.append("COMPARISON SUMMARY")
    lines.append("-" * 78)
    summary = (
        results.groupby("comparison", observed=True)
        .agg(
            n_tested=("p_value", lambda s: int(np.isfinite(s).sum())),
            n_significant=("significant_FDR_0_05", "sum"),
            min_p=("p_value", "min"),
            min_q=("p_adj_BH", "min"),
            median_blocks=("n_blocks", "median"),
        )
        .reset_index()
        .sort_values("comparison")
    )

    for _, row in summary.iterrows():
        lines.append(f"{row['comparison']}:")
        lines.append(f"  Tested features: {int(row['n_tested'])} / {len(rclr_features)}")
        lines.append(f"  Significant species at FDR < {FDR_ALPHA}: {int(row['n_significant'])}")
        lines.append(f"  Minimum raw p-value: {row['min_p']:.4e}" if pd.notna(row['min_p']) else "  Minimum raw p-value: NA")
        lines.append(f"  Minimum adjusted q-value: {row['min_q']:.4e}" if pd.notna(row['min_q']) else "  Minimum adjusted q-value: NA")
        lines.append(f"  Median contributing study blocks: {row['median_blocks']:.1f}")
        lines.append("")

        top = (
            results[(results["comparison"] == row["comparison"]) & results["p_adj_BH"].notna()]
            .sort_values(["p_adj_BH", "p_value"])
            .head(10)
        )
        if not top.empty:
            lines.append("  Top species by FDR-adjusted q-value:")
            for _, r in top.iterrows():
                q = r["p_adj_BH"]
                p = r["p_value"]
                md = r["median_difference_a_minus_b"]
                direction = r["direction"]
                lines.append(
                    f"    - {r['species']} | q={q:.4e}, p={p:.4e}, "
                    f"median_diff={md:.4f}, {direction}, blocks={int(r['n_blocks'])}"
                )
            lines.append("")

    lines.append("FILES WRITTEN")
    lines.append("-" * 78)
    lines.append("step3_1_blocked_wilcoxon_all_results.csv")
    lines.append("step3_1_blocked_wilcoxon_significant_results.csv")
    lines.append("step3_1_blocked_wilcoxon_summary.png")
    lines.append("step3_1_final_report.txt")
    lines.append("")

    lines.append("IMPORTANT NOTES")
    lines.append("-" * 78)
    lines.append("1. This script does not merge Adenoma with CRC or Control.")
    lines.append("2. Adenoma is analyzed as a separate biological/clinical condition.")
    lines.append("3. The blocked Wilcoxon test controls for study by comparing groups within study blocks.")
    lines.append("4. Features with too few valid study blocks are not interpretable; check n_blocks.")
    lines.append("5. These results should later be compared with SHAP features in Phase 5.")

    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not X_PATH.exists():
        raise FileNotFoundError(f"Could not find X_explore.csv at: {X_PATH}")
    if not METADATA_PATH.exists():
        raise FileNotFoundError(f"Could not find metadata.csv at: {METADATA_PATH}")

    x = pd.read_csv(X_PATH)
    metadata = pd.read_csv(METADATA_PATH)

    df, rclr_features = build_analysis_dataset(x, metadata)

    # Define comparisons for Step 3.1.
    # The first group listed is group_a, so median_difference_a_minus_b is interpreted as group_a - group_b.
    comparison_specs = [
        # CRC vs Control overall
        {
            "comparison_name": "CRC_vs_Control_all",
            "group_col": CONDITION_COL,
            "group_a": "CRC",
            "group_b": "Control",
        },
        # CRC vs Control stratified by sex
        {
            "comparison_name": "CRC_vs_Control_Male",
            "group_col": CONDITION_COL,
            "group_a": "CRC",
            "group_b": "Control",
            "filter_sex": "Male",
        },
        {
            "comparison_name": "CRC_vs_Control_Female",
            "group_col": CONDITION_COL,
            "group_a": "CRC",
            "group_b": "Control",
            "filter_sex": "Female",
        },
        # Male vs Female stratified by condition
        {
            "comparison_name": "Male_vs_Female_within_Control",
            "group_col": SEX_COL,
            "group_a": "Male",
            "group_b": "Female",
            "filter_condition": "Control",
        },
        {
            "comparison_name": "Male_vs_Female_within_Adenoma",
            "group_col": SEX_COL,
            "group_a": "Male",
            "group_b": "Female",
            "filter_condition": "Adenoma",
        },
        {
            "comparison_name": "Male_vs_Female_within_CRC",
            "group_col": SEX_COL,
            "group_a": "Male",
            "group_b": "Female",
            "filter_condition": "CRC",
        },
        # Adenoma comparisons
        {
            "comparison_name": "Adenoma_vs_Control",
            "group_col": CONDITION_COL,
            "group_a": "Adenoma",
            "group_b": "Control",
        },
        {
            "comparison_name": "CRC_vs_Adenoma",
            "group_col": CONDITION_COL,
            "group_a": "CRC",
            "group_b": "Adenoma",
        },
    ]

    all_results = []
    for spec in comparison_specs:
        res = run_comparison(
            df=df,
            rclr_features=rclr_features,
            comparison_name=spec["comparison_name"],
            group_col=spec["group_col"],
            group_a=spec["group_a"],
            group_b=spec["group_b"],
            filter_condition=spec.get("filter_condition"),
            filter_sex=spec.get("filter_sex"),
        )
        all_results.append(res)

    results = pd.concat(all_results, ignore_index=True)

    # Mark results with too few contributing study blocks as not tested.
    too_few_blocks = results["n_blocks"] < MIN_CONTRIBUTING_STUDIES
    results.loc[too_few_blocks, ["z", "p_value"]] = np.nan
    results.loc[too_few_blocks, "note"] = "too_few_contributing_study_blocks"
    results.loc[~too_few_blocks, "note"] = "ok"

    results = apply_fdr(results)
    results = results.sort_values(["comparison", "p_adj_BH", "p_value", "species"], na_position="last")

    all_results_path = OUTPUT_DIR / "step31_blocked_wilcoxon_all_results.csv"
    sig_results_path = OUTPUT_DIR / "step31_blocked_wilcoxon_significant_results.csv"
    plot_path = OUTPUT_DIR / "step31_blocked_wilcoxon_summary.png"
    report_path = OUTPUT_DIR / "step31_final_report.txt"

    results.to_csv(all_results_path, index=False)
    significant = results[results["significant_FDR_0_05"]].copy()
    significant.to_csv(sig_results_path, index=False)

    make_summary_plot(results, plot_path)
    write_report(report_path, df, rclr_features, results)

    print("Step 3.1 complete.")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Wrote: {all_results_path.name}")
    print(f"Wrote: {sig_results_path.name}")
    print(f"Wrote: {plot_path.name}")
    print(f"Wrote: {report_path.name}")


if __name__ == "__main__":
    main()
