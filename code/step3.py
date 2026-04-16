"""
STEP 3 — DIFFERENTIAL ABUNDANCE ANALYSIS

This script identifies microbial genera that differ significantly between
disease states (CRC vs Control) and between sexes, following the blueprint
analysis design.

Main functionalities:
- Loads merged dataset (metadata + CLR-transformed genus features)
- Filters samples to valid labels:
  - Condition: CRC, Control, Adenoma
  - Gender: Male, Female
- Defines biologically relevant contrasts:
  - CRC vs Control (within Male)
  - CRC vs Control (within Female)
  - Male vs Female (within Control)
  - Male vs Female (within CRC)
- Performs statistical testing:
  - Wilcoxon rank-sum test (non-parametric)
  - Benjamini–Hochberg FDR correction
- Computes effect sizes:
  - mean CLR difference
  - rank-biserial correlation
- Generates outputs:
  - full results tables per contrast
  - summary tables (number of significant taxa)
  - volcano plots (effect size vs significance)
- Identifies top discriminative genera per contrast
- Performs cross-contrast comparison:
  - builds heatmap of shared vs contrast-specific taxa
- Optional robustness check:
  - evaluates consistency of effects across studies

Blueprint steps covered:
- Statistical testing of microbiome differences
- Identification of disease-associated taxa
- Assessment of sex-specific microbiome effects
- Comparison across biological contrasts
- Validation of findings consistency across cohorts
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ranksums, mannwhitneyu

warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================
# Configuration
# ============================================================

ROOT = Path(".")
FEATURES_PATH = ROOT / "harmonized" / "unified_genera_clr.csv"
METADATA_PATH = ROOT / "harmonized" / "unified_metadata.csv"

OUTPUT_DIR = ROOT / "outputs" / "step3_differential_abundance"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_COL = "Sample"
STUDY_COL = "Study"
TARGET_COL = "Condition"
SEX_COL = "Gender"

RANDOM_STATE = 42
ALPHA = 0.05
TOP_N_VOLCANO_LABELS = 12


# ============================================================
# Helpers
# ============================================================

def print_header(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def normalize_labels(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.replace(r"\s+", " ", regex=True)
    )


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR correction.
    """
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)

    order = np.argsort(pvals)
    ranked_pvals = pvals[order]

    adjusted = np.empty(n, dtype=float)
    cumulative_min = 1.0

    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked_pvals[i] * n / rank
        cumulative_min = min(cumulative_min, val)
        adjusted[i] = cumulative_min

    adjusted = np.clip(adjusted, 0, 1)

    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return out


def rank_biserial_from_mwu(x1: np.ndarray, x0: np.ndarray) -> float:
    """
    Rank-biserial correlation from Mann-Whitney U.
    Positive means group1 tends to have larger values than group0.
    """
    n1 = len(x1)
    n0 = len(x0)
    if n1 == 0 or n0 == 0:
        return np.nan

    u_stat, _ = mannwhitneyu(x1, x0, alternative="two-sided")
    return (2 * u_stat) / (n1 * n0) - 1


def save_df(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)
    print(f"[saved] {path}")


def volcano_plot(
    results: pd.DataFrame,
    contrast_name: str,
    out_path: Path,
    alpha: float = 0.05,
    label_top_n: int = 12,
) -> None:
    df = results.copy()
    df["neg_log10_fdr"] = -np.log10(df["fdr_bh"].clip(lower=1e-300))

    sig_mask = (df["fdr_bh"] < alpha)

    plt.figure(figsize=(10, 8))
    plt.scatter(df["mean_diff_clr"], df["neg_log10_fdr"], s=18, alpha=0.7)

    if sig_mask.any():
        plt.scatter(
            df.loc[sig_mask, "mean_diff_clr"],
            df.loc[sig_mask, "neg_log10_fdr"],
            s=20,
            alpha=0.9,
        )

    plt.axhline(-np.log10(alpha), linestyle="--", linewidth=1)
    plt.xlabel("Mean CLR difference (group1 - group0)")
    plt.ylabel("-log10(FDR)")
    plt.title(f"Volcano plot: {contrast_name}")

    # Label top features by smallest FDR then largest absolute mean difference
    label_df = (
        df.sort_values(["fdr_bh", "abs_mean_diff"], ascending=[True, False])
        .head(label_top_n)
    )

    for _, row in label_df.iterrows():
        plt.text(
            row["mean_diff_clr"],
            row["neg_log10_fdr"],
            row["feature"],
            fontsize=8,
            alpha=0.85,
        )

    plt.tight_layout()
    plt.savefig(out_path, dpi=220)
    plt.close()
    print(f"[saved] {out_path}")


def comparison_heatmap(
    matrix: pd.DataFrame,
    title: str,
    out_path: Path,
) -> None:
    """
    Simple matplotlib heatmap without seaborn.
    """
    if matrix.empty:
        print(f"[skip] no data for heatmap: {title}")
        return

    fig, ax = plt.subplots(figsize=(10, max(6, 0.35 * len(matrix))))
    im = ax.imshow(matrix.values, aspect="auto")

    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index)

    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Mean CLR difference")

    plt.tight_layout()
    plt.savefig(out_path, dpi=220)
    plt.close()
    print(f"[saved] {out_path}")


# ============================================================
# Load data
# ============================================================

print_header("1. LOAD DATA")

features = pd.read_csv(FEATURES_PATH)
metadata = pd.read_csv(METADATA_PATH)

features[SAMPLE_COL] = normalize_labels(features[SAMPLE_COL])
metadata[SAMPLE_COL] = normalize_labels(metadata[SAMPLE_COL])

metadata[TARGET_COL] = normalize_labels(metadata[TARGET_COL])
metadata[SEX_COL] = normalize_labels(metadata[SEX_COL])
metadata[STUDY_COL] = normalize_labels(metadata[STUDY_COL])

merged = metadata.merge(features, on=SAMPLE_COL, how="inner")

feature_cols = [c for c in features.columns if c != SAMPLE_COL]
merged[feature_cols] = merged[feature_cols].apply(pd.to_numeric, errors="coerce")

if merged[feature_cols].isna().sum().sum() > 0:
    raise ValueError("Feature matrix contains NaNs after numeric coercion.")

print(f"Merged shape: {merged.shape}")
print(f"Number of genus features: {len(feature_cols)}")
print("Condition counts:")
print(merged[TARGET_COL].value_counts(dropna=False))
print("Gender counts:")
print(merged[SEX_COL].value_counts(dropna=False))


# ============================================================
# Restrict to labels of interest
# ============================================================

print_header("2. FILTER LABELS")

valid_conditions = {"CRC", "Control", "Adenoma"}
valid_sexes = {"Male", "Female"}

df = merged.loc[
    merged[TARGET_COL].isin(valid_conditions) &
    merged[SEX_COL].isin(valid_sexes)
].copy()

print(f"Filtered shape: {df.shape}")
print("Condition x Gender:")
print(pd.crosstab(df[TARGET_COL], df[SEX_COL], dropna=False))

save_df(
    pd.crosstab(df[TARGET_COL], df[SEX_COL], dropna=False).reset_index(),
    OUTPUT_DIR / "condition_by_gender_counts.csv"
)


# ============================================================
# Differential abundance engine
# ============================================================

def run_contrast(
    data: pd.DataFrame,
    group1_mask: pd.Series,
    group0_mask: pd.Series,
    contrast_name: str,
    group1_label: str,
    group0_label: str,
    feature_columns: List[str],
) -> pd.DataFrame:
    """
    Runs one differential abundance contrast in CLR space.
    """
    d1 = data.loc[group1_mask].copy()
    d0 = data.loc[group0_mask].copy()

    n1 = len(d1)
    n0 = len(d0)

    if n1 == 0 or n0 == 0:
        raise ValueError(f"{contrast_name}: one group has zero samples.")

    rows = []

    for feat in feature_columns:
        x1 = d1[feat].to_numpy(dtype=float)
        x0 = d0[feat].to_numpy(dtype=float)

        stat, pval = ranksums(x1, x0)
        rbc = rank_biserial_from_mwu(x1, x0)

        mean1 = np.mean(x1)
        mean0 = np.mean(x0)
        med1 = np.median(x1)
        med0 = np.median(x0)

        rows.append({
            "contrast": contrast_name,
            "group1_label": group1_label,
            "group0_label": group0_label,
            "n_group1": n1,
            "n_group0": n0,
            "feature": feat,
            "mean_group1": mean1,
            "mean_group0": mean0,
            "median_group1": med1,
            "median_group0": med0,
            "mean_diff_clr": mean1 - mean0,
            "median_diff_clr": med1 - med0,
            "abs_mean_diff": abs(mean1 - mean0),
            "ranksum_stat": stat,
            "p_value": pval,
            "rank_biserial": rbc,
        })

    res = pd.DataFrame(rows)
    res["fdr_bh"] = bh_fdr(res["p_value"].values)
    res["significant_fdr_0_05"] = res["fdr_bh"] < ALPHA
    res = res.sort_values(["fdr_bh", "abs_mean_diff"], ascending=[True, False]).reset_index(drop=True)

    return res


# ============================================================
# Define contrasts from blueprint
# ============================================================

print_header("3. RUN CONTRASTS")

contrasts: List[Tuple[str, pd.Series, pd.Series, str, str]] = []

# 1) CRC vs Control within Male
contrasts.append((
    "CRC_vs_Control__within_Male",
    (df[TARGET_COL] == "CRC") & (df[SEX_COL] == "Male"),
    (df[TARGET_COL] == "Control") & (df[SEX_COL] == "Male"),
    "CRC_Male",
    "Control_Male",
))

# 2) CRC vs Control within Female
contrasts.append((
    "CRC_vs_Control__within_Female",
    (df[TARGET_COL] == "CRC") & (df[SEX_COL] == "Female"),
    (df[TARGET_COL] == "Control") & (df[SEX_COL] == "Female"),
    "CRC_Female",
    "Control_Female",
))

# 3) Male vs Female within Control
contrasts.append((
    "Male_vs_Female__within_Control",
    (df[TARGET_COL] == "Control") & (df[SEX_COL] == "Male"),
    (df[TARGET_COL] == "Control") & (df[SEX_COL] == "Female"),
    "Male_Control",
    "Female_Control",
))

# 4) Male vs Female within CRC
contrasts.append((
    "Male_vs_Female__within_CRC",
    (df[TARGET_COL] == "CRC") & (df[SEX_COL] == "Male"),
    (df[TARGET_COL] == "CRC") & (df[SEX_COL] == "Female"),
    "Male_CRC",
    "Female_CRC",
))

all_results = []

for contrast_name, g1_mask, g0_mask, g1_label, g0_label in contrasts:
    print(f"\nRunning: {contrast_name}")
    print(f"  group1 n = {int(g1_mask.sum())}")
    print(f"  group0 n = {int(g0_mask.sum())}")

    res = run_contrast(
        data=df,
        group1_mask=g1_mask,
        group0_mask=g0_mask,
        contrast_name=contrast_name,
        group1_label=g1_label,
        group0_label=g0_label,
        feature_columns=feature_cols,
    )

    all_results.append(res)

    save_df(res, OUTPUT_DIR / f"{contrast_name}.csv")

    summary = pd.DataFrame({
        "contrast": [contrast_name],
        "n_group1": [int(g1_mask.sum())],
        "n_group0": [int(g0_mask.sum())],
        "n_features_tested": [len(res)],
        "n_significant_fdr_0_05": [int(res["significant_fdr_0_05"].sum())],
        "top_feature_by_fdr": [res.iloc[0]["feature"]],
        "top_feature_fdr": [res.iloc[0]["fdr_bh"]],
        "top_feature_mean_diff": [res.iloc[0]["mean_diff_clr"]],
    })
    save_df(summary, OUTPUT_DIR / f"{contrast_name}__summary.csv")

    volcano_plot(
        results=res,
        contrast_name=contrast_name,
        out_path=OUTPUT_DIR / f"{contrast_name}__volcano.png",
        alpha=ALPHA,
        label_top_n=TOP_N_VOLCANO_LABELS,
    )

all_results_df = pd.concat(all_results, ignore_index=True)
save_df(all_results_df, OUTPUT_DIR / "all_contrasts_combined.csv")


# ============================================================
# Summary tables
# ============================================================

print_header("4. SUMMARIES")

contrast_summary_rows = []
for contrast_name in all_results_df["contrast"].unique():
    tmp = all_results_df.loc[all_results_df["contrast"] == contrast_name].copy()
    tmp_sig = tmp.loc[tmp["significant_fdr_0_05"]].copy()

    contrast_summary_rows.append({
        "contrast": contrast_name,
        "n_tested_features": len(tmp),
        "n_significant_fdr_0_05": len(tmp_sig),
        "best_feature": tmp.iloc[0]["feature"],
        "best_fdr": tmp.iloc[0]["fdr_bh"],
        "best_mean_diff_clr": tmp.iloc[0]["mean_diff_clr"],
    })

contrast_summary = pd.DataFrame(contrast_summary_rows)
save_df(contrast_summary, OUTPUT_DIR / "contrast_summary.csv")
print(contrast_summary)


# Significant-only compact table
sig_table = (
    all_results_df.loc[all_results_df["significant_fdr_0_05"]]
    .sort_values(["contrast", "fdr_bh", "abs_mean_diff"], ascending=[True, True, False])
    .reset_index(drop=True)
)
save_df(sig_table, OUTPUT_DIR / "significant_results_fdr_0_05.csv")


# Top features per contrast
top_per_contrast = (
    all_results_df.sort_values(["contrast", "fdr_bh", "abs_mean_diff"], ascending=[True, True, False])
    .groupby("contrast", as_index=False)
    .head(20)
    .reset_index(drop=True)
)
save_df(top_per_contrast, OUTPUT_DIR / "top20_features_per_contrast.csv")


# ============================================================
# Cross-contrast comparison heatmap
# ============================================================

print_header("5. CROSS-CONTRAST COMPARISON")

# Collect union of top features across contrasts
top_features_union = (
    all_results_df.sort_values(["contrast", "fdr_bh", "abs_mean_diff"], ascending=[True, True, False])
    .groupby("contrast")
    .head(12)["feature"]
    .unique()
    .tolist()
)

heatmap_df = (
    all_results_df.loc[all_results_df["feature"].isin(top_features_union), ["contrast", "feature", "mean_diff_clr"]]
    .pivot(index="feature", columns="contrast", values="mean_diff_clr")
    .fillna(0.0)
)

# Order rows by overall magnitude
row_order = heatmap_df.abs().sum(axis=1).sort_values(ascending=False).index
heatmap_df = heatmap_df.loc[row_order]

save_df(
    heatmap_df.reset_index(),
    OUTPUT_DIR / "top_features_cross_contrast_mean_diff_matrix.csv"
)

comparison_heatmap(
    matrix=heatmap_df,
    title="Top genera across contrasts (mean CLR difference)",
    out_path=OUTPUT_DIR / "top_features_cross_contrast_heatmap.png",
)


# ============================================================
# Optional per-study robustness screen
# ============================================================

print_header("6. OPTIONAL STUDY-LEVEL ROBUSTNESS SCREEN")

robustness_rows = []

for contrast_name, g1_mask, g0_mask, g1_label, g0_label in contrasts:
    res = all_results_df.loc[all_results_df["contrast"] == contrast_name].copy()
    top_feats = res.head(10)["feature"].tolist()

    for feat in top_feats:
        for study, d_study in df.groupby(STUDY_COL):
            s_g1 = d_study.loc[g1_mask.loc[d_study.index], feat]
            s_g0 = d_study.loc[g0_mask.loc[d_study.index], feat]

            if len(s_g1) >= 5 and len(s_g0) >= 5:
                robustness_rows.append({
                    "contrast": contrast_name,
                    "study": study,
                    "feature": feat,
                    "n_group1_in_study": len(s_g1),
                    "n_group0_in_study": len(s_g0),
                    "mean_diff_clr_in_study": s_g1.mean() - s_g0.mean(),
                })

robustness_df = pd.DataFrame(robustness_rows)
save_df(robustness_df, OUTPUT_DIR / "study_level_effect_direction_screen.csv")


# ============================================================
# End
# ============================================================

print_header("DONE")

print("Saved outputs in:")
print(OUTPUT_DIR.resolve())

print("\nMain files to inspect next:")
print("- contrast_summary.csv")
print("- significant_results_fdr_0_05.csv")
print("- top20_features_per_contrast.csv")
print("- *_volcano.png")
print("- top_features_cross_contrast_heatmap.png")
