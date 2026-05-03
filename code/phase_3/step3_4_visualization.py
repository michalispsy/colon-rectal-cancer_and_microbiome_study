#!/usr/bin/env python3
"""
Step 3.4 — Visualization for Differential Abundance Analysis

This script is self-contained. It does NOT depend on outputs from Steps 3.1, 3.2, or 3.3.
It reads only:
    - X_explore.csv
    - metadata.csv

It produces a compact visualization/reporting set for Phase 3:
    - Volcano-style plots for rCLR features using simple effect sizes and Mann-Whitney p-values
    - Presence/absence volcano-style plots for PA features using prevalence differences and Fisher exact p-values
    - Heatmap of top differential species by condition and sex
    - Forest-style plot of per-study effects for top CRC-vs-Control species
    - CSV list of top differentially abundant species for later SHAP comparison
    - Final text report

Important:
    This is a visualization/reporting step. The formal blocked Wilcoxon, CMH, and random-effects
    meta-analysis scripts should still be treated as the main statistical evidence.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import mannwhitneyu, fisher_exact
from statsmodels.stats.multitest import multipletests


# =============================================================================
# Paths — project structure
# =============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
CODE_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = CODE_ROOT.parent

DATA_DIR = (PROJECT_ROOT / "data" / "crc_study_final_data" / "species_level").resolve()
OUTPUT_DIR = (PROJECT_ROOT / "OUTPUTS" / "phase_3" / "3.4").resolve()

X_PATH = (DATA_DIR / "X_explore.csv").resolve()
METADATA_PATH = (DATA_DIR / "metadata.csv").resolve()

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# User-adjustable parameters
# =============================================================================

ALPHA_FDR = 0.05
TOP_N_SPECIES = 25
MIN_GROUP_N = 10
MIN_STUDIES_FOR_FOREST = 3
EPS = 1e-12

# Comparisons used for volcano plots and top-species list.
# Tuple structure: (comparison_name, positive_group, negative_group)
CONDITION_COMPARISONS = [
    ("CRC_vs_Control", "CRC", "Control"),
    ("Adenoma_vs_Control", "Adenoma", "Control"),
    ("CRC_vs_Adenoma", "CRC", "Adenoma"),
]


# =============================================================================
# Utility functions
# =============================================================================

def clean_species_name(feature_name: str) -> str:
    """Remove PA_/rCLR_ prefix for cleaner labels."""
    if feature_name.startswith("rCLR_"):
        return feature_name.replace("rCLR_", "", 1)
    if feature_name.startswith("PA_"):
        return feature_name.replace("PA_", "", 1)
    return feature_name


def safe_neg_log10(q: float) -> float:
    if pd.isna(q):
        return np.nan
    return -np.log10(max(float(q), EPS))


def infer_column(df: pd.DataFrame, candidates: Iterable[str], required: bool = True) -> Optional[str]:
    """Find the first matching column from a list of possible names."""
    lower_to_original = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        if cand.lower() in lower_to_original:
            return lower_to_original[cand.lower()]
    if required:
        raise ValueError(f"Could not find required column. Tried: {list(candidates)}")
    return None


def standardize_condition_values(series: pd.Series) -> pd.Series:
    """Normalize common condition spellings."""
    mapping = {
        "crc": "CRC",
        "cancer": "CRC",
        "carcinoma": "CRC",
        "control": "Control",
        "healthy": "Control",
        "normal": "Control",
        "adenoma": "Adenoma",
        "adenomas": "Adenoma",
        "advanced adenoma": "Adenoma",
    }
    return series.astype(str).str.strip().str.lower().map(lambda x: mapping.get(x, x.title()))


def standardize_sex_values(series: pd.Series) -> pd.Series:
    """Normalize common sex spellings."""
    mapping = {
        "m": "Male",
        "male": "Male",
        "man": "Male",
        "f": "Female",
        "female": "Female",
        "woman": "Female",
    }
    return series.astype(str).str.strip().str.lower().map(lambda x: mapping.get(x, x.title()))


def load_data() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, str]]:
    """Load X_explore and metadata, infer important columns, merge them."""
    if not X_PATH.exists():
        raise FileNotFoundError(f"X_explore.csv not found at: {X_PATH}")
    if not METADATA_PATH.exists():
        raise FileNotFoundError(f"metadata.csv not found at: {METADATA_PATH}")

    X = pd.read_csv(X_PATH)
    metadata = pd.read_csv(METADATA_PATH)

    sample_col_x = infer_column(X, ["Sample", "sample", "sample_id", "Sample_ID"])
    sample_col_meta = infer_column(metadata, ["Sample", "sample", "sample_id", "Sample_ID"])
    condition_col = infer_column(metadata, ["study_condition", "condition", "Condition", "diagnosis"])
    study_col = infer_column(metadata, ["study_name", "study", "Study", "cohort", "dataset"])
    sex_col = infer_column(metadata, ["sex", "gender", "Sex", "Gender"], required=False)

    # Rename key columns consistently.
    if sample_col_x != "Sample":
        X = X.rename(columns={sample_col_x: "Sample"})
    if sample_col_meta != "Sample":
        metadata = metadata.rename(columns={sample_col_meta: "Sample"})
    if condition_col != "condition":
        metadata = metadata.rename(columns={condition_col: "condition"})
    if study_col != "study":
        metadata = metadata.rename(columns={study_col: "study"})
    if sex_col is not None and sex_col != "sex":
        metadata = metadata.rename(columns={sex_col: "sex"})

    metadata["condition"] = standardize_condition_values(metadata["condition"])
    if "sex" in metadata.columns:
        metadata["sex"] = standardize_sex_values(metadata["sex"])
    else:
        metadata["sex"] = "Unknown"

    merged = X.merge(metadata, on="Sample", how="inner")

    columns = {
        "sample": "Sample",
        "condition": "condition",
        "study": "study",
        "sex": "sex",
    }

    return X, metadata, merged, columns


def get_feature_columns(df: pd.DataFrame) -> Tuple[List[str], List[str]]:
    rclr_cols = [c for c in df.columns if c.startswith("rCLR_")]
    pa_cols = [c for c in df.columns if c.startswith("PA_")]

    if len(rclr_cols) == 0:
        raise ValueError("No rCLR_ columns found in X_explore.csv.")
    if len(pa_cols) == 0:
        raise ValueError("No PA_ columns found in X_explore.csv.")

    return rclr_cols, pa_cols


def benjamini_hochberg(results: pd.DataFrame, p_col: str = "p_value") -> pd.DataFrame:
    out = results.copy()
    valid = out[p_col].notna()
    out["q_value"] = np.nan
    if valid.sum() > 0:
        out.loc[valid, "q_value"] = multipletests(out.loc[valid, p_col], method="fdr_bh")[1]
    return out


# =============================================================================
# Differential summaries for visualization
# =============================================================================

def rclr_volcano_stats(
    df: pd.DataFrame,
    rclr_cols: List[str],
    comparison_name: str,
    positive_group: str,
    negative_group: str,
) -> pd.DataFrame:
    """
    Compute simple rCLR effect summaries for volcano visualization.

    x-axis: mean rCLR difference = mean(positive_group) - mean(negative_group)
    y-axis: -log10(FDR q-value)

    This is not the formal blocked Wilcoxon test. It is a visualization-oriented summary.
    """
    sub = df[df["condition"].isin([positive_group, negative_group])].copy()
    g1_mask = sub["condition"] == positive_group
    g0_mask = sub["condition"] == negative_group

    records = []
    for col in rclr_cols:
        x1 = pd.to_numeric(sub.loc[g1_mask, col], errors="coerce").dropna().to_numpy()
        x0 = pd.to_numeric(sub.loc[g0_mask, col], errors="coerce").dropna().to_numpy()

        if len(x1) < MIN_GROUP_N or len(x0) < MIN_GROUP_N:
            p_value = np.nan
            effect = np.nan
        else:
            effect = float(np.mean(x1) - np.mean(x0))
            try:
                p_value = float(mannwhitneyu(x1, x0, alternative="two-sided").pvalue)
            except ValueError:
                p_value = np.nan

        records.append(
            {
                "comparison": comparison_name,
                "feature_type": "rCLR",
                "feature": col,
                "species": clean_species_name(col),
                "positive_group": positive_group,
                "negative_group": negative_group,
                "n_positive": len(x1),
                "n_negative": len(x0),
                "effect": effect,
                "effect_name": "mean_rCLR_difference",
                "p_value": p_value,
            }
        )

    out = pd.DataFrame(records)
    out = benjamini_hochberg(out)
    out["minus_log10_q"] = out["q_value"].apply(safe_neg_log10)
    out["significant_fdr_0_05"] = out["q_value"] < ALPHA_FDR
    return out


def pa_volcano_stats(
    df: pd.DataFrame,
    pa_cols: List[str],
    comparison_name: str,
    positive_group: str,
    negative_group: str,
) -> pd.DataFrame:
    """
    Compute PA effect summaries for volcano visualization.

    x-axis: prevalence difference = prevalence(positive_group) - prevalence(negative_group)
    y-axis: -log10(FDR q-value)

    This is not the formal CMH test. It is a visualization-oriented summary.
    """
    sub = df[df["condition"].isin([positive_group, negative_group])].copy()
    g1_mask = sub["condition"] == positive_group
    g0_mask = sub["condition"] == negative_group

    records = []
    for col in pa_cols:
        x1 = pd.to_numeric(sub.loc[g1_mask, col], errors="coerce").fillna(0).astype(int).to_numpy()
        x0 = pd.to_numeric(sub.loc[g0_mask, col], errors="coerce").fillna(0).astype(int).to_numpy()

        if len(x1) < MIN_GROUP_N or len(x0) < MIN_GROUP_N:
            p_value = np.nan
            effect = np.nan
            prev1 = np.nan
            prev0 = np.nan
        else:
            present1 = int(np.sum(x1 > 0))
            absent1 = int(len(x1) - present1)
            present0 = int(np.sum(x0 > 0))
            absent0 = int(len(x0) - present0)

            prev1 = present1 / len(x1)
            prev0 = present0 / len(x0)
            effect = prev1 - prev0

            table = np.array([[present1, absent1], [present0, absent0]])
            try:
                p_value = float(fisher_exact(table, alternative="two-sided").pvalue)
            except ValueError:
                p_value = np.nan

        records.append(
            {
                "comparison": comparison_name,
                "feature_type": "PA",
                "feature": col,
                "species": clean_species_name(col),
                "positive_group": positive_group,
                "negative_group": negative_group,
                "n_positive": len(x1),
                "n_negative": len(x0),
                "prevalence_positive": prev1,
                "prevalence_negative": prev0,
                "effect": effect,
                "effect_name": "prevalence_difference",
                "p_value": p_value,
            }
        )

    out = pd.DataFrame(records)
    out = benjamini_hochberg(out)
    out["minus_log10_q"] = out["q_value"].apply(safe_neg_log10)
    out["significant_fdr_0_05"] = out["q_value"] < ALPHA_FDR
    return out


def compute_all_visual_stats(df: pd.DataFrame, rclr_cols: List[str], pa_cols: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rclr_results = []
    pa_results = []
    for comparison_name, positive, negative in CONDITION_COMPARISONS:
        rclr_results.append(rclr_volcano_stats(df, rclr_cols, comparison_name, positive, negative))
        pa_results.append(pa_volcano_stats(df, pa_cols, comparison_name, positive, negative))
    return pd.concat(rclr_results, ignore_index=True), pd.concat(pa_results, ignore_index=True)


# =============================================================================
# Plotting
# =============================================================================

def plot_volcano_grid(results: pd.DataFrame, title: str, x_label: str, output_path: Path) -> None:
    comparisons = [c[0] for c in CONDITION_COMPARISONS]
    n = len(comparisons)

    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5), constrained_layout=True)
    if n == 1:
        axes = [axes]

    for ax, comp in zip(axes, comparisons):
        sub = results[results["comparison"] == comp].copy()
        sub = sub.replace([np.inf, -np.inf], np.nan).dropna(subset=["effect", "minus_log10_q"])

        nonsig = sub[~sub["significant_fdr_0_05"]]
        sig = sub[sub["significant_fdr_0_05"]]

        ax.scatter(nonsig["effect"], nonsig["minus_log10_q"], s=20, alpha=0.45)
        ax.scatter(sig["effect"], sig["minus_log10_q"], s=25, alpha=0.85)

        if len(sub) > 0:
            ax.axhline(-np.log10(ALPHA_FDR), linestyle="--", linewidth=1)
            ax.axvline(0, linestyle="--", linewidth=1)

        # Label top significant species, or top by q-value if none significant.
        label_df = sig.sort_values("q_value").head(6)
        if label_df.empty:
            label_df = sub.sort_values("q_value").head(4)

        for _, row in label_df.iterrows():
            ax.text(row["effect"], row["minus_log10_q"], row["species"], fontsize=7)

        ax.set_title(comp.replace("_", " "))
        ax.set_xlabel(x_label)
        ax.set_ylabel("-log10(FDR q-value)")

    fig.suptitle(title, fontsize=14)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def select_top_species_for_heatmap(rclr_results: pd.DataFrame, pa_results: pd.DataFrame, top_n: int = TOP_N_SPECIES) -> List[str]:
    """Select top species across rCLR and PA visualization stats."""
    combined = pd.concat(
        [
            rclr_results[["species", "comparison", "q_value", "effect", "feature_type"]],
            pa_results[["species", "comparison", "q_value", "effect", "feature_type"]],
        ],
        ignore_index=True,
    )
    combined = combined.dropna(subset=["q_value"])
    combined["abs_effect"] = combined["effect"].abs()

    # One row per species: prioritize smallest q-value, then larger absolute effect.
    ranked = (
        combined.sort_values(["q_value", "abs_effect"], ascending=[True, False])
        .drop_duplicates("species")
        .head(top_n)
    )
    return ranked["species"].tolist()


def plot_condition_sex_heatmap(df: pd.DataFrame, species: List[str], output_path: Path) -> None:
    if not species:
        warnings.warn("No species available for heatmap. Skipping heatmap.")
        return

    rclr_map = {clean_species_name(c): c for c in df.columns if c.startswith("rCLR_")}
    selected_cols = [rclr_map[s] for s in species if s in rclr_map]
    selected_species = [clean_species_name(c) for c in selected_cols]

    if not selected_cols:
        warnings.warn("No selected rCLR columns available for heatmap. Skipping heatmap.")
        return

    temp = df[["condition", "sex"] + selected_cols].copy()
    temp["condition_sex"] = temp["condition"].astype(str) + " | " + temp["sex"].astype(str)

    matrix = temp.groupby("condition_sex")[selected_cols].mean().T
    matrix.index = selected_species

    # Z-score each species across condition-sex groups for visualization.
    mat = matrix.to_numpy(dtype=float)
    row_mean = np.nanmean(mat, axis=1, keepdims=True)
    row_std = np.nanstd(mat, axis=1, keepdims=True)
    row_std[row_std == 0] = 1.0
    mat_z = (mat - row_mean) / row_std

    fig_height = max(6, 0.35 * len(selected_species))
    fig_width = max(8, 0.65 * matrix.shape[1])
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    im = ax.imshow(mat_z, aspect="auto")

    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(selected_species)))
    ax.set_yticklabels(selected_species, fontsize=8)

    ax.set_title("Top differential species: mean rCLR by condition and sex")
    ax.set_xlabel("Condition | Sex")
    ax.set_ylabel("Species")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Row-scaled mean rCLR")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def hedges_g(x1: np.ndarray, x0: np.ndarray) -> Tuple[float, float]:
    """Return Hedges' g and approximate standard error for two independent groups."""
    n1, n0 = len(x1), len(x0)
    if n1 < 2 or n0 < 2:
        return np.nan, np.nan

    s1 = np.var(x1, ddof=1)
    s0 = np.var(x0, ddof=1)
    pooled_var = ((n1 - 1) * s1 + (n0 - 1) * s0) / (n1 + n0 - 2)
    if pooled_var <= 0 or pd.isna(pooled_var):
        return np.nan, np.nan

    d = (np.mean(x1) - np.mean(x0)) / np.sqrt(pooled_var)
    correction = 1 - (3 / (4 * (n1 + n0) - 9))
    g = correction * d

    # Approximate variance for standardized mean difference.
    var_g = ((n1 + n0) / (n1 * n0)) + (g ** 2 / (2 * (n1 + n0 - 2)))
    se_g = np.sqrt(max(var_g, EPS))
    return float(g), float(se_g)


def plot_forest_top_species(df: pd.DataFrame, top_species: List[str], output_path: Path) -> None:
    """Create forest-style plot for top CRC-vs-Control rCLR species across studies."""
    rclr_map = {clean_species_name(c): c for c in df.columns if c.startswith("rCLR_")}
    selected_species = [s for s in top_species if s in rclr_map][:8]

    records = []
    sub = df[df["condition"].isin(["CRC", "Control"])].copy()

    for sp in selected_species:
        col = rclr_map[sp]
        for study, sdf in sub.groupby("study"):
            crc = pd.to_numeric(sdf.loc[sdf["condition"] == "CRC", col], errors="coerce").dropna().to_numpy()
            ctrl = pd.to_numeric(sdf.loc[sdf["condition"] == "Control", col], errors="coerce").dropna().to_numpy()
            if len(crc) < MIN_GROUP_N or len(ctrl) < MIN_GROUP_N:
                continue
            g, se = hedges_g(crc, ctrl)
            if pd.isna(g) or pd.isna(se):
                continue
            records.append(
                {
                    "species": sp,
                    "study": study,
                    "hedges_g": g,
                    "se": se,
                    "ci_low": g - 1.96 * se,
                    "ci_high": g + 1.96 * se,
                    "n_crc": len(crc),
                    "n_control": len(ctrl),
                }
            )

    eff = pd.DataFrame(records)
    if eff.empty:
        warnings.warn("No forest plot data available. Skipping forest plot.")
        return

    # Keep species with enough contributing studies.
    counts = eff.groupby("species")["study"].nunique()
    keep_species = counts[counts >= MIN_STUDIES_FOR_FOREST].index.tolist()
    eff = eff[eff["species"].isin(keep_species)]

    if eff.empty:
        warnings.warn("No top species had enough contributing studies for forest plot. Skipping forest plot.")
        return

    # Limit to first few species for readability.
    keep_species = keep_species[:6]
    eff = eff[eff["species"].isin(keep_species)].copy()

    labels = []
    y_positions = []
    y = 0
    fig_height = max(6, 0.35 * len(eff) + 1.5 * len(keep_species))
    fig, ax = plt.subplots(figsize=(10, fig_height), constrained_layout=True)

    for sp in keep_species:
        sdf = eff[eff["species"] == sp].sort_values("hedges_g")
        labels.append(sp)
        y_positions.append(y)
        ax.text(ax.get_xlim()[0] if ax.get_xlim()[0] != 0 else -3, y, sp, fontsize=9, fontweight="bold")
        y += 1
        for _, row in sdf.iterrows():
            ax.plot([row["ci_low"], row["ci_high"]], [y, y], linewidth=1)
            ax.scatter(row["hedges_g"], y, s=25)
            labels.append(f"  {row['study']}  n={row['n_crc']}/{row['n_control']}")
            y_positions.append(y)
            y += 1
        y += 1

    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Per-study Hedges' g: CRC minus Control")
    ax.set_title("Forest-style plot for top CRC-vs-Control species")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_top_species_table(rclr_results: pd.DataFrame, pa_results: pd.DataFrame) -> pd.DataFrame:
    r = rclr_results.copy()
    p = pa_results.copy()

    r["ranking_score"] = r["q_value"]
    p["ranking_score"] = p["q_value"]

    shared_cols = [
        "comparison",
        "feature_type",
        "feature",
        "species",
        "positive_group",
        "negative_group",
        "effect_name",
        "effect",
        "p_value",
        "q_value",
        "significant_fdr_0_05",
    ]

    out = pd.concat([r[shared_cols], p[shared_cols]], ignore_index=True)
    out = out.sort_values(["q_value", "comparison", "feature_type"], ascending=[True, True, True])
    return out


# =============================================================================
# Report
# =============================================================================

def write_report(
    output_path: Path,
    df: pd.DataFrame,
    rclr_cols: List[str],
    pa_cols: List[str],
    rclr_results: pd.DataFrame,
    pa_results: pd.DataFrame,
    top_species_table: pd.DataFrame,
) -> None:
    lines = []
    lines.append("STEP 3.4 — VISUALIZATION")
    lines.append("=" * 80)
    lines.append("")
    lines.append("Purpose")
    lines.append("-------")
    lines.append("This script creates visualization-oriented summaries for Phase 3 differential abundance analysis.")
    lines.append("It is self-contained and uses only X_explore.csv and metadata.csv.")
    lines.append("")
    lines.append("Important note")
    lines.append("--------------")
    lines.append("The volcano plots in this script use simple visualization statistics:")
    lines.append("  - rCLR: mean rCLR difference with Mann-Whitney p-values")
    lines.append("  - PA: prevalence difference with Fisher exact p-values")
    lines.append("They are useful for plotting and biological orientation, but the formal evidence should come")
    lines.append("from the blocked Wilcoxon, CMH, and random-effects meta-analysis scripts.")
    lines.append("")

    lines.append("Input summary")
    lines.append("-------------")
    lines.append(f"Samples after merge: {len(df)}")
    lines.append(f"rCLR features: {len(rclr_cols)}")
    lines.append(f"PA features: {len(pa_cols)}")
    lines.append("")

    lines.append("Condition counts")
    lines.append("----------------")
    for condition, n in df["condition"].value_counts().items():
        lines.append(f"{condition}: {n}")
    lines.append("")

    lines.append("Sex counts")
    lines.append("----------")
    for sex, n in df["sex"].value_counts().items():
        lines.append(f"{sex}: {n}")
    lines.append("")

    lines.append("Study counts")
    lines.append("------------")
    for study, n in df["study"].value_counts().items():
        lines.append(f"{study}: {n}")
    lines.append("")

    lines.append("Visualization summary by comparison")
    lines.append("-----------------------------------")
    for comp, _, _ in CONDITION_COMPARISONS:
        rsub = rclr_results[rclr_results["comparison"] == comp]
        psub = pa_results[pa_results["comparison"] == comp]
        lines.append(f"{comp}:")
        lines.append(f"  rCLR significant at FDR < {ALPHA_FDR}: {int(rsub['significant_fdr_0_05'].sum())}")
        lines.append(f"  PA significant at FDR < {ALPHA_FDR}: {int(psub['significant_fdr_0_05'].sum())}")

        top_r = rsub.dropna(subset=["q_value"]).sort_values("q_value").head(5)
        top_p = psub.dropna(subset=["q_value"]).sort_values("q_value").head(5)

        lines.append("  Top rCLR species:")
        if top_r.empty:
            lines.append("    None available")
        else:
            for _, row in top_r.iterrows():
                lines.append(
                    f"    {row['species']} | effect={row['effect']:.4f} | "
                    f"p={row['p_value']:.3e} | q={row['q_value']:.3e}"
                )

        lines.append("  Top PA species:")
        if top_p.empty:
            lines.append("    None available")
        else:
            for _, row in top_p.iterrows():
                lines.append(
                    f"    {row['species']} | prevalence_difference={row['effect']:.4f} | "
                    f"p={row['p_value']:.3e} | q={row['q_value']:.3e}"
                )
        lines.append("")

    lines.append("Files written")
    lines.append("-------------")
    lines.append("step3_4_volcano_rclr.png")
    lines.append("step3_4_volcano_pa.png")
    lines.append("step3_4_heatmap_top_species_condition_by_sex.png")
    lines.append("step3_4_forest_top_crc_vs_control_species.png")
    lines.append("step3_4_top_differentially_abundant_species.csv")
    lines.append("step3_4_final_report.txt")
    lines.append("")

    lines.append("Recommended interpretation")
    lines.append("--------------------------")
    lines.append("Use these plots to visually summarize Phase 3 and to select candidate species for comparison")
    lines.append("with SHAP features in Phase 5. Do not use this visualization-only script as the sole")
    lines.append("basis for statistical claims; rely on the formal blocked and meta-analysis scripts for that.")
    lines.append("")

    output_path.write_text("\n".join(lines), encoding="utf-8")


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    print("Step 3.4 — Visualization")
    print(f"Reading X from:        {X_PATH}")
    print(f"Reading metadata from: {METADATA_PATH}")
    print(f"Writing outputs to:    {OUTPUT_DIR}")

    _, _, df, _ = load_data()
    rclr_cols, pa_cols = get_feature_columns(df)

    print(f"Merged samples: {len(df)}")
    print(f"rCLR features:  {len(rclr_cols)}")
    print(f"PA features:    {len(pa_cols)}")

    rclr_results, pa_results = compute_all_visual_stats(df, rclr_cols, pa_cols)

    top_species_table = build_top_species_table(rclr_results, pa_results)
    top_species_table.to_csv(OUTPUT_DIR / "step3_4_top_differentially_abundant_species.csv", index=False)

    plot_volcano_grid(
        rclr_results,
        title="rCLR volcano summaries",
        x_label="Mean rCLR difference",
        output_path=OUTPUT_DIR / "step3_4_volcano_rclr.png",
    )

    plot_volcano_grid(
        pa_results,
        title="Presence/absence volcano summaries",
        x_label="Prevalence difference",
        output_path=OUTPUT_DIR / "step3_4_volcano_pa.png",
    )

    top_species = select_top_species_for_heatmap(rclr_results, pa_results, top_n=TOP_N_SPECIES)
    plot_condition_sex_heatmap(
        df,
        top_species,
        output_path=OUTPUT_DIR / "step3_4_heatmap_top_species_condition_by_sex.png",
    )

    # Forest-style plot is focused on CRC vs Control because it is usually the main DA contrast.
    crc_control_top = (
        rclr_results[rclr_results["comparison"] == "CRC_vs_Control"]
        .dropna(subset=["q_value"])
        .sort_values("q_value")
        ["species"]
        .drop_duplicates()
        .head(12)
        .tolist()
    )
    plot_forest_top_species(
        df,
        crc_control_top,
        output_path=OUTPUT_DIR / "step3_4_forest_top_crc_vs_control_species.png",
    )

    write_report(
        output_path=OUTPUT_DIR / "step3_4_final_report.txt",
        df=df,
        rclr_cols=rclr_cols,
        pa_cols=pa_cols,
        rclr_results=rclr_results,
        pa_results=pa_results,
        top_species_table=top_species_table,
    )

    print("Done.")
    print("Written files:")
    for path in sorted(OUTPUT_DIR.glob("step3_4_*")):
        print(f"  - {path.name}")


if __name__ == "__main__":
    main()
