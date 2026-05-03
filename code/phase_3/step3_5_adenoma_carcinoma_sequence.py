#!/usr/bin/env python3
"""
Step 3.5 — Adenoma-Carcinoma Sequence Analysis

Purpose
-------
Analyze whether species-level microbiome signals follow an ordered
Control -> Adenoma -> CRC pattern.

This script is self-contained. It reads only:
    - X_explore.csv
    - metadata.csv

It produces:
    - Venn-like overlap plot of biomarkers
    - Effect-size gradient plot for top ordered species
    - Heatmap of mean rCLR abundance across Control/Adenoma/CRC
    - CSV with trend/effect-size results
    - CSV with biomarker categories
    - Final text report

Run from:
    project/code/phase_3/

Expected project tree:
    project/
    ├── code/
    │   └── phase_3/
    │       └── step3_5_adenoma_carcinoma_sequence.py
    ├── data/
    │   └── crc_study_final_data/
    │       └── species_level/
    │           ├── X_explore.csv
    │           └── metadata.csv
    └── OUTPUTS/
        └── phase_3/
            └── 3.5/

Notes
-----
The Jonckheere-Terpstra test is implemented directly in Python.
The ordered alternative tested is:
    Control < Adenoma < CRC

Because rCLR values can be positive or negative, this tests whether a species
increases along the adenoma-carcinoma sequence. The script also reports species
that decrease along the sequence by running the same test on -abundance.
"""

from __future__ import annotations

from pathlib import Path
from itertools import combinations
from math import erf, sqrt
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy import stats
except ImportError as exc:
    raise SystemExit("This script requires scipy. Install with: pip install scipy") from exc

try:
    from statsmodels.stats.multitest import multipletests
except ImportError as exc:
    raise SystemExit("This script requires statsmodels. Install with: pip install statsmodels") from exc


# =============================================================================
# Paths
# =============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
CODE_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = CODE_ROOT.parent

DATA_DIR = (PROJECT_ROOT / "data" / "crc_study_final_data" / "species_level").resolve()
OUTPUT_DIR = (PROJECT_ROOT / "OUTPUTS" / "phase_3" / "3.5").resolve()

X_PATH = (DATA_DIR / "X_explore.csv").resolve()
METADATA_PATH = (DATA_DIR / "metadata.csv").resolve()

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Configuration
# =============================================================================

CONDITION_ORDER = ["Control", "Adenoma", "CRC"]
CONDITION_TO_SCORE = {"Control": 0, "Adenoma": 1, "CRC": 2}
ALPHA = 0.05
MIN_GROUP_N = 5
TOP_N_PLOTS = 25


# =============================================================================
# Utility functions
# =============================================================================

def normal_sf(z: float) -> float:
    """Survival function for standard normal using erf only."""
    return 0.5 * (1.0 - erf(z / sqrt(2.0)))


def clean_species_name(feature: str) -> str:
    """Remove rCLR_/PA_ prefix for readable species names."""
    if feature.startswith("rCLR_"):
        return feature[5:]
    if feature.startswith("PA_"):
        return feature[3:]
    return feature


def find_column(df: pd.DataFrame, candidates: list[str], required: bool = True) -> str | None:
    """Find a column by trying common candidate names case-insensitively."""
    exact = {c: c for c in df.columns}
    lowered = {c.lower(): c for c in df.columns}

    for cand in candidates:
        if cand in exact:
            return exact[cand]
        if cand.lower() in lowered:
            return lowered[cand.lower()]

    if required:
        raise ValueError(
            f"Could not find any of these columns: {candidates}\n"
            f"Available columns: {list(df.columns)}"
        )
    return None


def load_and_merge() -> tuple[pd.DataFrame, list[str], str, str, str | None]:
    """Load X_explore and metadata, then merge by Sample."""
    if not X_PATH.exists():
        raise FileNotFoundError(f"Missing X file: {X_PATH}")
    if not METADATA_PATH.exists():
        raise FileNotFoundError(f"Missing metadata file: {METADATA_PATH}")

    x = pd.read_csv(X_PATH)
    meta = pd.read_csv(METADATA_PATH)

    sample_x = find_column(x, ["Sample", "sample", "sample_id", "SampleID", "sampleID"])
    sample_meta = find_column(meta, ["Sample", "sample", "sample_id", "SampleID", "sampleID"])
    condition_col = find_column(meta, ["study_condition", "condition", "Condition", "diagnosis", "group"])
    study_col = find_column(meta, ["study_name", "study", "Study", "cohort", "dataset"])
    sex_col = find_column(meta, ["sex", "gender", "Sex", "Gender"], required=False)

    if sample_x != "Sample":
        x = x.rename(columns={sample_x: "Sample"})
    if sample_meta != "Sample":
        meta = meta.rename(columns={sample_meta: "Sample"})

    rclr_cols = [c for c in x.columns if c.startswith("rCLR_")]
    if not rclr_cols:
        raise ValueError("No rCLR_ columns found in X_explore.csv. Step 3.5 needs rCLR features.")

    merged = x[["Sample"] + rclr_cols].merge(meta, on="Sample", how="inner")
    merged = merged[merged[condition_col].isin(CONDITION_ORDER)].copy()

    if merged.empty:
        raise ValueError(
            "After filtering to Control/Adenoma/CRC, no samples remained. "
            f"Check condition column: {condition_col}"
        )

    return merged, rclr_cols, condition_col, study_col, sex_col


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """
    Cliff's delta for x vs y.
    Positive means x tends to be larger than y.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan

    # Efficient enough for this dataset size and 107 species.
    greater = 0
    less = 0
    for val in x:
        greater += np.sum(val > y)
        less += np.sum(val < y)
    return (greater - less) / (len(x) * len(y))


def mann_whitney_effect_p(df: pd.DataFrame, feature: str, condition_col: str, group_hi: str, group_lo: str) -> tuple[float, float]:
    """
    Mann-Whitney p-value and Cliff's delta for group_hi vs group_lo.
    Positive delta means group_hi > group_lo.
    """
    x = df.loc[df[condition_col] == group_hi, feature].astype(float).to_numpy()
    y = df.loc[df[condition_col] == group_lo, feature].astype(float).to_numpy()
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    if len(x) < MIN_GROUP_N or len(y) < MIN_GROUP_N:
        return np.nan, np.nan
    try:
        p = stats.mannwhitneyu(x, y, alternative="two-sided").pvalue
    except ValueError:
        p = np.nan
    delta = cliffs_delta(x, y)
    return p, delta


def jonckheere_terpstra(values_by_group: list[np.ndarray], alternative: str = "increasing") -> tuple[float, float, float]:
    """
    Jonckheere-Terpstra trend test for ordered groups.

    Parameters
    ----------
    values_by_group:
        List of arrays in the hypothesized order.
    alternative:
        "increasing" tests group1 < group2 < group3.
        "decreasing" tests group1 > group2 > group3 by reversing the sign.

    Returns
    -------
    J statistic, z statistic, one-sided p-value.

    Implementation notes
    --------------------
    J counts pairwise wins in the ordered direction across all group pairs.
    Ties contribute 0.5. The variance formula is the standard no-tie large-sample
    approximation. With continuous rCLR values, ties are limited but possible.
    """
    clean_groups = []
    for arr in values_by_group:
        arr = np.asarray(arr, dtype=float)
        arr = arr[~np.isnan(arr)]
        clean_groups.append(arr)

    if any(len(arr) < MIN_GROUP_N for arr in clean_groups):
        return np.nan, np.nan, np.nan

    if alternative == "decreasing":
        clean_groups = [-arr for arr in clean_groups]
    elif alternative != "increasing":
        raise ValueError("alternative must be 'increasing' or 'decreasing'")

    ns = [len(arr) for arr in clean_groups]
    J = 0.0
    for i, j in combinations(range(len(clean_groups)), 2):
        lower = clean_groups[i]
        higher = clean_groups[j]
        for x in lower:
            J += np.sum(higher > x)
            J += 0.5 * np.sum(higher == x)

    mean_J = 0.25 * sum(ns[i] * ns[j] for i, j in combinations(range(len(ns)), 2))
    var_J = (1.0 / 72.0) * (
        sum(ns[i] * ns[j] * (ns[i] + ns[j] + 1) for i, j in combinations(range(len(ns)), 2))
    )

    if var_J <= 0:
        return J, np.nan, np.nan

    z = (J - mean_J) / sqrt(var_J)
    p_one_sided = normal_sf(z)
    return J, z, p_one_sided


def analyze_sequence(df: pd.DataFrame, rclr_cols: list[str], condition_col: str) -> pd.DataFrame:
    """Run pairwise effects and JT trend tests for every rCLR species."""
    rows = []
    for feature in rclr_cols:
        by_group = {
            cond: df.loc[df[condition_col] == cond, feature].astype(float).to_numpy()
            for cond in CONDITION_ORDER
        }

        means = {cond: np.nanmean(vals) for cond, vals in by_group.items()}
        medians = {cond: np.nanmedian(vals) for cond, vals in by_group.items()}
        ns = {cond: int(np.sum(~np.isnan(vals))) for cond, vals in by_group.items()}

        j_inc, z_inc, p_inc = jonckheere_terpstra(
            [by_group[cond] for cond in CONDITION_ORDER], alternative="increasing"
        )
        j_dec, z_dec, p_dec = jonckheere_terpstra(
            [by_group[cond] for cond in CONDITION_ORDER], alternative="decreasing"
        )

        p_crc_ctrl, delta_crc_ctrl = mann_whitney_effect_p(df, feature, condition_col, "CRC", "Control")
        p_ade_ctrl, delta_ade_ctrl = mann_whitney_effect_p(df, feature, condition_col, "Adenoma", "Control")
        p_crc_ade, delta_crc_ade = mann_whitney_effect_p(df, feature, condition_col, "CRC", "Adenoma")

        gradient_mean = means["CRC"] - means["Control"]
        adenoma_position = np.nan
        denom = means["CRC"] - means["Control"]
        if np.isfinite(denom) and abs(denom) > 1e-12:
            adenoma_position = (means["Adenoma"] - means["Control"]) / denom

        monotonic_increasing_means = means["Control"] <= means["Adenoma"] <= means["CRC"]
        monotonic_decreasing_means = means["Control"] >= means["Adenoma"] >= means["CRC"]

        rows.append(
            {
                "feature": feature,
                "species": clean_species_name(feature),
                "n_control": ns["Control"],
                "n_adenoma": ns["Adenoma"],
                "n_crc": ns["CRC"],
                "mean_control": means["Control"],
                "mean_adenoma": means["Adenoma"],
                "mean_crc": means["CRC"],
                "median_control": medians["Control"],
                "median_adenoma": medians["Adenoma"],
                "median_crc": medians["CRC"],
                "gradient_mean_crc_minus_control": gradient_mean,
                "adenoma_position_between_control_crc": adenoma_position,
                "jt_J_increasing": j_inc,
                "jt_z_increasing": z_inc,
                "jt_p_increasing": p_inc,
                "jt_J_decreasing": j_dec,
                "jt_z_decreasing": z_dec,
                "jt_p_decreasing": p_dec,
                "p_crc_vs_control": p_crc_ctrl,
                "cliffs_delta_crc_vs_control": delta_crc_ctrl,
                "p_adenoma_vs_control": p_ade_ctrl,
                "cliffs_delta_adenoma_vs_control": delta_ade_ctrl,
                "p_crc_vs_adenoma": p_crc_ade,
                "cliffs_delta_crc_vs_adenoma": delta_crc_ade,
                "monotonic_increasing_means": bool(monotonic_increasing_means),
                "monotonic_decreasing_means": bool(monotonic_decreasing_means),
            }
        )

    res = pd.DataFrame(rows)

    for p_col, q_col in [
        ("jt_p_increasing", "jt_q_increasing"),
        ("jt_p_decreasing", "jt_q_decreasing"),
        ("p_crc_vs_control", "q_crc_vs_control"),
        ("p_adenoma_vs_control", "q_adenoma_vs_control"),
        ("p_crc_vs_adenoma", "q_crc_vs_adenoma"),
    ]:
        mask = res[p_col].notna()
        res[q_col] = np.nan
        if mask.any():
            res.loc[mask, q_col] = multipletests(res.loc[mask, p_col], method="fdr_bh")[1]

    # Best ordered direction for convenience.
    res["best_trend_direction"] = np.where(
        res["jt_q_increasing"].fillna(1.0) <= res["jt_q_decreasing"].fillna(1.0),
        "increasing_Control_to_Adenoma_to_CRC",
        "decreasing_Control_to_Adenoma_to_CRC",
    )
    res["best_trend_q"] = np.minimum(
        res["jt_q_increasing"].fillna(1.0),
        res["jt_q_decreasing"].fillna(1.0),
    )

    return res


def classify_biomarkers(res: pd.DataFrame) -> pd.DataFrame:
    """
    Categorize species into early, late, progression, and ordered-trend biomarkers.

    Categories are heuristic summaries for reporting and Phase 5 comparison:
    - EARLY: Adenoma vs Control significant and CRC vs Control significant.
    - LATE: CRC vs Control significant, Adenoma vs Control not significant.
    - PROGRESSION: CRC vs Adenoma significant.
    - ORDERED_INCREASING: JT increasing significant and means monotonic increasing.
    - ORDERED_DECREASING: JT decreasing significant and means monotonic decreasing.
    """
    rows = []
    for _, row in res.iterrows():
        categories = []

        crc_ctrl_sig = row.get("q_crc_vs_control", np.nan) < ALPHA
        ade_ctrl_sig = row.get("q_adenoma_vs_control", np.nan) < ALPHA
        crc_ade_sig = row.get("q_crc_vs_adenoma", np.nan) < ALPHA
        inc_sig = (row.get("jt_q_increasing", np.nan) < ALPHA) and bool(row.get("monotonic_increasing_means", False))
        dec_sig = (row.get("jt_q_decreasing", np.nan) < ALPHA) and bool(row.get("monotonic_decreasing_means", False))

        if ade_ctrl_sig and crc_ctrl_sig:
            categories.append("EARLY")
        if crc_ctrl_sig and not ade_ctrl_sig:
            categories.append("LATE")
        if crc_ade_sig:
            categories.append("PROGRESSION")
        if inc_sig:
            categories.append("ORDERED_INCREASING")
        if dec_sig:
            categories.append("ORDERED_DECREASING")

        if not categories:
            categories = ["NOT_SIGNIFICANT_BY_FDR"]

        rows.append(
            {
                "feature": row["feature"],
                "species": row["species"],
                "categories": ";".join(categories),
                "is_early": "EARLY" in categories,
                "is_late": "LATE" in categories,
                "is_progression": "PROGRESSION" in categories,
                "is_ordered_increasing": "ORDERED_INCREASING" in categories,
                "is_ordered_decreasing": "ORDERED_DECREASING" in categories,
                "best_trend_direction": row["best_trend_direction"],
                "best_trend_q": row["best_trend_q"],
                "q_crc_vs_control": row["q_crc_vs_control"],
                "q_adenoma_vs_control": row["q_adenoma_vs_control"],
                "q_crc_vs_adenoma": row["q_crc_vs_adenoma"],
                "mean_control": row["mean_control"],
                "mean_adenoma": row["mean_adenoma"],
                "mean_crc": row["mean_crc"],
                "adenoma_position_between_control_crc": row["adenoma_position_between_control_crc"],
            }
        )
    return pd.DataFrame(rows)


# =============================================================================
# Plotting
# =============================================================================

def plot_biomarker_overlap(categories: pd.DataFrame) -> None:
    """Create a compact Venn/UpSet-like bar plot of biomarker category overlaps."""
    cat_cols = ["is_early", "is_late", "is_progression", "is_ordered_increasing", "is_ordered_decreasing"]
    label_map = {
        "is_early": "Early",
        "is_late": "Late",
        "is_progression": "Progression",
        "is_ordered_increasing": "Ordered inc.",
        "is_ordered_decreasing": "Ordered dec.",
    }

    combo_counts = {}
    for _, row in categories.iterrows():
        active = [label_map[c] for c in cat_cols if bool(row[c])]
        if not active:
            continue
        label = " + ".join(active)
        combo_counts[label] = combo_counts.get(label, 0) + 1

    if not combo_counts:
        combo_counts = {"No FDR-significant biomarker categories": 0}

    items = sorted(combo_counts.items(), key=lambda x: x[1], reverse=True)
    labels = [x[0] for x in items]
    values = [x[1] for x in items]

    fig_h = max(4, 0.45 * len(labels) + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_h))
    y = np.arange(len(labels))
    ax.barh(y, values)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Number of species")
    ax.set_title("Step 3.5 — Shared vs unique biomarker categories")
    for idx, val in enumerate(values):
        ax.text(val + 0.05, idx, str(val), va="center")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "step3_5_biomarker_overlap.png", dpi=300)
    plt.close(fig)


def plot_effect_gradient(res: pd.DataFrame) -> None:
    """Line plot of Control -> Adenoma -> CRC means for top trend species."""
    plot_df = res.copy()
    plot_df = plot_df.sort_values(["best_trend_q", "jt_q_increasing"], ascending=True).head(TOP_N_PLOTS)

    if plot_df.empty:
        return

    fig, ax = plt.subplots(figsize=(11, max(6, 0.28 * len(plot_df))))
    x = np.arange(len(CONDITION_ORDER))

    for _, row in plot_df.iterrows():
        y = [row["mean_control"], row["mean_adenoma"], row["mean_crc"]]
        label = row["species"]
        ax.plot(x, y, marker="o", linewidth=1.2, alpha=0.75, label=label)

    ax.set_xticks(x)
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_ylabel("Mean rCLR abundance")
    ax.set_title(f"Step 3.5 — Effect-size gradient for top {len(plot_df)} trend species")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "step3_5_effect_size_gradient_top_species.png", dpi=300)
    plt.close(fig)


def plot_heatmap(res: pd.DataFrame) -> None:
    """Heatmap of mean rCLR abundance for top ordered/trend species."""
    plot_df = res.sort_values("best_trend_q", ascending=True).head(TOP_N_PLOTS).copy()
    if plot_df.empty:
        return

    mat = plot_df[["mean_control", "mean_adenoma", "mean_crc"]].to_numpy(dtype=float)
    species = plot_df["species"].tolist()

    # Row-standardize for pattern visualization.
    row_mean = np.nanmean(mat, axis=1, keepdims=True)
    row_sd = np.nanstd(mat, axis=1, keepdims=True)
    row_sd[row_sd == 0] = 1.0
    zmat = (mat - row_mean) / row_sd

    fig, ax = plt.subplots(figsize=(7, max(6, 0.28 * len(species))))
    im = ax.imshow(zmat, aspect="auto")
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_yticks(np.arange(len(species)))
    ax.set_yticklabels(species, fontsize=7)
    ax.set_title(f"Step 3.5 — Control → Adenoma → CRC heatmap, top {len(species)} species")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Row-scaled mean rCLR")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "step3_5_heatmap_ordered_sequence_top_species.png", dpi=300)
    plt.close(fig)


def plot_adenoma_position(res: pd.DataFrame) -> None:
    """Plot adenoma position between Control and CRC for species with CRC-Control gradient."""
    plot_df = res.copy()
    plot_df = plot_df[np.isfinite(plot_df["adenoma_position_between_control_crc"])]
    plot_df = plot_df.sort_values("best_trend_q", ascending=True).head(TOP_N_PLOTS)
    if plot_df.empty:
        return

    fig, ax = plt.subplots(figsize=(10, max(6, 0.28 * len(plot_df))))
    y = np.arange(len(plot_df))
    x = plot_df["adenoma_position_between_control_crc"].to_numpy(dtype=float)
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.axvline(0.5, linestyle="--", linewidth=1)
    ax.axvline(1, linestyle="--", linewidth=1)
    ax.scatter(x, y)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["species"].tolist(), fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Adenoma position: 0=Control-like, 1=CRC-like")
    ax.set_title("Step 3.5 — Adenoma position along Control-to-CRC mean gradient")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "step3_5_adenoma_position_top_species.png", dpi=300)
    plt.close(fig)


# =============================================================================
# Reporting
# =============================================================================

def write_report(
    df: pd.DataFrame,
    rclr_cols: list[str],
    condition_col: str,
    study_col: str,
    sex_col: str | None,
    res: pd.DataFrame,
    categories: pd.DataFrame,
) -> None:
    """Write final text report."""
    counts = df[condition_col].value_counts().reindex(CONDITION_ORDER).fillna(0).astype(int)
    n_studies = df[study_col].nunique()
    adenoma_studies = df.loc[df[condition_col] == "Adenoma", study_col].nunique()

    inc_sig = res[(res["jt_q_increasing"] < ALPHA) & (res["monotonic_increasing_means"])]
    dec_sig = res[(res["jt_q_decreasing"] < ALPHA) & (res["monotonic_decreasing_means"])]

    early = categories[categories["is_early"]]
    late = categories[categories["is_late"]]
    progression = categories[categories["is_progression"]]
    ordered_inc = categories[categories["is_ordered_increasing"]]
    ordered_dec = categories[categories["is_ordered_decreasing"]]

    top_inc = inc_sig.sort_values("jt_q_increasing", ascending=True).head(10)
    top_dec = dec_sig.sort_values("jt_q_decreasing", ascending=True).head(10)
    top_best = res.sort_values("best_trend_q", ascending=True).head(15)

    lines = []
    lines.append("STEP 3.5 — ADENOMA-CARCINOMA SEQUENCE ANALYSIS")
    lines.append("=" * 72)
    lines.append("")
    lines.append("INPUTS")
    lines.append(f"X_explore: {X_PATH}")
    lines.append(f"metadata:  {METADATA_PATH}")
    lines.append("")
    lines.append("DATA SUMMARY")
    lines.append(f"Samples included: {len(df)}")
    lines.append(f"rCLR species tested: {len(rclr_cols)}")
    lines.append(f"Condition column: {condition_col}")
    lines.append(f"Study column: {study_col}")
    lines.append(f"Sex column: {sex_col if sex_col else 'not found / not used'}")
    lines.append(f"Number of studies: {n_studies}")
    lines.append(f"Number of studies containing Adenoma: {adenoma_studies}")
    lines.append("")
    lines.append("Condition counts:")
    for cond in CONDITION_ORDER:
        lines.append(f"  {cond}: {counts.loc[cond]}")
    lines.append("")

    lines.append("METHODS")
    lines.append("For each rCLR species, the script computed mean and median abundance in Control, Adenoma, and CRC.")
    lines.append("It then tested the ordered trend Control < Adenoma < CRC using a Jonckheere-Terpstra test.")
    lines.append("The reverse ordered trend Control > Adenoma > CRC was also tested to detect decreasing sequence markers.")
    lines.append("Pairwise Mann-Whitney tests and Cliff's delta were computed for CRC vs Control, Adenoma vs Control, and CRC vs Adenoma.")
    lines.append("All p-values were corrected with Benjamini-Hochberg FDR correction.")
    lines.append("")

    lines.append("KEY FINDINGS")
    lines.append(f"Ordered increasing markers, FDR < {ALPHA}: {len(ordered_inc)}")
    lines.append(f"Ordered decreasing markers, FDR < {ALPHA}: {len(ordered_dec)}")
    lines.append(f"EARLY markers: {len(early)}")
    lines.append(f"LATE markers: {len(late)}")
    lines.append(f"PROGRESSION markers: {len(progression)}")
    lines.append("")

    lines.append("CATEGORY DEFINITIONS")
    lines.append("EARLY: significant in Adenoma vs Control and CRC vs Control.")
    lines.append("LATE: significant in CRC vs Control but not Adenoma vs Control.")
    lines.append("PROGRESSION: significant in CRC vs Adenoma.")
    lines.append("ORDERED_INCREASING: significant Jonckheere-Terpstra trend Control < Adenoma < CRC and monotonic group means.")
    lines.append("ORDERED_DECREASING: significant Jonckheere-Terpstra trend Control > Adenoma > CRC and monotonic group means.")
    lines.append("")

    def add_table(title: str, table: pd.DataFrame, q_col: str) -> None:
        lines.append(title)
        if table.empty:
            lines.append("  None found at FDR threshold.")
            lines.append("")
            return
        for _, row in table.iterrows():
            lines.append(
                "  "
                f"{row['species']} | "
                f"q={row[q_col]:.4g} | "
                f"means: Control={row['mean_control']:.4f}, "
                f"Adenoma={row['mean_adenoma']:.4f}, "
                f"CRC={row['mean_crc']:.4f} | "
                f"adenoma_position={row['adenoma_position_between_control_crc']:.3f}"
            )
        lines.append("")

    add_table("TOP ORDERED INCREASING SPECIES", top_inc, "jt_q_increasing")
    add_table("TOP ORDERED DECREASING SPECIES", top_dec, "jt_q_decreasing")

    lines.append("TOP SPECIES BY BEST TREND Q-VALUE")
    for _, row in top_best.iterrows():
        lines.append(
            "  "
            f"{row['species']} | "
            f"direction={row['best_trend_direction']} | "
            f"best_q={row['best_trend_q']:.4g} | "
            f"means: Control={row['mean_control']:.4f}, "
            f"Adenoma={row['mean_adenoma']:.4f}, "
            f"CRC={row['mean_crc']:.4f}"
        )
    lines.append("")

    lines.append("OUTPUT FILES")
    lines.append("  step3_5_sequence_analysis_all_results.csv")
    lines.append("  step3_5_biomarker_categories.csv")
    lines.append("  step3_5_biomarker_overlap.png")
    lines.append("  step3_5_effect_size_gradient_top_species.png")
    lines.append("  step3_5_heatmap_ordered_sequence_top_species.png")
    lines.append("  step3_5_adenoma_position_top_species.png")
    lines.append("  step3_5_final_report.txt")
    lines.append("")

    lines.append("INTERPRETATION GUIDANCE")
    lines.append("If many ORDERED_INCREASING or ORDERED_DECREASING markers are found, this supports a gradual microbiome shift across Control -> Adenoma -> CRC.")
    lines.append("If LATE markers dominate, the microbiome signal appears stronger at the carcinoma stage than at the adenoma stage.")
    lines.append("If EARLY markers dominate, adenoma already carries CRC-relevant microbial changes.")
    lines.append("If PROGRESSION markers are rare, CRC and Adenoma may be hard to separate at species-level rCLR abundance.")
    lines.append("Because Adenoma samples usually occur in fewer cohorts, interpret sequence findings as exploratory unless replicated by study-blocked or LODO analyses.")
    lines.append("")

    (OUTPUT_DIR / "step3_5_final_report.txt").write_text("\n".join(lines), encoding="utf-8")


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    df, rclr_cols, condition_col, study_col, sex_col = load_and_merge()

    results = analyze_sequence(df, rclr_cols, condition_col)
    categories = classify_biomarkers(results)

    results = results.sort_values(["best_trend_q", "species"], ascending=[True, True])
    categories = categories.sort_values(["best_trend_q", "species"], ascending=[True, True])

    results.to_csv(OUTPUT_DIR / "step3_5_sequence_analysis_all_results.csv", index=False)
    categories.to_csv(OUTPUT_DIR / "step3_5_biomarker_categories.csv", index=False)

    plot_biomarker_overlap(categories)
    plot_effect_gradient(results)
    plot_heatmap(results)
    plot_adenoma_position(results)

    write_report(df, rclr_cols, condition_col, study_col, sex_col, results, categories)

    print("Step 3.5 complete.")
    print(f"Outputs written to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
