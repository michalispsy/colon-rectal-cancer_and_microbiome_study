#!/usr/bin/env python3
"""
Step 3.3 — Random Effects Meta-Analysis for rCLR Species Features

Purpose
-------
Performs per-species random-effects meta-analysis across studies using
standardized mean differences (Hedges' g) computed within each study.

This is the Python analogue of:
    metafor::escalc(..., measure = "SMD")
    metafor::rma(..., method = "REML")

Inputs
------
Expected project structure:

project/
├── code/
│   └── step3_3_random_effects_meta_analysis.py
├── data/
│   └── crc_study_final_data/
│       └── species_level/
│           ├── X_explore.csv
│           └── metadata.csv
└── OUTPUTS/
    └── phase_3/
        └── 3.3/

Outputs
-------
OUTPUTS/phase_3/3.3/
├── step3_3_meta_analysis_all_results.csv
├── step3_3_meta_analysis_significant_results.csv
├── step3_3_per_study_effects.csv
├── step3_3_meta_analysis_summary.png
├── step3_3_forest_top_species.png
└── step3_3_final_report.txt

Notes
-----
- Uses only rCLR_ features from X_explore.csv.
- Uses Study as the meta-analysis unit.
- Main comparisons are disease-stage comparisons:
    1. CRC vs Control
    2. Adenoma vs Control
    3. CRC vs Adenoma
- Each comparison is run in:
    - all samples
    - Male only
    - Female only
- The positive effect direction is always Group 1 minus Group 0.
  Example: CRC vs Control means positive Hedges' g = higher in CRC.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import norm


# ---------------------------------------------------------------------
# Paths requested by the user
# ---------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
CODE_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = CODE_ROOT.parent

DATA_DIR = (PROJECT_ROOT / "data" / "crc_study_final_data" / "species_level").resolve()
OUTPUT_DIR = (PROJECT_ROOT / "OUTPUTS" / "phase_3" / "3.3").resolve()

X_PATH = (DATA_DIR / "X_explore.csv").resolve()
METADATA_PATH = (DATA_DIR / "metadata.csv").resolve()


# ---------------------------------------------------------------------
# Analysis parameters
# ---------------------------------------------------------------------
MIN_N_PER_GROUP_PER_STUDY = 3
ALPHA_FDR = 0.05
TOP_N_FOR_REPORT = 15
TOP_N_FOREST = 6

# Positive direction is group1 - group0
COMPARISONS = [
    {
        "comparison": "CRC_vs_Control",
        "group1": "CRC",
        "group0": "Control",
        "label": "CRC vs Control",
    },
    {
        "comparison": "Adenoma_vs_Control",
        "group1": "Adenoma",
        "group0": "Control",
        "label": "Adenoma vs Control",
    },
    {
        "comparison": "CRC_vs_Adenoma",
        "group1": "CRC",
        "group0": "Adenoma",
        "label": "CRC vs Adenoma",
    },
]

STRATA = [
    {"stratum": "All", "gender": None, "label": "All samples"},
    {"stratum": "Male", "gender": "Male", "label": "Males only"},
    {"stratum": "Female", "gender": "Female", "label": "Females only"},
]


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------
def ensure_output_dir() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def normalize_condition(value: object) -> str:
    text = str(value).strip()
    low = text.lower()
    if low in {"crc", "cancer", "colorectal_cancer", "case", "cases"}:
        return "CRC"
    if low in {"control", "controls", "healthy", "normal"}:
        return "Control"
    if low in {"adenoma", "adenomas", "ade"}:
        return "Adenoma"
    return text


def normalize_gender(value: object) -> str:
    text = str(value).strip()
    low = text.lower()
    if low in {"m", "male", "man", "men"}:
        return "Male"
    if low in {"f", "female", "woman", "women"}:
        return "Female"
    return text


def clean_species_name(feature: str) -> str:
    return feature.replace("rCLR_", "")


def benjamini_hochberg(pvalues: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR correction.

    Returns q-values in the original order.
    NaN p-values remain NaN.
    """
    pvalues = np.asarray(list(pvalues), dtype=float)
    qvalues = np.full_like(pvalues, np.nan, dtype=float)
    valid = np.isfinite(pvalues)
    p = pvalues[valid]
    m = len(p)
    if m == 0:
        return qvalues

    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    temp = np.empty(m, dtype=float)
    temp[order] = adjusted
    qvalues[valid] = temp
    return qvalues


@dataclass
class StudyEffect:
    study: str
    n_group1: int
    n_group0: int
    mean_group1: float
    mean_group0: float
    sd_group1: float
    sd_group0: float
    hedges_g: float
    variance: float
    se: float


@dataclass
class MetaResult:
    effect_random: float
    se_random: float
    z_value: float
    p_value: float
    ci_low: float
    ci_high: float
    tau2_reml: float
    q_statistic: float
    q_df: int
    q_p_value: float
    i2_percent: float
    n_studies: int
    total_n_group1: int
    total_n_group0: int


def hedges_g_for_two_groups(x1: np.ndarray, x0: np.ndarray) -> Optional[Tuple[float, float, float, float, float, float, float]]:
    """Compute Hedges' g and sampling variance for two independent groups.

    Returns:
        hedges_g, variance, mean1, mean0, sd1, sd0, pooled_sd
    """
    x1 = np.asarray(x1, dtype=float)
    x0 = np.asarray(x0, dtype=float)
    x1 = x1[np.isfinite(x1)]
    x0 = x0[np.isfinite(x0)]
    n1 = len(x1)
    n0 = len(x0)

    if n1 < MIN_N_PER_GROUP_PER_STUDY or n0 < MIN_N_PER_GROUP_PER_STUDY:
        return None

    mean1 = float(np.mean(x1))
    mean0 = float(np.mean(x0))
    sd1 = float(np.std(x1, ddof=1))
    sd0 = float(np.std(x0, ddof=1))

    df = n1 + n0 - 2
    if df <= 0:
        return None

    pooled_var = ((n1 - 1) * sd1**2 + (n0 - 1) * sd0**2) / df
    if not np.isfinite(pooled_var) or pooled_var <= 0:
        return None

    pooled_sd = math.sqrt(pooled_var)
    cohens_d = (mean1 - mean0) / pooled_sd

    # Hedges small-sample correction.
    # J = 1 - 3/(4df - 1)
    correction = 1.0 - (3.0 / (4.0 * df - 1.0)) if (4.0 * df - 1.0) > 0 else 1.0
    hedges_g = correction * cohens_d

    # Approximate sampling variance for Hedges' g.
    # Common SMD approximation used in meta-analysis.
    variance = ((n1 + n0) / (n1 * n0)) + ((hedges_g**2) / (2 * (n1 + n0 - 2)))
    if not np.isfinite(variance) or variance <= 0:
        return None

    return hedges_g, variance, mean1, mean0, sd1, sd0, pooled_sd


def compute_study_effects(
    data: pd.DataFrame,
    feature: str,
    group1: str,
    group0: str,
) -> List[StudyEffect]:
    effects: List[StudyEffect] = []

    for study, study_df in data.groupby("Study", sort=True):
        x1 = study_df.loc[study_df["Condition"] == group1, feature].to_numpy(dtype=float)
        x0 = study_df.loc[study_df["Condition"] == group0, feature].to_numpy(dtype=float)
        result = hedges_g_for_two_groups(x1, x0)
        if result is None:
            continue

        hedges_g, variance, mean1, mean0, sd1, sd0, _pooled_sd = result
        effects.append(
            StudyEffect(
                study=str(study),
                n_group1=int(np.isfinite(x1).sum()),
                n_group0=int(np.isfinite(x0).sum()),
                mean_group1=mean1,
                mean_group0=mean0,
                sd_group1=sd1,
                sd_group0=sd0,
                hedges_g=float(hedges_g),
                variance=float(variance),
                se=float(math.sqrt(variance)),
            )
        )

    return effects


def reml_tau2_intercept_only(yi: np.ndarray, vi: np.ndarray) -> float:
    """Estimate tau^2 by REML for an intercept-only random-effects model.

    The restricted negative log-likelihood is minimized over tau^2 >= 0.
    Constants are omitted because they do not affect the optimum.
    """
    yi = np.asarray(yi, dtype=float)
    vi = np.asarray(vi, dtype=float)
    k = len(yi)

    if k <= 1:
        return 0.0

    valid = np.isfinite(yi) & np.isfinite(vi) & (vi > 0)
    yi = yi[valid]
    vi = vi[valid]
    k = len(yi)
    if k <= 1:
        return 0.0

    def restricted_nll(tau2: float) -> float:
        v = vi + tau2
        if np.any(v <= 0):
            return np.inf
        w = 1.0 / v
        w_sum = np.sum(w)
        mu_hat = np.sum(w * yi) / w_sum
        resid = yi - mu_hat
        # REML profile likelihood for intercept-only model.
        return 0.5 * (np.sum(np.log(v)) + np.log(w_sum) + np.sum(w * resid**2))

    # Conservative upper bound. If between-study heterogeneity is large,
    # this still gives the optimizer enough room for microbiome SMDs.
    q75, q25 = np.percentile(yi, [75, 25])
    iqr = max(q75 - q25, 1e-6)
    upper = max(10.0, float(np.var(yi, ddof=1) * 10.0 + iqr**2 * 10.0))

    opt = minimize_scalar(restricted_nll, bounds=(0.0, upper), method="bounded", options={"xatol": 1e-10})

    if not opt.success or not np.isfinite(opt.x):
        return 0.0

    tau2 = max(0.0, float(opt.x))

    # If the optimum is extremely close to zero, report exactly zero.
    if tau2 < 1e-10:
        tau2 = 0.0
    return tau2


def random_effects_meta_analysis(effects: List[StudyEffect]) -> Optional[MetaResult]:
    if len(effects) < 2:
        return None

    yi = np.asarray([e.hedges_g for e in effects], dtype=float)
    vi = np.asarray([e.variance for e in effects], dtype=float)
    valid = np.isfinite(yi) & np.isfinite(vi) & (vi > 0)
    yi = yi[valid]
    vi = vi[valid]
    kept_effects = [e for e, keep in zip(effects, valid) if keep]
    k = len(yi)

    if k < 2:
        return None

    # Fixed-effect Q for heterogeneity.
    w_fixed = 1.0 / vi
    mu_fixed = np.sum(w_fixed * yi) / np.sum(w_fixed)
    q_statistic = float(np.sum(w_fixed * (yi - mu_fixed) ** 2))
    q_df = k - 1
    # Chi-square survival function without importing chi2 at top? import here.
    from scipy.stats import chi2

    q_p_value = float(chi2.sf(q_statistic, q_df)) if q_df > 0 else np.nan
    i2_percent = float(max(0.0, (q_statistic - q_df) / q_statistic) * 100.0) if q_statistic > 0 else 0.0

    tau2 = reml_tau2_intercept_only(yi, vi)
    w_random = 1.0 / (vi + tau2)
    mu_random = float(np.sum(w_random * yi) / np.sum(w_random))
    se_random = float(math.sqrt(1.0 / np.sum(w_random)))

    if se_random <= 0 or not np.isfinite(se_random):
        return None

    z_value = mu_random / se_random
    p_value = float(2.0 * norm.sf(abs(z_value)))
    ci_low = float(mu_random - 1.96 * se_random)
    ci_high = float(mu_random + 1.96 * se_random)

    return MetaResult(
        effect_random=mu_random,
        se_random=se_random,
        z_value=float(z_value),
        p_value=p_value,
        ci_low=ci_low,
        ci_high=ci_high,
        tau2_reml=tau2,
        q_statistic=q_statistic,
        q_df=q_df,
        q_p_value=q_p_value,
        i2_percent=i2_percent,
        n_studies=k,
        total_n_group1=int(sum(e.n_group1 for e in kept_effects)),
        total_n_group0=int(sum(e.n_group0 for e in kept_effects)),
    )


def load_and_merge() -> Tuple[pd.DataFrame, List[str]]:
    if not X_PATH.exists():
        raise FileNotFoundError(f"Missing X_explore file: {X_PATH}")
    if not METADATA_PATH.exists():
        raise FileNotFoundError(f"Missing metadata file: {METADATA_PATH}")

    x = pd.read_csv(X_PATH)
    meta = pd.read_csv(METADATA_PATH)

    required_meta = {"Sample", "Study", "Condition", "Gender"}
    missing_meta = required_meta - set(meta.columns)
    if missing_meta:
        raise ValueError(f"metadata.csv is missing required columns: {sorted(missing_meta)}")
    if "Sample" not in x.columns:
        raise ValueError("X_explore.csv must contain a 'Sample' column.")

    rclr_features = [c for c in x.columns if c.startswith("rCLR_")]
    if not rclr_features:
        raise ValueError("No rCLR_ features found in X_explore.csv.")

    x = x[["Sample"] + rclr_features].copy()
    meta = meta.copy()
    meta["Condition"] = meta["Condition"].map(normalize_condition)
    meta["Gender"] = meta["Gender"].map(normalize_gender)

    merged = meta.merge(x, on="Sample", how="inner")
    if merged.empty:
        raise ValueError("After merging metadata with X_explore by Sample, zero rows remain.")

    # Ensure rCLR features are numeric.
    for col in rclr_features:
        merged[col] = pd.to_numeric(merged[col], errors="coerce")

    return merged, rclr_features


def run_meta_analysis(data: pd.DataFrame, rclr_features: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    all_meta_rows: List[Dict[str, object]] = []
    all_study_rows: List[Dict[str, object]] = []

    for comp in COMPARISONS:
        comparison = comp["comparison"]
        group1 = comp["group1"]
        group0 = comp["group0"]

        for stratum in STRATA:
            stratum_name = stratum["stratum"]
            gender = stratum["gender"]

            subset = data[data["Condition"].isin([group1, group0])].copy()
            if gender is not None:
                subset = subset[subset["Gender"] == gender].copy()

            if subset.empty:
                continue

            for feature in rclr_features:
                study_effects = compute_study_effects(subset, feature, group1, group0)

                for e in study_effects:
                    all_study_rows.append(
                        {
                            "comparison": comparison,
                            "stratum": stratum_name,
                            "species": clean_species_name(feature),
                            "feature": feature,
                            "study": e.study,
                            "group1": group1,
                            "group0": group0,
                            "n_group1": e.n_group1,
                            "n_group0": e.n_group0,
                            "mean_group1": e.mean_group1,
                            "mean_group0": e.mean_group0,
                            "sd_group1": e.sd_group1,
                            "sd_group0": e.sd_group0,
                            "hedges_g": e.hedges_g,
                            "variance": e.variance,
                            "se": e.se,
                        }
                    )

                meta = random_effects_meta_analysis(study_effects)
                if meta is None:
                    all_meta_rows.append(
                        {
                            "comparison": comparison,
                            "stratum": stratum_name,
                            "species": clean_species_name(feature),
                            "feature": feature,
                            "group1": group1,
                            "group0": group0,
                            "n_studies": len(study_effects),
                            "total_n_group1": sum(e.n_group1 for e in study_effects),
                            "total_n_group0": sum(e.n_group0 for e in study_effects),
                            "effect_random": np.nan,
                            "se_random": np.nan,
                            "ci_low": np.nan,
                            "ci_high": np.nan,
                            "z_value": np.nan,
                            "p_value": np.nan,
                            "tau2_reml": np.nan,
                            "q_statistic": np.nan,
                            "q_df": np.nan,
                            "q_p_value": np.nan,
                            "i2_percent": np.nan,
                            "direction": "not_tested",
                            "note": "Fewer than two studies with valid per-study SMD.",
                        }
                    )
                    continue

                direction = (
                    f"Higher in {group1}"
                    if meta.effect_random > 0
                    else f"Higher in {group0}"
                    if meta.effect_random < 0
                    else "No direction"
                )
                all_meta_rows.append(
                    {
                        "comparison": comparison,
                        "stratum": stratum_name,
                        "species": clean_species_name(feature),
                        "feature": feature,
                        "group1": group1,
                        "group0": group0,
                        "n_studies": meta.n_studies,
                        "total_n_group1": meta.total_n_group1,
                        "total_n_group0": meta.total_n_group0,
                        "effect_random": meta.effect_random,
                        "se_random": meta.se_random,
                        "ci_low": meta.ci_low,
                        "ci_high": meta.ci_high,
                        "z_value": meta.z_value,
                        "p_value": meta.p_value,
                        "tau2_reml": meta.tau2_reml,
                        "q_statistic": meta.q_statistic,
                        "q_df": meta.q_df,
                        "q_p_value": meta.q_p_value,
                        "i2_percent": meta.i2_percent,
                        "direction": direction,
                        "note": "",
                    }
                )

    meta_df = pd.DataFrame(all_meta_rows)
    study_df = pd.DataFrame(all_study_rows)

    if not meta_df.empty:
        meta_df["q_value"] = np.nan
        for (comparison, stratum), idx in meta_df.groupby(["comparison", "stratum"]).groups.items():
            meta_df.loc[idx, "q_value"] = benjamini_hochberg(meta_df.loc[idx, "p_value"].to_numpy())
        meta_df["significant_fdr_0_05"] = meta_df["q_value"] < ALPHA_FDR
        meta_df["abs_effect_random"] = meta_df["effect_random"].abs()
        meta_df = meta_df.sort_values(
            ["comparison", "stratum", "q_value", "p_value", "abs_effect_random"],
            ascending=[True, True, True, True, False],
        )

    return meta_df, study_df


def make_summary_plot(meta_df: pd.DataFrame) -> None:
    valid = meta_df[np.isfinite(meta_df["p_value"])]
    if valid.empty:
        return

    summary = (
        valid.groupby(["comparison", "stratum"], as_index=False)
        .agg(
            n_tested=("species", "count"),
            n_significant=("significant_fdr_0_05", "sum"),
            median_abs_effect=("abs_effect_random", "median"),
            median_i2=("i2_percent", "median"),
        )
        .sort_values(["comparison", "stratum"])
    )

    labels = [f"{r.comparison}\n{r.stratum}" for r in summary.itertuples()]
    x = np.arange(len(summary))

    fig, ax1 = plt.subplots(figsize=(max(10, len(summary) * 1.35), 6))
    bars = ax1.bar(x, summary["n_significant"].to_numpy())
    ax1.set_ylabel("Significant species (FDR < 0.05)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=45, ha="right")
    ax1.set_title("Step 3.3 Random-effects meta-analysis summary")

    for bar, value in zip(bars, summary["n_significant"].to_numpy()):
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(int(value)), ha="center", va="bottom", fontsize=9)

    ax2 = ax1.twinx()
    ax2.plot(x, summary["median_i2"].to_numpy(), marker="o")
    ax2.set_ylabel("Median I² (%)")

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "step3_3_meta_analysis_summary.png", dpi=200)
    plt.close(fig)


def select_forest_species(meta_df: pd.DataFrame) -> pd.DataFrame:
    sig = meta_df[meta_df["significant_fdr_0_05"] & np.isfinite(meta_df["q_value"])].copy()
    if sig.empty:
        valid = meta_df[np.isfinite(meta_df["p_value"])].copy()
        return valid.sort_values(["q_value", "p_value"]).head(TOP_N_FOREST)
    return sig.sort_values(["q_value", "p_value", "abs_effect_random"], ascending=[True, True, False]).head(TOP_N_FOREST)


def make_forest_plot(meta_df: pd.DataFrame, study_df: pd.DataFrame) -> None:
    selected = select_forest_species(meta_df)
    if selected.empty or study_df.empty:
        return

    n_panels = len(selected)
    fig_height = max(4, 2.2 * n_panels)
    fig, axes = plt.subplots(n_panels, 1, figsize=(12, fig_height), squeeze=False)
    axes = axes[:, 0]

    for ax, row in zip(axes, selected.itertuples(index=False)):
        sub = study_df[
            (study_df["comparison"] == row.comparison)
            & (study_df["stratum"] == row.stratum)
            & (study_df["species"] == row.species)
        ].copy()
        sub = sub.sort_values("study")

        if sub.empty:
            ax.axis("off")
            continue

        y = np.arange(len(sub))
        effects = sub["hedges_g"].to_numpy(dtype=float)
        ses = sub["se"].to_numpy(dtype=float)
        ci_low = effects - 1.96 * ses
        ci_high = effects + 1.96 * ses
        xerr = np.vstack([effects - ci_low, ci_high - effects])

        ax.errorbar(effects, y, xerr=xerr, fmt="o", capsize=3, label="Study SMD")
        ax.axvline(0, linestyle="--", linewidth=1)
        ax.axvline(row.effect_random, linewidth=2, label="REML pooled SMD")
        ax.fill_betweenx(
            [-0.5, len(sub) - 0.5],
            row.ci_low,
            row.ci_high,
            alpha=0.15,
            label="Pooled 95% CI",
        )
        ax.set_yticks(y)
        ax.set_yticklabels(sub["study"].tolist(), fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Hedges' g")
        ax.set_title(
            f"{row.comparison} | {row.stratum} | {row.species}\n"
            f"pooled g={row.effect_random:.3f}, q={row.q_value:.3g}, I²={row.i2_percent:.1f}%",
            fontsize=10,
        )

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right")
    fig.tight_layout(rect=[0, 0, 0.92, 1])
    fig.savefig(OUTPUT_DIR / "step3_3_forest_top_species.png", dpi=200)
    plt.close(fig)


def write_report(data: pd.DataFrame, rclr_features: List[str], meta_df: pd.DataFrame, study_df: pd.DataFrame) -> None:
    report_path = OUTPUT_DIR / "step3_3_final_report.txt"

    lines: List[str] = []
    lines.append("STEP 3.3 — RANDOM EFFECTS META-ANALYSIS")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Input paths")
    lines.append("-" * 72)
    lines.append(f"X_explore: {X_PATH}")
    lines.append(f"metadata:  {METADATA_PATH}")
    lines.append(f"output:    {OUTPUT_DIR}")
    lines.append("")

    lines.append("Data sanity check")
    lines.append("-" * 72)
    lines.append(f"Merged samples: {len(data)}")
    lines.append(f"rCLR species features analyzed: {len(rclr_features)}")
    lines.append(f"Studies: {data['Study'].nunique()}")
    lines.append("")
    lines.append("Condition counts:")
    for key, value in data["Condition"].value_counts().sort_index().items():
        lines.append(f"  {key}: {value}")
    lines.append("")
    lines.append("Gender counts:")
    for key, value in data["Gender"].value_counts().sort_index().items():
        lines.append(f"  {key}: {value}")
    lines.append("")

    lines.append("Method")
    lines.append("-" * 72)
    lines.append("For each rCLR species and each study, the script computes Hedges' g")
    lines.append("for Group 1 minus Group 0. Hedges' g is a standardized mean difference")
    lines.append("with a small-sample correction. Per-study SMDs are then combined across")
    lines.append("studies using an intercept-only random-effects model with REML tau²")
    lines.append("estimation. The report includes pooled effect size, standard error, 95% CI,")
    lines.append("p-value, Benjamini-Hochberg q-value, Cochran's Q, Q p-value, tau², and I².")
    lines.append("")
    lines.append(f"Minimum per-study group size: {MIN_N_PER_GROUP_PER_STUDY} per group")
    lines.append("A species/comparison requires at least two valid studies to be meta-analyzed.")
    lines.append("Positive effect means higher rCLR abundance in Group 1.")
    lines.append("")

    if meta_df.empty:
        lines.append("No valid meta-analysis results were produced.")
        report_path.write_text("\n".join(lines), encoding="utf-8")
        return

    valid = meta_df[np.isfinite(meta_df["p_value"])].copy()
    sig = valid[valid["significant_fdr_0_05"]].copy()

    lines.append("Global summary by comparison and stratum")
    lines.append("-" * 72)
    if valid.empty:
        lines.append("No species had enough valid studies for meta-analysis.")
    else:
        summary = (
            valid.groupby(["comparison", "stratum"], as_index=False)
            .agg(
                n_tested=("species", "count"),
                n_significant=("significant_fdr_0_05", "sum"),
                median_abs_effect=("abs_effect_random", "median"),
                median_i2=("i2_percent", "median"),
                median_n_studies=("n_studies", "median"),
            )
            .sort_values(["comparison", "stratum"])
        )
        for r in summary.itertuples(index=False):
            lines.append(
                f"{r.comparison} | {r.stratum}: "
                f"tested={int(r.n_tested)}, significant_FDR<0.05={int(r.n_significant)}, "
                f"median|g|={r.median_abs_effect:.3f}, median_I2={r.median_i2:.1f}%, "
                f"median_valid_studies={r.median_n_studies:.1f}"
            )
    lines.append("")

    lines.append("Top results per comparison/stratum")
    lines.append("-" * 72)
    if valid.empty:
        lines.append("No valid p-values available.")
    else:
        for (comparison, stratum), group in valid.groupby(["comparison", "stratum"], sort=True):
            lines.append("")
            lines.append(f"{comparison} | {stratum}")
            lines.append("~" * 72)
            top = group.sort_values(["q_value", "p_value", "abs_effect_random"], ascending=[True, True, False]).head(TOP_N_FOR_REPORT)
            if top.empty:
                lines.append("No valid species.")
                continue
            for r in top.itertuples(index=False):
                sig_marker = "*" if bool(r.significant_fdr_0_05) else " "
                lines.append(
                    f"{sig_marker} {r.species}: "
                    f"g={r.effect_random:.3f} [{r.ci_low:.3f}, {r.ci_high:.3f}], "
                    f"p={r.p_value:.3g}, q={r.q_value:.3g}, "
                    f"I2={r.i2_percent:.1f}%, tau2={r.tau2_reml:.4f}, "
                    f"studies={int(r.n_studies)}, direction={r.direction}"
                )
    lines.append("")

    lines.append("Most heterogeneous significant species")
    lines.append("-" * 72)
    if sig.empty:
        lines.append("No FDR-significant species were found, so heterogeneity among significant species is not applicable.")
    else:
        top_hetero = sig.sort_values("i2_percent", ascending=False).head(TOP_N_FOR_REPORT)
        for r in top_hetero.itertuples(index=False):
            lines.append(
                f"{r.comparison} | {r.stratum} | {r.species}: "
                f"I2={r.i2_percent:.1f}%, Q={r.q_statistic:.3f}, Q_p={r.q_p_value:.3g}, "
                f"g={r.effect_random:.3f}, q={r.q_value:.3g}"
            )
    lines.append("")

    lines.append("Interpretation guide")
    lines.append("-" * 72)
    lines.append("- Hedges' g > 0 means higher rCLR abundance in the first group of the comparison.")
    lines.append("  Example: CRC_vs_Control with g > 0 means higher in CRC.")
    lines.append("- Hedges' g < 0 means higher in the second group.")
    lines.append("- I² describes the fraction of between-study heterogeneity not explained by sampling error.")
    lines.append("  High I² means the species effect is cohort-dependent and should be interpreted cautiously.")
    lines.append("- Sex-stratified analyses are exploratory if few studies contain enough male/female samples")
    lines.append("  in both comparison groups.")
    lines.append("")

    lines.append("Generated files")
    lines.append("-" * 72)
    lines.append("step3_3_meta_analysis_all_results.csv")
    lines.append("step3_3_meta_analysis_significant_results.csv")
    lines.append("step3_3_per_study_effects.csv")
    lines.append("step3_3_meta_analysis_summary.png")
    lines.append("step3_3_forest_top_species.png")
    lines.append("step3_3_final_report.txt")
    lines.append("")

    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    ensure_output_dir()

    data, rclr_features = load_and_merge()
    meta_df, study_df = run_meta_analysis(data, rclr_features)

    all_results_path = OUTPUT_DIR / "step3_3_meta_analysis_all_results.csv"
    significant_path = OUTPUT_DIR / "step3_3_meta_analysis_significant_results.csv"
    per_study_path = OUTPUT_DIR / "step3_3_per_study_effects.csv"

    meta_df.to_csv(all_results_path, index=False)
    if meta_df.empty:
        pd.DataFrame().to_csv(significant_path, index=False)
    else:
        meta_df[meta_df["significant_fdr_0_05"]].to_csv(significant_path, index=False)
    study_df.to_csv(per_study_path, index=False)

    make_summary_plot(meta_df)
    make_forest_plot(meta_df, study_df)
    write_report(data, rclr_features, meta_df, study_df)

    print("Step 3.3 completed successfully.")
    print(f"Outputs written to: {OUTPUT_DIR}")
    print(f"All results: {all_results_path.name}")
    print(f"Final report: step3_3_final_report.txt")


if __name__ == "__main__":
    main()
