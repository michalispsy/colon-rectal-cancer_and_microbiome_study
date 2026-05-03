#!/usr/bin/env python3
"""
Step 2.3 — PCA centroid distances, minimal reproducible script

Run location
------------
Place this file inside your project's code/ directory and run it from there:

    cd /absolute/path/to/project/code
    python step23_centroid_distances_minimal.py

Expected project layout
-----------------------
project/
├── code/
│   └── step23_centroid_distances_minimal.py
├── data/
│   ├── X_explore.csv
│   └── metadata.csv
└── outputs/
    └── step23_centroid_distances/
        ├── step23_pca_centroids_PC12.png
        └── step23_final_report.txt

Only two files are written:
1. step23_pca_centroids_PC12.png
2. step23_final_report.txt

Purpose
-------
Quantify whether Adenoma samples are closer to CRC or Control in PCA space.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# =============================================================================
# ABSOLUTE PATHS, DERIVED FROM THE LOCATION OF THIS SCRIPT
# =============================================================================
# The script is meant to live in: project/code/
# Therefore PROJECT_ROOT is:       project/
# These .resolve() calls produce absolute paths, so the script does not depend on
# the current working directory.

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

DATA_DIR = (PROJECT_ROOT / "data" / "crc_study_final_data" / "species_level").resolve()
OUTPUT_DIR = (PROJECT_ROOT / "OUTPUTS" / "phase_2" / "2.3").resolve()

X_PATH = (DATA_DIR / "X_explore.csv").resolve()
METADATA_PATH = (DATA_DIR / "metadata.csv").resolve()

FIGURE_PATH = (OUTPUT_DIR / "step2_3_pca_centroids_PC12.png").resolve()
REPORT_PATH = (OUTPUT_DIR / "step2_3_report.txt").resolve()


# =============================================================================
# CONFIGURATION
# =============================================================================

SAMPLE_COL = "Sample"
CONDITION_COL = "Condition"
STUDY_COL = "Study"

CONDITION_ORDER = ["Control", "Adenoma", "CRC"]
RANDOM_STATE = 42
N_BOOTSTRAP = 1000
MIN_PER_CONDITION_PER_STUDY = 3


# =============================================================================
# FUNCTIONS
# =============================================================================

def load_and_merge() -> tuple[pd.DataFrame, list[str]]:
    """Load X_explore and metadata, merge by Sample, and return feature columns."""
    if not X_PATH.exists():
        raise FileNotFoundError(f"Cannot find X_explore.csv at: {X_PATH}")
    if not METADATA_PATH.exists():
        raise FileNotFoundError(f"Cannot find metadata.csv at: {METADATA_PATH}")

    x = pd.read_csv(X_PATH)
    meta = pd.read_csv(METADATA_PATH)

    required_x = {SAMPLE_COL}
    required_meta = {SAMPLE_COL, CONDITION_COL, STUDY_COL}

    missing_x = required_x - set(x.columns)
    missing_meta = required_meta - set(meta.columns)

    if missing_x:
        raise ValueError(f"Missing columns in X_explore.csv: {sorted(missing_x)}")
    if missing_meta:
        raise ValueError(f"Missing columns in metadata.csv: {sorted(missing_meta)}")

    if x[SAMPLE_COL].duplicated().any():
        duplicates = x.loc[x[SAMPLE_COL].duplicated(), SAMPLE_COL].head().tolist()
        raise ValueError(f"Duplicate Sample IDs in X_explore.csv. Examples: {duplicates}")

    if meta[SAMPLE_COL].duplicated().any():
        duplicates = meta.loc[meta[SAMPLE_COL].duplicated(), SAMPLE_COL].head().tolist()
        raise ValueError(f"Duplicate Sample IDs in metadata.csv. Examples: {duplicates}")

    merged = x.merge(meta, on=SAMPLE_COL, how="inner")
    if merged.empty:
        raise ValueError("Merge produced zero rows. Check whether Sample IDs match.")

    dropped_x = x.shape[0] - merged.shape[0]
    dropped_meta = meta.shape[0] - merged.shape[0]
    if dropped_x or dropped_meta:
        warnings.warn(
            f"Merge did not keep all rows: dropped {dropped_x} X rows and "
            f"{dropped_meta} metadata rows."
        )

    feature_cols = [c for c in x.columns if c != SAMPLE_COL]
    non_numeric = [c for c in feature_cols if not pd.api.types.is_numeric_dtype(x[c])]
    if non_numeric:
        warnings.warn(f"Dropping non-numeric feature columns: {non_numeric[:10]}")
        feature_cols = [c for c in feature_cols if c not in non_numeric]

    before = merged.shape[0]
    merged = merged.dropna(subset=[CONDITION_COL, STUDY_COL]).copy()
    after = merged.shape[0]
    if before != after:
        warnings.warn(f"Dropped {before - after} rows with missing condition/study.")

    n_missing_features = int(merged[feature_cols].isna().sum().sum())
    if n_missing_features > 0:
        warnings.warn(f"Found {n_missing_features} missing feature values. Filling them with 0.")
        merged[feature_cols] = merged[feature_cols].fillna(0)

    available_conditions = set(merged[CONDITION_COL].astype(str))
    missing_conditions = [c for c in CONDITION_ORDER if c not in available_conditions]
    if missing_conditions:
        raise ValueError(
            f"Missing required conditions: {missing_conditions}. "
            f"Available conditions: {sorted(available_conditions)}"
        )

    return merged, feature_cols


def run_pca(df: pd.DataFrame, feature_cols: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Z-score feature matrix and run PCA."""
    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(df[feature_cols].to_numpy())

    n_components = min(x_scaled.shape[0], x_scaled.shape[1])
    pca = PCA(n_components=n_components, random_state=RANDOM_STATE)
    pcs = pca.fit_transform(x_scaled)

    pc_cols = [f"PC{i + 1}" for i in range(pcs.shape[1])]

    scores = pd.DataFrame(pcs, columns=pc_cols)
    scores.insert(0, SAMPLE_COL, df[SAMPLE_COL].values)
    scores[CONDITION_COL] = df[CONDITION_COL].values
    scores[STUDY_COL] = df[STUDY_COL].values

    explained = pd.DataFrame({
        "PC": pc_cols,
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "cumulative_explained_variance": np.cumsum(pca.explained_variance_ratio_),
    })

    return scores, explained


def pcs_for_variance(explained: pd.DataFrame, threshold: float) -> list[str]:
    """Return PC columns needed to reach a cumulative variance threshold."""
    idx = int(np.searchsorted(
        explained["cumulative_explained_variance"].to_numpy(),
        threshold,
        side="left",
    ))
    return [f"PC{i + 1}" for i in range(idx + 1)]


def condition_centroids(scores: pd.DataFrame, pc_cols: list[str]) -> pd.DataFrame:
    """Compute one PCA centroid per condition."""
    return scores.groupby(CONDITION_COL, observed=True)[pc_cols].mean().reset_index()


def euclidean(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b))


def centroid_distances(scores: pd.DataFrame, pc_cols: list[str], space_label: str) -> dict[str, float | int | str]:
    """Compute centroid distances and CRC-likeness in a given PCA space."""
    centroids = condition_centroids(scores, pc_cols)
    centroid_vectors = {
        row[CONDITION_COL]: row[pc_cols].to_numpy(dtype=float)
        for _, row in centroids.iterrows()
    }

    d_ade_crc = euclidean(centroid_vectors["Adenoma"], centroid_vectors["CRC"])
    d_ade_ctrl = euclidean(centroid_vectors["Adenoma"], centroid_vectors["Control"])
    d_crc_ctrl = euclidean(centroid_vectors["CRC"], centroid_vectors["Control"])
    crc_likeness = d_ade_ctrl / (d_ade_ctrl + d_ade_crc)

    if np.isclose(d_ade_crc, d_ade_ctrl):
        closer = "Tie"
    elif d_ade_crc < d_ade_ctrl:
        closer = "CRC"
    else:
        closer = "Control"

    return {
        "space": space_label,
        "n_pcs": len(pc_cols),
        "d_Adenoma_CRC": d_ade_crc,
        "d_Adenoma_Control": d_ade_ctrl,
        "d_CRC_Control": d_crc_ctrl,
        "CRC_likeness": crc_likeness,
        "Adenoma_closer_to": closer,
    }


def bootstrap_pc80(scores: pd.DataFrame, pc80_cols: list[str]) -> pd.DataFrame:
    """Bootstrap PC80 centroid distances by resampling within condition."""
    rng = np.random.default_rng(RANDOM_STATE)

    arrays = {
        condition: scores.loc[scores[CONDITION_COL] == condition, pc80_cols].to_numpy()
        for condition in CONDITION_ORDER
    }

    rows = []
    for _ in range(N_BOOTSTRAP):
        sampled = {}
        for condition, arr in arrays.items():
            idx = rng.integers(0, arr.shape[0], size=arr.shape[0])
            sampled[condition] = arr[idx].mean(axis=0)

        d_ade_crc = euclidean(sampled["Adenoma"], sampled["CRC"])
        d_ade_ctrl = euclidean(sampled["Adenoma"], sampled["Control"])
        d_crc_ctrl = euclidean(sampled["CRC"], sampled["Control"])
        crc_likeness = d_ade_ctrl / (d_ade_ctrl + d_ade_crc)

        rows.extend([
            {"metric": "d_Adenoma_CRC", "value": d_ade_crc},
            {"metric": "d_Adenoma_Control", "value": d_ade_ctrl},
            {"metric": "d_CRC_Control", "value": d_crc_ctrl},
            {"metric": "CRC_likeness", "value": crc_likeness},
        ])

    boot = pd.DataFrame(rows)
    return (
        boot.groupby("metric", observed=True)["value"]
        .agg(
            mean="mean",
            median="median",
            sd="std",
            ci_low=lambda x: np.quantile(x, 0.025),
            ci_high=lambda x: np.quantile(x, 0.975),
        )
        .reset_index()
    )


def per_study_pc80(scores: pd.DataFrame, pc80_cols: list[str]) -> pd.DataFrame:
    """Repeat centroid distance analysis within studies that contain all three groups."""
    rows = []

    for study, sdf in scores.groupby(STUDY_COL, observed=True):
        counts = sdf[CONDITION_COL].value_counts().to_dict()

        if not all(condition in counts for condition in CONDITION_ORDER):
            continue

        smallest_group = min(counts["Control"], counts["Adenoma"], counts["CRC"])
        if smallest_group < MIN_PER_CONDITION_PER_STUDY:
            continue

        result = centroid_distances(sdf, pc80_cols, "PC80")
        rows.append({
            "Study": study,
            "n_Control": counts["Control"],
            "n_Adenoma": counts["Adenoma"],
            "n_CRC": counts["CRC"],
            "d_Adenoma_CRC": result["d_Adenoma_CRC"],
            "d_Adenoma_Control": result["d_Adenoma_Control"],
            "d_CRC_Control": result["d_CRC_Control"],
            "CRC_likeness": result["CRC_likeness"],
            "Adenoma_closer_to": result["Adenoma_closer_to"],
        })

    return pd.DataFrame(rows).sort_values("Study") if rows else pd.DataFrame()


def plot_pc12(scores: pd.DataFrame, explained: pd.DataFrame) -> None:
    """Create the only figure: samples and condition centroids in PC1-PC2."""
    fig, ax = plt.subplots(figsize=(8, 6))

    for condition in CONDITION_ORDER:
        subset = scores[scores[CONDITION_COL] == condition]
        ax.scatter(
            subset["PC1"],
            subset["PC2"],
            alpha=0.35,
            s=18,
            label=f"{condition} samples",
        )

    centroids = condition_centroids(scores, ["PC1", "PC2"])
    for _, row in centroids.iterrows():
        ax.scatter(
            row["PC1"],
            row["PC2"],
            s=240,
            marker="X",
            edgecolor="black",
            linewidth=1.2,
            label=f"{row[CONDITION_COL]} centroid",
        )
        ax.text(
            row["PC1"],
            row["PC2"],
            f"  {row[CONDITION_COL]}",
            fontsize=10,
            weight="bold",
        )

    pc1_var = explained.loc[0, "explained_variance_ratio"] * 100
    pc2_var = explained.loc[1, "explained_variance_ratio"] * 100

    ax.set_xlabel(f"PC1 ({pc1_var:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pc2_var:.1f}% variance)")
    ax.set_title("Step 2.3: PCA centroids by condition")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURE_PATH, dpi=300)
    plt.close(fig)


def format_table(df: pd.DataFrame, float_digits: int = 4) -> str:
    """Format a pandas table for a plain text report."""
    if df.empty:
        return "No rows."
    return df.to_string(index=False, float_format=lambda x: f"{x:.{float_digits}f}")


def write_report(
    merged: pd.DataFrame,
    feature_cols: list[str],
    explained: pd.DataFrame,
    global_distances: pd.DataFrame,
    boot_summary: pd.DataFrame,
    per_study: pd.DataFrame,
) -> None:
    """Write every numeric result into one final text file."""
    pc80 = global_distances[global_distances["space"] == "PC80"].iloc[0]

    if pc80["CRC_likeness"] > 0.5:
        conclusion = "Adenoma is closer to CRC than to Control in PC80 centroid space."
    elif pc80["CRC_likeness"] < 0.5:
        conclusion = "Adenoma is closer to Control than to CRC in PC80 centroid space."
    else:
        conclusion = "Adenoma is exactly midway between CRC and Control in PC80 centroid space."

    condition_counts = merged[CONDITION_COL].value_counts().rename_axis("condition").reset_index(name="n")
    study_counts = merged[STUDY_COL].value_counts().rename_axis("study").reset_index(name="n")

    n_pc80 = int(global_distances.loc[global_distances["space"] == "PC80", "n_pcs"].iloc[0])
    n_pc90 = int(global_distances.loc[global_distances["space"] == "PC90", "n_pcs"].iloc[0])

    variance_lines = [
        f"PC1 explained variance: {explained.loc[0, 'explained_variance_ratio']:.6f}",
        f"PC2 explained variance: {explained.loc[1, 'explained_variance_ratio']:.6f}",
        f"PC1 + PC2 cumulative variance: {explained.loc[1, 'cumulative_explained_variance']:.6f}",
        f"Number of PCs for >=80% variance: {n_pc80}",
        f"Cumulative variance at PC80 cutoff: {explained.loc[n_pc80 - 1, 'cumulative_explained_variance']:.6f}",
        f"Number of PCs for >=90% variance: {n_pc90}",
        f"Cumulative variance at PC90 cutoff: {explained.loc[n_pc90 - 1, 'cumulative_explained_variance']:.6f}",
    ]

    report = f"""
STEP 2.3 — PCA CENTROID DISTANCES
==================================

Paths
-----
Project root: {PROJECT_ROOT}
Input X_explore: {X_PATH}
Input metadata: {METADATA_PATH}
Output figure: {FIGURE_PATH}
Output report: {REPORT_PATH}

Input sanity check
------------------
Number of merged samples: {merged.shape[0]}
Number of feature columns used: {len(feature_cols)}
Number of studies: {merged[STUDY_COL].nunique()}

Condition counts:
{format_table(condition_counts)}

Study counts:
{format_table(study_counts)}

PCA variance
------------
{chr(10).join(variance_lines)}

Global centroid distances
-------------------------
CRC-likeness definition:
    d(Adenoma, Control) / [d(Adenoma, Control) + d(Adenoma, CRC)]

Interpretation of CRC-likeness:
    < 0.5  Adenoma is closer to Control
    = 0.5  Adenoma is exactly intermediate
    > 0.5  Adenoma is closer to CRC

{format_table(global_distances)}

Bootstrap uncertainty for PC80 distances
----------------------------------------
Bootstrap iterations: {N_BOOTSTRAP}
Bootstrap method: resample samples within each condition, recompute condition centroids, recompute distances.

{format_table(boot_summary)}

Per-study PC80 centroid distances
---------------------------------
Only studies containing Control, Adenoma, and CRC are included.
Minimum required samples per condition per study: {MIN_PER_CONDITION_PER_STUDY}

{format_table(per_study)}

Main conclusion
---------------
Primary space: PC80, meaning the PCA coordinates explaining at least 80% cumulative variance.

d(Adenoma, CRC)     = {pc80['d_Adenoma_CRC']:.6f}
d(Adenoma, Control) = {pc80['d_Adenoma_Control']:.6f}
CRC-likeness        = {pc80['CRC_likeness']:.6f}

{conclusion}

Important caution
-----------------
This is exploratory geometry in PCA space, not a supervised classifier and not a formal causal or biological proof.
Because Adenoma samples are available only in a subset of studies, the global centroid result should be interpreted together with the per-study table.
""".strip()

    REPORT_PATH.write_text(report + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    merged, feature_cols = load_and_merge()
    scores, explained = run_pca(merged, feature_cols)

    plot_pc12(scores, explained)

    pc12_cols = ["PC1", "PC2"]
    pc80_cols = pcs_for_variance(explained, 0.80)
    pc90_cols = pcs_for_variance(explained, 0.90)

    global_distances = pd.DataFrame([
        centroid_distances(scores, pc12_cols, "PC1_PC2"),
        centroid_distances(scores, pc80_cols, "PC80"),
        centroid_distances(scores, pc90_cols, "PC90"),
    ])

    boot_summary = bootstrap_pc80(scores, pc80_cols)
    per_study = per_study_pc80(scores, pc80_cols)

    write_report(
        merged=merged,
        feature_cols=feature_cols,
        explained=explained,
        global_distances=global_distances,
        boot_summary=boot_summary,
        per_study=per_study,
    )

    print("Step 2.3 complete.")
    print(f"Figure written to: {FIGURE_PATH}")
    print(f"Report written to: {REPORT_PATH}")


if __name__ == "__main__":
    main()
