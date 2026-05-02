"""
STEP 1–2 — DATA PREPARATION & UNSUPERVISED ANALYSIS

This script prepares the harmonized dataset and performs exploratory analysis
to understand the structure of the microbiome data before supervised modeling.

Main functionalities:
- Loads and merges:
  - CLR-transformed genus abundance matrix
  - metadata table (condition, sex, study)
- Performs sanity checks:
  - duplicate sample IDs
  - missing values
  - dataset dimensions and distributions
- Constructs binary classification label (CRC vs Control)
- Prepares feature matrix (X) for analysis
- Applies dimensionality reduction:
  - PCA (variance explanation + visualization)
  - t-SNE (nonlinear structure)
  - UMAP (optional)
- Performs clustering:
  - KMeans (k=2–6)
  - DBSCAN
- Evaluates clustering quality:
  - Silhouette score
  - Adjusted Rand Index (ARI) vs:
    - study (batch effects)
    - condition (disease signal)
    - gender (sex effects)
- Generates plots:
  - PCA / t-SNE colored by study, condition, gender
- Audits dataset composition:
  - study sizes
  - sex balance
  - condition balance

Blueprint steps covered:
- Data cleaning and harmonization
- Feature matrix preparation
- Unsupervised exploration (PCA, t-SNE, clustering)
- Detection of batch effects (study-driven structure)
- Preliminary assessment of condition and sex effects
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.preprocessing import StandardScaler

try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False

warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================
# Configuration
# ============================================================

ROOT = Path(".")
FEATURES_PATH = ROOT / "harmonized" / "unified_genera_clr.csv"
METADATA_PATH = ROOT / "harmonized" / "unified_metadata.csv"
OUTPUT_DIR = ROOT / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True)

SAMPLE_COL = "Sample"
STUDY_COL = "Study"
TARGET_COL = "Condition"
SEX_COL = "Gender"

RANDOM_STATE = 42


# ============================================================
# Utilities
# ============================================================

def save_df(df: pd.DataFrame, filename: str) -> None:
    path = OUTPUT_DIR / filename
    df.to_csv(path, index=False)
    print(f"[saved] {path}")


def print_header(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def normalize_labels(series: pd.Series) -> pd.Series:
    """Basic cleaning for labels."""
    return (
        series.astype(str)
        .str.strip()
        .str.replace(r"\s+", " ", regex=True)
    )


def infer_binary_crc_label(condition: pd.Series) -> pd.Series:
    """
    Build a binary target:
    1 = CRC/cancer-like
    0 = control/healthy-like
    Other labels become NaN and can be excluded later.
    """
    c = normalize_labels(condition).str.lower()

    positive_patterns = [
        "crc", "cancer", "carcinoma", "colorectal cancer", "tumor", "tumour"
    ]
    negative_patterns = [
        "control", "healthy", "normal"
    ]

    y = pd.Series(np.nan, index=condition.index, dtype="float")

    pos_mask = pd.Series(False, index=condition.index)
    neg_mask = pd.Series(False, index=condition.index)

    for p in positive_patterns:
        pos_mask = pos_mask | c.str.contains(p, na=False)
    for p in negative_patterns:
        neg_mask = neg_mask | c.str.contains(p, na=False)

    y.loc[pos_mask] = 1
    y.loc[neg_mask] = 0
    return y


def plot_scatter(
    df_plot: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: str,
    title: str,
    filename: str,
    max_legend_items: int = 30,
) -> None:
    plt.figure(figsize=(10, 8))

    labels = df_plot[color_col].astype(str).fillna("Missing")
    unique_labels = labels.unique()

    if len(unique_labels) <= max_legend_items:
        for label in unique_labels:
            mask = labels == label
            plt.scatter(
                df_plot.loc[mask, x_col],
                df_plot.loc[mask, y_col],
                s=20,
                alpha=0.75,
                label=label,
            )
        plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    else:
        codes = pd.Categorical(labels).codes
        plt.scatter(
            df_plot[x_col],
            df_plot[y_col],
            c=codes,
            s=20,
            alpha=0.75,
        )

    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(title)
    plt.tight_layout()
    out = OUTPUT_DIR / filename
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"[saved] {out}")


def evaluate_clustering(
    X: pd.DataFrame,
    cluster_labels: np.ndarray,
    metadata: pd.DataFrame,
    label_cols: Iterable[str],
    method_name: str,
) -> pd.DataFrame:
    rows = []

    valid_cluster_mask = cluster_labels != -1
    n_clusters = len(set(cluster_labels[valid_cluster_mask])) if valid_cluster_mask.any() else 0

    sil = np.nan
    if n_clusters >= 2 and valid_cluster_mask.sum() > n_clusters:
        sil = silhouette_score(X.loc[valid_cluster_mask], cluster_labels[valid_cluster_mask])

    for col in label_cols:
        valid = metadata[col].notna()
        ari = adjusted_rand_score(
            metadata.loc[valid, col].astype(str),
            pd.Series(cluster_labels, index=metadata.index).loc[valid].astype(str),
        )
        rows.append({
            "method": method_name,
            "label_reference": col,
            "n_clusters_excluding_noise": n_clusters,
            "silhouette": sil,
            "ARI": ari,
        })

    return pd.DataFrame(rows)


# ============================================================
# Loading and merging
# ============================================================

print_header("1. LOAD DATA")

features = pd.read_csv(FEATURES_PATH)
metadata = pd.read_csv(METADATA_PATH)

print(f"Features shape: {features.shape}")
print(f"Metadata shape: {metadata.shape}")

if SAMPLE_COL not in features.columns:
    raise ValueError(f"{SAMPLE_COL!r} not found in features file")
if SAMPLE_COL not in metadata.columns:
    raise ValueError(f"{SAMPLE_COL!r} not found in metadata file")

features[SAMPLE_COL] = normalize_labels(features[SAMPLE_COL])
metadata[SAMPLE_COL] = normalize_labels(metadata[SAMPLE_COL])

dup_features = features[SAMPLE_COL].duplicated().sum()
dup_metadata = metadata[SAMPLE_COL].duplicated().sum()

print(f"Duplicate Sample IDs in features: {dup_features}")
print(f"Duplicate Sample IDs in metadata: {dup_metadata}")

if dup_features > 0 or dup_metadata > 0:
    raise ValueError("Duplicate Sample IDs detected. Resolve duplicates before proceeding.")

merged = metadata.merge(features, on=SAMPLE_COL, how="inner")
print(f"Merged shape: {merged.shape}")

missing_in_metadata = set(features[SAMPLE_COL]) - set(metadata[SAMPLE_COL])
missing_in_features = set(metadata[SAMPLE_COL]) - set(features[SAMPLE_COL])

print(f"Samples in features but not metadata: {len(missing_in_metadata)}")
print(f"Samples in metadata but not features: {len(missing_in_features)}")

feature_cols = [c for c in features.columns if c != SAMPLE_COL]

# Ensure numeric features
merged[feature_cols] = merged[feature_cols].apply(pd.to_numeric, errors="coerce")

n_missing_feature_values = merged[feature_cols].isna().sum().sum()
print(f"Total missing numeric feature values after coercion: {n_missing_feature_values}")

# Save merged copy
save_df(merged, "merged_harmonized_dataset.csv")


# ============================================================
# Sanity checks and dataset summary
# ============================================================

print_header("2. SANITY CHECKS")

for col in [STUDY_COL, TARGET_COL, SEX_COL]:
    if col not in merged.columns:
        raise ValueError(f"Required metadata column {col!r} not found")

merged[STUDY_COL] = normalize_labels(merged[STUDY_COL])
merged[TARGET_COL] = normalize_labels(merged[TARGET_COL])
merged[SEX_COL] = normalize_labels(merged[SEX_COL])

summary = {
    "n_samples": len(merged),
    "n_features": len(feature_cols),
    "n_studies": merged[STUDY_COL].nunique(dropna=True),
    "studies": sorted(merged[STUDY_COL].dropna().unique().tolist()),
    "condition_counts": merged[TARGET_COL].value_counts(dropna=False).to_dict(),
    "gender_counts": merged[SEX_COL].value_counts(dropna=False).to_dict(),
    "age_missing": int(merged["Age"].isna().sum()) if "Age" in merged.columns else None,
    "bmi_missing": int(merged["BMI"].isna().sum()) if "BMI" in merged.columns else None,
}

for k, v in summary.items():
    print(f"{k}: {v}")

study_table = (
    merged.groupby([STUDY_COL, TARGET_COL, SEX_COL], dropna=False)
    .size()
    .reset_index(name="n")
    .sort_values([STUDY_COL, TARGET_COL, SEX_COL])
)
save_df(study_table, "study_condition_gender_counts.csv")

feature_summary = pd.DataFrame({
    "feature": feature_cols,
    "mean": merged[feature_cols].mean(axis=0).values,
    "std": merged[feature_cols].std(axis=0).values,
    "min": merged[feature_cols].min(axis=0).values,
    "max": merged[feature_cols].max(axis=0).values,
})
save_df(feature_summary, "feature_summary_stats.csv")


# ============================================================
# Binary CRC label creation
# ============================================================

print_header("3. BUILD CRC-vs-CONTROL LABEL")

merged["crc_binary"] = infer_binary_crc_label(merged[TARGET_COL])

label_counts = merged["crc_binary"].value_counts(dropna=False).to_dict()
print("crc_binary counts:", label_counts)

analysis_df = merged.loc[merged["crc_binary"].notna()].copy()
analysis_df["crc_binary"] = analysis_df["crc_binary"].astype(int)

print(f"Samples retained for binary CRC/control analysis: {len(analysis_df)}")
print("Condition values retained:")
print(analysis_df[TARGET_COL].value_counts())

save_df(analysis_df[[SAMPLE_COL, STUDY_COL, TARGET_COL, SEX_COL, "crc_binary"]], "analysis_labels_binary_crc.csv")


# ============================================================
# Matrix for unsupervised analysis
# ============================================================

print_header("4. PREPARE ANALYSIS MATRIX")

X = analysis_df[feature_cols].copy()

if X.isna().sum().sum() > 0:
    # CLR-transformed data ideally should not have missing values, but guard anyway.
    X = X.fillna(X.median(axis=0))

print(f"Final X shape: {X.shape}")
print(f"Feature matrix NaN count: {int(X.isna().sum().sum())}")

# Optional scaling before some methods.
# CLR is already centered/log-ratio transformed, but scaling can still help for t-SNE / KMeans.
X_scaled = pd.DataFrame(
    StandardScaler().fit_transform(X),
    index=X.index,
    columns=X.columns,
)


# ============================================================
# PCA
# ============================================================

print_header("5. PCA")

pca = PCA(n_components=10, random_state=RANDOM_STATE)
X_pca = pca.fit_transform(X_scaled)

explained = pd.DataFrame({
    "PC": [f"PC{i+1}" for i in range(X_pca.shape[1])],
    "explained_variance_ratio": pca.explained_variance_ratio_,
    "cumulative_explained_variance": np.cumsum(pca.explained_variance_ratio_),
})
save_df(explained, "pca_explained_variance.csv")

pca_df = analysis_df[[SAMPLE_COL, STUDY_COL, TARGET_COL, SEX_COL]].copy()
for i in range(X_pca.shape[1]):
    pca_df[f"PC{i+1}"] = X_pca[:, i]

save_df(pca_df, "pca_coordinates.csv")

print(explained.head(10))

plot_scatter(
    pca_df, "PC1", "PC2", STUDY_COL,
    "PCA (PC1 vs PC2) colored by Study",
    "pca_by_study.png"
)
plot_scatter(
    pca_df, "PC1", "PC2", TARGET_COL,
    "PCA (PC1 vs PC2) colored by Condition",
    "pca_by_condition.png"
)
plot_scatter(
    pca_df, "PC1", "PC2", SEX_COL,
    "PCA (PC1 vs PC2) colored by Gender",
    "pca_by_gender.png"
)


# ============================================================
# t-SNE
# ============================================================

print_header("6. T-SNE")

tsne = TSNE(
    n_components=2,
    perplexity=30,
    learning_rate="auto",
    init="pca",
    random_state=RANDOM_STATE,
)
X_tsne = tsne.fit_transform(X_scaled)

tsne_df = analysis_df[[SAMPLE_COL, STUDY_COL, TARGET_COL, SEX_COL]].copy()
tsne_df["TSNE1"] = X_tsne[:, 0]
tsne_df["TSNE2"] = X_tsne[:, 1]

save_df(tsne_df, "tsne_coordinates.csv")

plot_scatter(
    tsne_df, "TSNE1", "TSNE2", STUDY_COL,
    "t-SNE colored by Study",
    "tsne_by_study.png"
)
plot_scatter(
    tsne_df, "TSNE1", "TSNE2", TARGET_COL,
    "t-SNE colored by Condition",
    "tsne_by_condition.png"
)
plot_scatter(
    tsne_df, "TSNE1", "TSNE2", SEX_COL,
    "t-SNE colored by Gender",
    "tsne_by_gender.png"
)


# ============================================================
# UMAP
# ============================================================

print_header("7. UMAP")

if HAS_UMAP:
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=15,
        min_dist=0.1,
        metric="euclidean",
        random_state=RANDOM_STATE,
    )
    X_umap = reducer.fit_transform(X_scaled)

    umap_df = analysis_df[[SAMPLE_COL, STUDY_COL, TARGET_COL, SEX_COL]].copy()
    umap_df["UMAP1"] = X_umap[:, 0]
    umap_df["UMAP2"] = X_umap[:, 1]

    save_df(umap_df, "umap_coordinates.csv")

    plot_scatter(
        umap_df, "UMAP1", "UMAP2", STUDY_COL,
        "UMAP colored by Study",
        "umap_by_study.png"
    )
    plot_scatter(
        umap_df, "UMAP1", "UMAP2", TARGET_COL,
        "UMAP colored by Condition",
        "umap_by_condition.png"
    )
    plot_scatter(
        umap_df, "UMAP1", "UMAP2", SEX_COL,
        "UMAP colored by Gender",
        "umap_by_gender.png"
    )
else:
    print("UMAP not installed. Install with: pip install umap-learn")


# ============================================================
# Clustering and alignment with study/sex/condition
# ============================================================

print_header("8. CLUSTERING")

cluster_results = []

# KMeans with k=2..6
for k in range(2, 7):
    km = KMeans(n_clusters=k, random_state=RANDOM_STATE, n_init=20)
    km_labels = km.fit_predict(X_scaled)

    km_eval = evaluate_clustering(
        X=X_scaled,
        cluster_labels=km_labels,
        metadata=analysis_df,
        label_cols=[STUDY_COL, SEX_COL, TARGET_COL],
        method_name=f"KMeans_k{k}",
    )
    cluster_results.append(km_eval)

# DBSCAN
db = DBSCAN(eps=8.0, min_samples=10)
db_labels = db.fit_predict(X_scaled)

db_eval = evaluate_clustering(
    X=X_scaled,
    cluster_labels=db_labels,
    metadata=analysis_df,
    label_cols=[STUDY_COL, SEX_COL, TARGET_COL],
    method_name="DBSCAN",
)
cluster_results.append(db_eval)

cluster_results_df = pd.concat(cluster_results, ignore_index=True)
save_df(cluster_results_df, "clustering_evaluation.csv")
print(cluster_results_df.sort_values(["label_reference", "ARI"], ascending=[True, False]))


# ============================================================
# Quick study composition audit
# ============================================================

print_header("9. STUDY AUDIT")

study_sizes = analysis_df[STUDY_COL].value_counts().rename_axis(STUDY_COL).reset_index(name="n_samples")
save_df(study_sizes, "study_sizes_binary_crc.csv")
print(study_sizes)

sex_by_study = pd.crosstab(analysis_df[STUDY_COL], analysis_df[SEX_COL], dropna=False)
sex_by_study.to_csv(OUTPUT_DIR / "sex_by_study.csv")
print("[saved]", OUTPUT_DIR / "sex_by_study.csv")

condition_by_study = pd.crosstab(analysis_df[STUDY_COL], analysis_df[TARGET_COL], dropna=False)
condition_by_study.to_csv(OUTPUT_DIR / "condition_by_study.csv")
print("[saved]", OUTPUT_DIR / "condition_by_study.csv")


# ============================================================
# End
# ============================================================

print_header("DONE")

print("Outputs written to:", OUTPUT_DIR.resolve())
print("Next step after reviewing plots:")
print("1) confirm the binary label mapping is correct")
print("2) confirm whether any studies are too imbalanced by sex or condition")
print("3) proceed to differential abundance and LOSO supervised modeling")
