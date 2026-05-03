#!/usr/bin/env python3
"""
Step 2.4 — Unsupervised clustering for CRC microbiome exploratory analysis

Run location
------------
Place this file inside your project's code/ directory and run it from there:

    cd /absolute/path/to/project/code
    python step24_unsupervised_clustering.py

Expected project layout
-----------------------
project/
├── code/
│   └── step24_unsupervised_clustering.py
├── data/
│   └── crc_study_final_data/
│       └── species_level/
│           ├── X_explore.csv
│           └── metadata.csv
└── OUTPUTS/
    └── step_2_4/
        ├── step24_clustering_overview.png
        ├── step24_adenoma_k2_assignment.png
        └── step24_final_report.txt

Only important files are written:
1. step24_clustering_overview.png
2. step24_adenoma_k2_assignment.png
3. step24_final_report.txt

Purpose
-------
Run Step 2.4 from the project plan:
- K-Means with k=2, k=3, and k=number_of_studies
- DBSCAN with automatic eps sweep
- Evaluate clustering using silhouette score and ARI vs study, sex, condition
- Specifically answer: with k=2, do Adenoma samples fall with CRC or Control?
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN, KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


# =============================================================================
# ABSOLUTE PATHS, DERIVED FROM THE LOCATION OF THIS SCRIPT
# =============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

DATA_DIR = (PROJECT_ROOT / "data" / "crc_study_final_data" / "species_level").resolve()
OUTPUT_DIR = (PROJECT_ROOT / "OUTPUTS" / "phase_2" / "2.4").resolve()

X_PATH = (DATA_DIR / "X_explore.csv").resolve()
METADATA_PATH = (DATA_DIR / "metadata.csv").resolve()

OVERVIEW_FIGURE_PATH = (OUTPUT_DIR / "step2_4_clustering_overview.png").resolve()
ADENOMA_FIGURE_PATH = (OUTPUT_DIR / "step2_4_adenoma_k2_assignment.png").resolve()
REPORT_PATH = (OUTPUT_DIR / "step2_4_final_report.txt").resolve()


# =============================================================================
# CONFIGURATION
# =============================================================================

SAMPLE_COL = "Sample"
CONDITION_COL = "Condition"
STUDY_COL = "Study"
SEX_COL = "Gender"   # Your metadata uses Gender, not Sex.

RANDOM_STATE = 42
KMEANS_N_INIT = 1  # Reproducible k-means++ run. Increase to 10 for a more exhaustive but slower run.
DBSCAN_MIN_SAMPLES = 10
DBSCAN_EPS_QUANTILES = np.linspace(0.60, 0.98, 9)
PC_VARIANCE_FOR_DBSCAN = 0.80
DBSCAN_MAX_PCS = 20

CONDITION_ORDER = ["Control", "Adenoma", "CRC"]


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def load_and_prepare() -> tuple[pd.DataFrame, list[str], np.ndarray, np.ndarray, pd.DataFrame]:
    """Load X_explore + metadata, merge, scale feature matrix, and compute PCA scores."""
    if not X_PATH.exists():
        raise FileNotFoundError(f"Cannot find X_explore.csv at: {X_PATH}")
    if not METADATA_PATH.exists():
        raise FileNotFoundError(f"Cannot find metadata.csv at: {METADATA_PATH}")

    x = pd.read_csv(X_PATH)
    meta = pd.read_csv(METADATA_PATH)

    required_x = {SAMPLE_COL}
    required_meta = {SAMPLE_COL, CONDITION_COL, STUDY_COL, SEX_COL}

    missing_x = required_x - set(x.columns)
    missing_meta = required_meta - set(meta.columns)
    if missing_x:
        raise ValueError(f"Missing columns in X_explore.csv: {sorted(missing_x)}")
    if missing_meta:
        raise ValueError(f"Missing columns in metadata.csv: {sorted(missing_meta)}")

    if x[SAMPLE_COL].duplicated().any():
        examples = x.loc[x[SAMPLE_COL].duplicated(), SAMPLE_COL].head().tolist()
        raise ValueError(f"Duplicate Sample IDs in X_explore.csv. Examples: {examples}")
    if meta[SAMPLE_COL].duplicated().any():
        examples = meta.loc[meta[SAMPLE_COL].duplicated(), SAMPLE_COL].head().tolist()
        raise ValueError(f"Duplicate Sample IDs in metadata.csv. Examples: {examples}")

    feature_cols = [c for c in x.columns if c != SAMPLE_COL]
    non_numeric = [c for c in feature_cols if not pd.api.types.is_numeric_dtype(x[c])]
    if non_numeric:
        warnings.warn(f"Dropping non-numeric feature columns: {non_numeric[:10]}")
        feature_cols = [c for c in feature_cols if c not in non_numeric]

    merged = x.merge(meta, on=SAMPLE_COL, how="inner")
    if merged.empty:
        raise ValueError("Merge produced zero rows. Check Sample IDs.")

    dropped_x = x.shape[0] - merged.shape[0]
    dropped_meta = meta.shape[0] - merged.shape[0]
    if dropped_x or dropped_meta:
        warnings.warn(
            f"Merge did not keep all rows: dropped {dropped_x} X rows and "
            f"{dropped_meta} metadata rows."
        )

    before = merged.shape[0]
    merged = merged.dropna(subset=[CONDITION_COL, STUDY_COL, SEX_COL]).copy()
    after = merged.shape[0]
    if before != after:
        warnings.warn(f"Dropped {before - after} rows with missing condition/study/sex.")

    missing_features = int(merged[feature_cols].isna().sum().sum())
    if missing_features:
        warnings.warn(f"Found {missing_features} missing feature values. Filling with 0.")
        merged[feature_cols] = merged[feature_cols].fillna(0)

    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(merged[feature_cols].to_numpy())

    pca = PCA(random_state=RANDOM_STATE)
    pca_scores_array = pca.fit_transform(x_scaled)
    pc_cols = [f"PC{i + 1}" for i in range(pca_scores_array.shape[1])]
    pca_scores = pd.DataFrame(pca_scores_array, columns=pc_cols)
    pca_scores[SAMPLE_COL] = merged[SAMPLE_COL].values
    pca_scores[CONDITION_COL] = merged[CONDITION_COL].values
    pca_scores[STUDY_COL] = merged[STUDY_COL].values
    pca_scores[SEX_COL] = merged[SEX_COL].values

    explained = pd.DataFrame({
        "PC": pc_cols,
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "cumulative_explained_variance": np.cumsum(pca.explained_variance_ratio_),
    })

    return merged, feature_cols, x_scaled, pca_scores_array, pca_scores, explained


def ari_against(labels: np.ndarray, truth: pd.Series) -> float:
    """Adjusted Rand Index, ignoring missing truth labels."""
    valid = truth.notna().to_numpy()
    if valid.sum() < 2:
        return np.nan
    if len(np.unique(truth[valid])) < 2:
        return np.nan
    if len(np.unique(labels[valid])) < 2:
        return np.nan
    return float(adjusted_rand_score(truth[valid].astype(str), labels[valid]))


def safe_silhouette(x_matrix: np.ndarray, labels: np.ndarray, ignore_noise: bool = False) -> float:
    """Silhouette score with safeguards. For DBSCAN, optionally ignore noise label -1."""
    labels_for_score = labels.copy()
    x_for_score = x_matrix

    if ignore_noise:
        keep = labels_for_score != -1
        labels_for_score = labels_for_score[keep]
        x_for_score = x_matrix[keep]

    n_labels = len(np.unique(labels_for_score))
    n_samples = len(labels_for_score)

    if n_samples < 3 or n_labels < 2 or n_labels >= n_samples:
        return np.nan

    return float(silhouette_score(
        x_for_score,
        labels_for_score,
        sample_size=min(500, n_samples),
        random_state=RANDOM_STATE,
    ))


def fit_kmeans_models(df: pd.DataFrame, x_scaled: np.ndarray) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    """Run K-Means for k=2, k=3, and k=number of studies."""
    n_studies = int(df[STUDY_COL].nunique())
    k_values = {
        "KMeans_k2": 2,
        "KMeans_k3": 3,
        f"KMeans_k{n_studies}_n_studies": n_studies,
    }

    label_dict = {}
    metric_rows = []

    for model_name, k in k_values.items():
        model = KMeans(n_clusters=k, random_state=RANDOM_STATE, n_init=KMEANS_N_INIT, max_iter=100, algorithm='lloyd')
        labels = model.fit_predict(x_scaled)
        label_dict[model_name] = labels

        metric_rows.append({
            "model": model_name,
            "n_clusters": k,
            "noise_fraction": 0.0,
            "silhouette": safe_silhouette(x_scaled, labels),
            "ARI_vs_study": ari_against(labels, df[STUDY_COL]),
            "ARI_vs_condition": ari_against(labels, df[CONDITION_COL]),
            "ARI_vs_sex": ari_against(labels, df[SEX_COL]),
        })

    return label_dict, pd.DataFrame(metric_rows)


def select_pc80_matrix(pca_scores_array: np.ndarray, explained: pd.DataFrame) -> tuple[np.ndarray, int, float]:
    """Select PCs explaining at least PC_VARIANCE_FOR_DBSCAN of variance."""
    cumulative = explained["cumulative_explained_variance"].to_numpy()
    n_needed = int(np.searchsorted(cumulative, PC_VARIANCE_FOR_DBSCAN, side="left") + 1)
    n_pc80 = min(n_needed, DBSCAN_MAX_PCS)
    achieved = float(cumulative[n_pc80 - 1])
    return pca_scores_array[:, :n_pc80], n_pc80, achieved


def fit_dbscan_with_sweep(
    df: pd.DataFrame,
    x_pc80: np.ndarray,
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    """
    Run DBSCAN over an eps grid derived from k-nearest-neighbor distances.

    Selection rule:
    - prefer solutions with >=2 non-noise clusters
    - avoid trivial solutions with >80% noise
    - maximize silhouette among non-noise samples
    - if no valid solution exists, use the candidate with the largest number of non-noise clusters,
      then lowest noise fraction
    """
    nn = NearestNeighbors(n_neighbors=DBSCAN_MIN_SAMPLES, algorithm='ball_tree')
    nn.fit(x_pc80)
    distances, _ = nn.kneighbors(x_pc80)
    kth_distances = distances[:, -1]
    eps_values = np.quantile(kth_distances, DBSCAN_EPS_QUANTILES)
    eps_values = np.unique(np.round(eps_values, 6))

    rows = []
    label_candidates = {}

    for eps in eps_values:
        model = DBSCAN(eps=float(eps), min_samples=DBSCAN_MIN_SAMPLES, algorithm='ball_tree')
        labels = model.fit_predict(x_pc80)
        label_candidates[float(eps)] = labels

        non_noise = labels != -1
        non_noise_labels = labels[non_noise]
        n_non_noise_clusters = len(set(non_noise_labels))
        noise_fraction = float(np.mean(labels == -1))

        rows.append({
            "eps": float(eps),
            "min_samples": DBSCAN_MIN_SAMPLES,
            "n_clusters_excluding_noise": n_non_noise_clusters,
            "noise_fraction": noise_fraction,
            "silhouette_non_noise": safe_silhouette(x_pc80, labels, ignore_noise=True),
            "ARI_vs_study": ari_against(labels, df[STUDY_COL]),
            "ARI_vs_condition": ari_against(labels, df[CONDITION_COL]),
            "ARI_vs_sex": ari_against(labels, df[SEX_COL]),
        })

    sweep = pd.DataFrame(rows)
    valid = sweep[
        (sweep["n_clusters_excluding_noise"] >= 2)
        & (sweep["noise_fraction"] <= 0.80)
        & (sweep["silhouette_non_noise"].notna())
    ].copy()

    if not valid.empty:
        selected = valid.sort_values(
            ["silhouette_non_noise", "noise_fraction"],
            ascending=[False, True],
        ).iloc[0]
        selection_note = "Selected eps by maximizing non-noise silhouette among non-trivial DBSCAN solutions."
    else:
        selected = sweep.sort_values(
            ["n_clusters_excluding_noise", "noise_fraction"],
            ascending=[False, True],
        ).iloc[0]
        selection_note = "No strong non-trivial DBSCAN solution found; selected least-degenerate candidate."

    selected_eps = float(selected["eps"])
    selected_labels = label_candidates[selected_eps]

    selected_metrics = pd.DataFrame([{
        "model": f"DBSCAN_eps_{selected_eps:.4f}_min{DBSCAN_MIN_SAMPLES}",
        "n_clusters": int(selected["n_clusters_excluding_noise"]),
        "noise_fraction": float(selected["noise_fraction"]),
        "silhouette": float(selected["silhouette_non_noise"]) if pd.notna(selected["silhouette_non_noise"]) else np.nan,
        "ARI_vs_study": ari_against(selected_labels, df[STUDY_COL]),
        "ARI_vs_condition": ari_against(selected_labels, df[CONDITION_COL]),
        "ARI_vs_sex": ari_against(selected_labels, df[SEX_COL]),
        "selection_note": selection_note,
    }])

    return selected_labels, selected_metrics, sweep


def label_k2_clusters_by_crc_control(df: pd.DataFrame, k2_labels: np.ndarray) -> pd.DataFrame:
    """
    Label each k=2 cluster as CRC-like or Control-like using only CRC and Control samples.
    Adenoma is not used to name the clusters; it is evaluated afterwards.
    """
    temp = df[[SAMPLE_COL, CONDITION_COL]].copy()
    temp["cluster"] = k2_labels

    rows = []
    for cluster in sorted(temp["cluster"].unique()):
        subset = temp[temp["cluster"] == cluster]
        crc_n = int(((subset[CONDITION_COL] == "CRC")).sum())
        control_n = int(((subset[CONDITION_COL] == "Control")).sum())
        adenoma_n = int(((subset[CONDITION_COL] == "Adenoma")).sum())
        total_n = int(subset.shape[0])
        denom = crc_n + control_n
        crc_fraction_among_crc_control = crc_n / denom if denom else np.nan
        control_fraction_among_crc_control = control_n / denom if denom else np.nan

        margin = abs(crc_n - control_n) / denom if denom else np.nan
        if denom and margin < 0.05:
            if crc_n > control_n:
                biological_label = "Weakly CRC-like"
            elif control_n > crc_n:
                biological_label = "Weakly Control-like"
            else:
                biological_label = "Mixed/tie"
        elif crc_n > control_n:
            biological_label = "CRC-like"
        elif control_n > crc_n:
            biological_label = "Control-like"
        else:
            biological_label = "Mixed/tie"

        rows.append({
            "cluster": int(cluster),
            "biological_label_without_adenoma": biological_label,
            "total_n": total_n,
            "CRC_n": crc_n,
            "Control_n": control_n,
            "Adenoma_n": adenoma_n,
            "CRC_fraction_among_CRC_Control": crc_fraction_among_crc_control,
            "Control_fraction_among_CRC_Control": control_fraction_among_crc_control,
            "CRC_Control_margin": margin,
        })

    return pd.DataFrame(rows)


def adenoma_assignment_summary(k2_cluster_summary: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    """Summarize where Adenoma falls in k=2."""
    adenoma_total = int(k2_cluster_summary["Adenoma_n"].sum())
    out = k2_cluster_summary.copy()
    out["Adenoma_fraction"] = out["Adenoma_n"] / adenoma_total if adenoma_total else np.nan

    by_label = out.groupby("biological_label_without_adenoma", as_index=False)["Adenoma_n"].sum()
    by_label["Adenoma_fraction"] = by_label["Adenoma_n"] / adenoma_total if adenoma_total else np.nan

    if by_label.empty or adenoma_total == 0:
        conclusion = "No Adenoma samples available."
    else:
        top = by_label.sort_values("Adenoma_n", ascending=False).iloc[0]
        conclusion = (
            f"Most Adenoma samples fall in the {top['biological_label_without_adenoma']} "
            f"k=2 cluster: {int(top['Adenoma_n'])}/{adenoma_total} "
            f"({100 * float(top['Adenoma_fraction']):.1f}%)."
        )

    return by_label, conclusion


def cluster_composition_table(df: pd.DataFrame, labels: np.ndarray, model_name: str) -> pd.DataFrame:
    """Condition composition per cluster."""
    temp = df[[CONDITION_COL]].copy()
    temp["cluster"] = labels
    counts = pd.crosstab(temp["cluster"], temp[CONDITION_COL])
    for condition in CONDITION_ORDER:
        if condition not in counts.columns:
            counts[condition] = 0
    counts = counts[CONDITION_ORDER]
    counts["Total"] = counts.sum(axis=1)
    props = counts[CONDITION_ORDER].div(counts["Total"], axis=0).add_suffix("_fraction")
    out = pd.concat([counts, props], axis=1).reset_index()
    out.insert(0, "model", model_name)
    return out


def plot_overview(
    pca_scores: pd.DataFrame,
    labels_dict: dict[str, np.ndarray],
    dbscan_labels: np.ndarray,
    explained: pd.DataFrame,
) -> None:
    """Create one compact overview figure with KMeans and DBSCAN clusters in PC1-PC2 space."""
    plot_specs = [
        ("KMeans k=2", labels_dict["KMeans_k2"]),
        ("KMeans k=3", labels_dict["KMeans_k3"]),
        ("DBSCAN", dbscan_labels),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)
    x = pca_scores["PC1"].to_numpy()
    y = pca_scores["PC2"].to_numpy()

    for ax, (title, labels) in zip(axes, plot_specs):
        scatter = ax.scatter(x, y, c=labels, s=14, alpha=0.75, cmap="tab20")
        ax.set_title(title)
        ax.set_xlabel(f"PC1 ({100 * explained.loc[0, 'explained_variance_ratio']:.1f}% var.)")
        ax.set_ylabel(f"PC2 ({100 * explained.loc[1, 'explained_variance_ratio']:.1f}% var.)")

        # Overlay Adenoma samples as larger open circles to answer the new question visually.
        adenoma_mask = pca_scores[CONDITION_COL].astype(str).eq("Adenoma").to_numpy()
        ax.scatter(
            x[adenoma_mask],
            y[adenoma_mask],
            facecolors="none",
            edgecolors="black",
            s=40,
            linewidths=0.8,
            label="Adenoma",
        )
        ax.legend(loc="best", frameon=True)

    fig.suptitle("Step 2.4 clustering overview in PC1-PC2 space; Adenoma highlighted with open circles")
    fig.savefig(OVERVIEW_FIGURE_PATH, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_adenoma_k2(adenoma_by_label: pd.DataFrame) -> None:
    """Bar plot showing whether Adenoma falls with CRC-like or Control-like k=2 cluster."""
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    if adenoma_by_label.empty:
        ax.text(0.5, 0.5, "No Adenoma samples", ha="center", va="center")
        ax.set_axis_off()
    else:
        x = adenoma_by_label["biological_label_without_adenoma"].astype(str).tolist()
        y = adenoma_by_label["Adenoma_n"].to_numpy()
        bars = ax.bar(x, y)
        total = int(y.sum())
        for bar, count in zip(bars, y):
            pct = 100 * count / total if total else 0
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                f"{int(count)}\n({pct:.1f}%)",
                ha="center",
                va="bottom",
            )
        ax.set_ylabel("Number of Adenoma samples")
        ax.set_title("Where do Adenoma samples fall when K-Means k=2?")
        ax.set_ylim(0, max(y) * 1.20 if len(y) and max(y) > 0 else 1)

    fig.savefig(ADENOMA_FIGURE_PATH, dpi=300, bbox_inches="tight")
    plt.close(fig)


def format_table(df: pd.DataFrame, float_digits: int = 4) -> str:
    """Format a DataFrame for a plain-text report."""
    if df.empty:
        return "[empty table]"
    out = df.copy()
    for col in out.select_dtypes(include=[float]).columns:
        out[col] = out[col].map(lambda v: "NA" if pd.isna(v) else f"{v:.{float_digits}f}")
    return out.to_string(index=False)


def write_report(
    df: pd.DataFrame,
    feature_cols: list[str],
    explained: pd.DataFrame,
    metrics: pd.DataFrame,
    dbscan_sweep: pd.DataFrame,
    k2_summary: pd.DataFrame,
    adenoma_by_label: pd.DataFrame,
    adenoma_conclusion: str,
    composition_tables: list[pd.DataFrame],
    n_pc80: int,
    pc80_variance: float,
) -> None:
    """Write the final report containing all important numbers and interpretation."""
    n_studies = df[STUDY_COL].nunique()
    condition_counts = df[CONDITION_COL].value_counts().reindex(CONDITION_ORDER).fillna(0).astype(int)
    sex_counts = df[SEX_COL].value_counts(dropna=False)
    study_counts = df[STUDY_COL].value_counts().sort_index()

    best_study_model = metrics.sort_values("ARI_vs_study", ascending=False).iloc[0]
    best_condition_model = metrics.sort_values("ARI_vs_condition", ascending=False).iloc[0]
    best_sex_model = metrics.sort_values("ARI_vs_sex", ascending=False).iloc[0]

    dbscan_selected = metrics[metrics["model"].str.startswith("DBSCAN")].iloc[0]
    composition = pd.concat(composition_tables, ignore_index=True)

    lines = []
    lines.append("STEP 2.4 — UNSUPERVISED CLUSTERING REPORT")
    lines.append("=" * 72)
    lines.append("")
    lines.append("INPUTS")
    lines.append("------")
    lines.append(f"X_explore path: {X_PATH}")
    lines.append(f"metadata path:  {METADATA_PATH}")
    lines.append(f"Output directory: {OUTPUT_DIR}")
    lines.append("")
    lines.append("DATA SANITY CHECK")
    lines.append("-----------------")
    lines.append(f"Merged samples used: {df.shape[0]}")
    lines.append(f"Feature columns used: {len(feature_cols)}")
    lines.append(f"Studies: {n_studies}")
    lines.append("Condition counts:")
    lines.append(condition_counts.to_string())
    lines.append("Sex counts:")
    lines.append(sex_counts.to_string())
    lines.append("Samples per study:")
    lines.append(study_counts.to_string())
    lines.append("")
    lines.append("PCA USED FOR VISUALIZATION AND DBSCAN")
    lines.append("-------------------------------------")
    lines.append(f"PC1 variance: {100 * explained.loc[0, 'explained_variance_ratio']:.3f}%")
    lines.append(f"PC2 variance: {100 * explained.loc[1, 'explained_variance_ratio']:.3f}%")
    lines.append(f"PCs used for DBSCAN: {n_pc80}")
    lines.append(f"Variance explained by DBSCAN PC space: {100 * pc80_variance:.3f}%")
    lines.append("Note: DBSCAN is capped at the first 20 PCs for speed and density-estimation stability; this is an EDA visualization/structure check, not a supervised pipeline step.")
    lines.append("")
    lines.append("CLUSTERING METRICS")
    lines.append("------------------")
    lines.append("Silhouette measures compact/separated clusters. ARI measures agreement with known labels.")
    lines.append("ARI close to 0 means little/no agreement; ARI close to 1 means strong agreement.")
    lines.append(format_table(metrics.drop(columns=["selection_note"], errors="ignore")))
    lines.append("")
    lines.append("MAIN EDA INTERPRETATION")
    lines.append("-----------------------")
    lines.append(
        f"Highest ARI vs study: {best_study_model['model']} "
        f"(ARI={best_study_model['ARI_vs_study']:.4f})."
    )
    lines.append(
        f"Highest ARI vs condition: {best_condition_model['model']} "
        f"(ARI={best_condition_model['ARI_vs_condition']:.4f})."
    )
    lines.append(
        f"Highest ARI vs sex: {best_sex_model['model']} "
        f"(ARI={best_sex_model['ARI_vs_sex']:.4f})."
    )

    if best_study_model["ARI_vs_study"] > best_condition_model["ARI_vs_condition"]:
        lines.append(
            "The unsupervised clustering aligns more strongly with study/batch than with disease condition. "
            "This is expected in multi-cohort microbiome data and supports keeping LODO evaluation later."
        )
    else:
        lines.append(
            "The unsupervised clustering aligns at least as strongly with condition as with study. "
            "This would suggest a relatively visible disease signal in exploratory space."
        )
    lines.append("")
    lines.append("K-MEANS k=2 — WHERE DO ADENOMA SAMPLES FALL?")
    lines.append("------------------------------------------------")
    lines.append("Clusters were named CRC-like or Control-like using only CRC and Control samples.")
    lines.append("Adenoma samples were then assigned afterwards, so Adenoma did not define the cluster names.")
    lines.append(format_table(k2_summary))
    lines.append("")
    lines.append("Adenoma assignment by biological cluster label:")
    lines.append(format_table(adenoma_by_label))
    lines.append(adenoma_conclusion)
    lines.append("")
    lines.append("CONDITION COMPOSITION BY CLUSTER")
    lines.append("--------------------------------")
    lines.append(format_table(composition))
    lines.append("")
    lines.append("DBSCAN PARAMETER SELECTION")
    lines.append("--------------------------")
    lines.append(f"min_samples: {DBSCAN_MIN_SAMPLES}")
    lines.append(f"Selected model: {dbscan_selected['model']}")
    lines.append(f"Selected DBSCAN noise fraction: {dbscan_selected['noise_fraction']:.4f}")
    lines.append(f"Selected DBSCAN non-noise clusters: {int(dbscan_selected['n_clusters'])}")
    if "selection_note" in dbscan_selected:
        lines.append(str(dbscan_selected["selection_note"]))
    lines.append("")
    lines.append("Top 10 DBSCAN eps candidates by non-noise silhouette:")
    top_dbscan = dbscan_sweep.sort_values("silhouette_non_noise", ascending=False).head(10)
    lines.append(format_table(top_dbscan))
    lines.append("")
    lines.append("OUTPUT FILES")
    lines.append("------------")
    lines.append(f"Clustering overview figure: {OVERVIEW_FIGURE_PATH}")
    lines.append(f"Adenoma k=2 assignment figure: {ADENOMA_FIGURE_PATH}")
    lines.append(f"Final report: {REPORT_PATH}")
    lines.append("")
    lines.append("REPORT-READY SUMMARY")
    lines.append("--------------------")
    lines.append(
        "We performed unsupervised clustering on the z-score scaled hurdle-encoded feature matrix. "
        "K-Means was evaluated with k=2, k=3, and k equal to the number of studies. "
        "DBSCAN was run in a capped PCA space, with eps selected "
        "by a nearest-neighbor distance sweep. Cluster quality was summarized using silhouette score, "
        "and agreement with study, condition, and sex labels was measured using adjusted Rand index."
    )
    lines.append(
        f"In k=2 K-Means, {adenoma_conclusion} "
        "However, the k=2 clusters have very similar CRC/control composition, so the CRC-like vs Control-like naming should be treated as weak. "
        "Because this is unsupervised EDA, these clusters should not be interpreted as predictive model "
        "performance; they mainly indicate whether the strongest structure in the data corresponds to "
        "batch/study, sex, or disease condition."
    )

    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df, feature_cols, x_scaled, pca_scores_array, pca_scores, explained = load_and_prepare()

    kmeans_labels, kmeans_metrics = fit_kmeans_models(df, x_scaled)

    x_pc80, n_pc80, pc80_variance = select_pc80_matrix(pca_scores_array, explained)
    dbscan_labels, dbscan_metrics, dbscan_sweep = fit_dbscan_with_sweep(df, x_pc80)

    all_metrics = pd.concat([kmeans_metrics, dbscan_metrics], ignore_index=True)

    k2_summary = label_k2_clusters_by_crc_control(df, kmeans_labels["KMeans_k2"])
    adenoma_by_label, adenoma_conclusion = adenoma_assignment_summary(k2_summary)

    composition_tables = []
    for name, labels in kmeans_labels.items():
        composition_tables.append(cluster_composition_table(df, labels, name))
    composition_tables.append(cluster_composition_table(df, dbscan_labels, "DBSCAN_selected"))

    plot_overview(pca_scores, kmeans_labels, dbscan_labels, explained)
    plot_adenoma_k2(adenoma_by_label)

    write_report(
        df=df,
        feature_cols=feature_cols,
        explained=explained,
        metrics=all_metrics,
        dbscan_sweep=dbscan_sweep,
        k2_summary=k2_summary,
        adenoma_by_label=adenoma_by_label,
        adenoma_conclusion=adenoma_conclusion,
        composition_tables=composition_tables,
        n_pc80=n_pc80,
        pc80_variance=pc80_variance,
    )

    print("Step 2.4 completed successfully.")
    print(f"Wrote: {OVERVIEW_FIGURE_PATH}")
    print(f"Wrote: {ADENOMA_FIGURE_PATH}")
    print(f"Wrote: {REPORT_PATH}")


if __name__ == "__main__":
    main()
