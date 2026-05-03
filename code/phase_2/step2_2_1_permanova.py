"""
ΒΗΜΑ 2.2 — PERMANOVA: Ποσοτικοποίηση πηγών διακύμανσης
=======================================================

Εκτελεί PERMANOVA πάνω στο Hurdle-encoded feature matrix (X_explore)
για να ποσοτικοποιήσει πόση διακύμανση εξηγεί κάθε παράγοντας:
  - Study     → batch effect
  - Condition → CRC signal
  - Gender    → sex variation
  - Age       → age effect (binned)

Αναλύσεις:
  A) Marginal — ένας παράγοντας τη φορά → ανεξάρτητο R² ανά παράγοντα
  B) Partial  — Condition controlling for Study (condition | study)
                → Εξηγεί η Condition κάτι πέρα από το batch effect;

Distance metric: Euclidean πάνω σε StandardScaled Hurdle features.
"""

from __future__ import annotations

import os
import warnings
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler
from skbio import DistanceMatrix
from skbio.stats.distance import permanova

warnings.filterwarnings("ignore")

# ============================================================
# Paths
# ============================================================

ROOT_DIR        = "/home/michalis/Documents/ece_ntua/8th/health"
X_EXPLORE_PATH  = os.path.join(ROOT_DIR, "data/crc_study_final_data/species_level/X_explore.csv")
METADATA_PATH   = os.path.join(ROOT_DIR, "data/crc_study_final_data/species_level/metadata.csv")
OUT_DIR         = os.path.join(ROOT_DIR, "OUTPUTS/phase_2/2.2")
os.makedirs(OUT_DIR, exist_ok=True)

N_PERMUTATIONS = 999

# ============================================================
# 1. Load & Merge
# ============================================================

print("=" * 70)
print("ΒΗΜΑ 2.2 — PERMANOVA")
print("=" * 70)

print("\n[1] Loading data...")
X_df    = pd.read_csv(X_EXPLORE_PATH)
meta_df = pd.read_csv(METADATA_PATH)

df = pd.merge(X_df, meta_df, on="Sample", how="inner")
print(f"    Merged shape: {df.shape}")

feature_cols = [c for c in X_df.columns if c.startswith("PA_") or c.startswith("rCLR_")]
print(f"    Features: {len(feature_cols)}")

# ============================================================
# 2. Handle missing values
# ============================================================

print("\n[2] Handling missing values...")
factors = ["Study", "Condition", "Gender", "Age"]

n_before = len(df)
df_clean = df.dropna(subset=factors).copy()

# Bin Age into decades for categorical PERMANOVA
df_clean["Age_bin"] = pd.cut(
    df_clean["Age"],
    bins=[0, 40, 50, 60, 70, 120],
    labels=["<40", "40-49", "50-59", "60-69", "70+"],
)

df_clean = df_clean.set_index("Sample")
n_after = len(df_clean)
print(f"    Dropped {n_before - n_after} rows with NaN in {factors}")
print(f"    Retained: {n_after} samples")

# ============================================================
# 3. Distance matrix
# ============================================================

print("\n[3] Computing Euclidean distance matrix...")
X = df_clean[feature_cols].values
X_scaled = StandardScaler().fit_transform(X)

dist_condensed = pdist(X_scaled, metric="euclidean")
dist_square    = squareform(dist_condensed)

sample_ids = df_clean.index.astype(str).tolist()
dm = DistanceMatrix(dist_square, ids=sample_ids)
print(f"    Distance matrix: {dm.shape[0]}×{dm.shape[1]}")

# ============================================================
# 4A. Marginal PERMANOVA — one factor at a time
# ============================================================

print("\n[4A] Marginal PERMANOVA (single-factor tests)...")
print(f"     Permutations: {N_PERMUTATIONS}")
print("-" * 65)

test_factors = {
    "Study":     df_clean["Study"].astype(str),
    "Condition": df_clean["Condition"].astype(str),
    "Gender":    df_clean["Gender"].astype(str),
    "Age_bin":   df_clean["Age_bin"].astype(str),
}

marginal_rows = []
for name, grouping in test_factors.items():
    n_groups = grouping.nunique()

    t0 = time.time()
    result = permanova(dm, grouping, permutations=N_PERMUTATIONS)
    elapsed = time.time() - t0

    F_stat  = result["test statistic"]
    p_val   = result["p-value"]
    # R² = SS_between / SS_total = F*(k-1) / (F*(k-1) + N-k)
    r2 = (F_stat * (n_groups - 1)) / (F_stat * (n_groups - 1) + (n_after - n_groups))

    row = {
        "Factor":      name,
        "n_groups":    n_groups,
        "pseudo_F":    round(F_stat, 4),
        "p_value":     p_val,
        "R2":          round(r2, 6),
        "R2_pct":      f"{r2 * 100:.2f}%",
        "permutations": N_PERMUTATIONS,
    }
    marginal_rows.append(row)
    sig = "***" if p_val <= 0.001 else ("**" if p_val <= 0.01 else ("*" if p_val <= 0.05 else "n.s."))
    print(f"  {name:12s}  R²={r2:.4f} ({r2*100:.2f}%)  F={F_stat:8.2f}  p={p_val:.3f} {sig}  [{elapsed:.1f}s]")

marginal_df = pd.DataFrame(marginal_rows)
marginal_path = os.path.join(OUT_DIR, "permanova_marginal.csv")
marginal_df.to_csv(marginal_path, index=False)
print(f"\n  → Saved: {marginal_path}")

# ============================================================
# 4B. Partial test — Condition WITHIN each Study
#     Uses Study as 'strata' by testing within-study distances
# ============================================================

print("\n[4B] Partial PERMANOVA — Condition | Study")
print("     (testing condition effect within each study)")
print("-" * 65)

# Test condition within each study that has both CRC and Control
partial_rows = []
studies = df_clean["Study"].unique()

for study in sorted(studies):
    mask = df_clean["Study"] == study
    sub = df_clean.loc[mask]
    conditions = sub["Condition"].unique()

    # Need at least 2 conditions and enough samples per group
    if len(conditions) < 2 or len(sub) < 10:
        print(f"  {study:30s}  SKIPPED (conditions={list(conditions)}, n={len(sub)})")
        continue

    sub_X = X_scaled[np.where(mask.values)[0]]
    sub_dist = squareform(pdist(sub_X, metric="euclidean"))
    sub_ids  = sub.index.astype(str).tolist()
    sub_dm   = DistanceMatrix(sub_dist, ids=sub_ids)
    sub_grouping = sub["Condition"].astype(str)

    n_groups = sub_grouping.nunique()
    n_sub = len(sub)

    t0 = time.time()
    result = permanova(sub_dm, sub_grouping, permutations=N_PERMUTATIONS)
    elapsed = time.time() - t0

    F_stat = result["test statistic"]
    p_val  = result["p-value"]
    r2 = (F_stat * (n_groups - 1)) / (F_stat * (n_groups - 1) + (n_sub - n_groups))

    sig = "***" if p_val <= 0.001 else ("**" if p_val <= 0.01 else ("*" if p_val <= 0.05 else "n.s."))
    print(f"  {study:30s}  n={n_sub:4d}  R²={r2:.4f} ({r2*100:.2f}%)  F={F_stat:6.2f}  p={p_val:.3f} {sig}")

    partial_rows.append({
        "Study":      study,
        "n_samples":  n_sub,
        "n_conditions": n_groups,
        "conditions": ", ".join(sorted(conditions)),
        "pseudo_F":   round(F_stat, 4),
        "p_value":    p_val,
        "R2":         round(r2, 6),
        "R2_pct":     f"{r2 * 100:.2f}%",
    })

partial_df = pd.DataFrame(partial_rows)
partial_path = os.path.join(OUT_DIR, "permanova_condition_per_study.csv")
partial_df.to_csv(partial_path, index=False)
print(f"\n  → Saved: {partial_path}")

# ============================================================
# 5. Bar plot — Marginal R² per factor
# ============================================================

print("\n[5] Generating R² bar plot...")

fig, ax = plt.subplots(figsize=(8, 5))

colors = {
    "Study":     "#e74c3c",
    "Condition": "#2ecc71",
    "Gender":    "#3498db",
    "Age_bin":   "#f39c12",
}

bars = ax.barh(
    marginal_df["Factor"],
    marginal_df["R2"] * 100,
    color=[colors.get(f, "#95a5a6") for f in marginal_df["Factor"]],
    edgecolor="white",
    linewidth=1.5,
    height=0.55,
)

for bar, row in zip(bars, marginal_df.itertuples()):
    ax.text(
        bar.get_width() + 0.15,
        bar.get_y() + bar.get_height() / 2,
        f"{row.R2 * 100:.2f}%  (p={row.p_value:.3f})",
        va="center",
        fontsize=10,
        fontweight="bold",
    )

ax.set_xlabel("Explained Variance (R² %)", fontsize=12)
ax.set_title("PERMANOVA — Marginal R² per Factor\n(Euclidean distance, Hurdle-encoded features)", fontsize=13, pad=15)
ax.set_xlim(0, max(marginal_df["R2"] * 100) * 1.5)
ax.invert_yaxis()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()

barplot_path = os.path.join(OUT_DIR, "permanova_marginal_R2_barplot.png")
plt.savefig(barplot_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"  → Saved: {barplot_path}")

# ============================================================
# 6. Bar plot — Condition R² within each study
# ============================================================

print("[6] Generating per-study condition R² plot...")

fig, ax = plt.subplots(figsize=(10, 6))

partial_sorted = partial_df.sort_values("R2", ascending=True)
bar_colors = ["#2ecc71" if p <= 0.05 else "#bdc3c7" for p in partial_sorted["p_value"]]

bars = ax.barh(
    partial_sorted["Study"],
    partial_sorted["R2"] * 100,
    color=bar_colors,
    edgecolor="white",
    linewidth=1.5,
    height=0.55,
)

for bar, (_, row) in zip(bars, partial_sorted.iterrows()):
    sig = "***" if row["p_value"] <= 0.001 else ("**" if row["p_value"] <= 0.01 else ("*" if row["p_value"] <= 0.05 else "n.s."))
    ax.text(
        bar.get_width() + 0.1,
        bar.get_y() + bar.get_height() / 2,
        f"{row['R2'] * 100:.2f}%  p={row['p_value']:.3f} {sig}",
        va="center",
        fontsize=9,
    )

ax.set_xlabel("Explained Variance (R² %)", fontsize=12)
ax.set_title("PERMANOVA — Condition Effect Within Each Study\n(green = significant at p≤0.05)", fontsize=13, pad=15)
ax.set_xlim(0, max(partial_sorted["R2"] * 100) * 1.6 + 1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()

per_study_path = os.path.join(OUT_DIR, "permanova_condition_per_study_barplot.png")
plt.savefig(per_study_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"  → Saved: {per_study_path}")

# ============================================================
# 7. Summary
# ============================================================

print("\n" + "=" * 70)
print("SUMMARY — Marginal PERMANOVA")
print("=" * 70)
print(marginal_df[["Factor", "R2_pct", "pseudo_F", "p_value"]].to_string(index=False))

print("\n" + "=" * 70)
print("SUMMARY — Condition | Study (per-study)")
print("=" * 70)
print(partial_df[["Study", "n_samples", "R2_pct", "p_value"]].to_string(index=False))

print(f"\n✓ Step 2.2 complete. Outputs in: {OUT_DIR}")
