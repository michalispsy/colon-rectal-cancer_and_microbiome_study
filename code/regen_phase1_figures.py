"""
Regenerate Phase 1 figures from the current (post-Gupta) dataset.
Overwrites existing figures in OUTPUTS/phase_1/graph_and_md_files/figures/
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = "/home/michalis/Documents/ece_ntua/8th/health"
META_PATH = os.path.join(ROOT, "data/crc_study_final_data/species_level/metadata.csv")
SPECIES_PATH = os.path.join(ROOT, "data/crc_study_final_data/species_level/species_filtered_crossstudy.csv")
FIG_DIR = os.path.join(ROOT, "OUTPUTS/phase_1/graph_and_md_files/figures")
os.makedirs(FIG_DIR, exist_ok=True)

meta = pd.read_csv(META_PATH)
species = pd.read_csv(SPECIES_PATH)

study_order = meta["Study"].value_counts().index.tolist()
n_total = len(meta)

# ── Colour palettes ──
study_colors = plt.cm.tab20(np.linspace(0, 1, len(study_order)))
cond_colors = {"CRC": "#e74c3c", "Control": "#2ecc71", "Adenoma": "#f39c12"}
sex_colors  = {"Male": "#3498db", "Female": "#e91e63"}

# ════════════════════════════════════════════════════════════════
# Fig 1 — Samples per Study
# ════════════════════════════════════════════════════════════════
print("fig1 — samples per study")
fig, ax = plt.subplots(figsize=(10, 5))
counts = meta["Study"].value_counts().reindex(study_order)
bars = ax.barh(counts.index, counts.values, color=study_colors, edgecolor="white", linewidth=1.2)
for bar, v in zip(bars, counts.values):
    ax.text(v + 5, bar.get_y() + bar.get_height() / 2, str(v), va="center", fontsize=9, fontweight="bold")
ax.set_xlabel("Αριθμός Δειγμάτων", fontsize=11)
ax.set_title(f"Δείγματα ανά Μελέτη  (n = {n_total})", fontsize=13, pad=12)
ax.invert_yaxis()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig1_samples_per_study.png"), dpi=300)
plt.close()

# ════════════════════════════════════════════════════════════════
# Fig 2 — Condition distribution (pie + stacked bar per study)
# ════════════════════════════════════════════════════════════════
print("fig2 — condition distribution")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Pie
cond_counts = meta["Condition"].value_counts()
axes[0].pie(cond_counts.values,
            labels=[f"{c}\n({v}, {v/n_total*100:.1f}%)" for c, v in cond_counts.items()],
            colors=[cond_colors[c] for c in cond_counts.index],
            startangle=90, wedgeprops={"edgecolor": "white", "linewidth": 1.5})
axes[0].set_title("Κατανομή Condition", fontsize=12)

# Stacked bar per study
ct = pd.crosstab(meta["Study"], meta["Condition"]).reindex(study_order)
ct[["CRC", "Control", "Adenoma"]].plot.barh(
    stacked=True, ax=axes[1],
    color=[cond_colors["CRC"], cond_colors["Control"], cond_colors["Adenoma"]],
    edgecolor="white", linewidth=0.8)
axes[1].invert_yaxis()
axes[1].set_xlabel("Δείγματα")
axes[1].set_title("Condition ανά Study", fontsize=12)
axes[1].legend(title="Condition", loc="lower right")
axes[1].spines["top"].set_visible(False)
axes[1].spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig2_condition_distribution.png"), dpi=300)
plt.close()

# ════════════════════════════════════════════════════════════════
# Fig 3 — Gender distribution (pie + stacked bar per study)
# ════════════════════════════════════════════════════════════════
print("fig3 — gender distribution")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sex_counts = meta["Gender"].value_counts()
axes[0].pie(sex_counts.values,
            labels=[f"{s}\n({v}, {v/n_total*100:.1f}%)" for s, v in sex_counts.items()],
            colors=[sex_colors[s] for s in sex_counts.index],
            startangle=90, wedgeprops={"edgecolor": "white", "linewidth": 1.5})
axes[0].set_title("Κατανομή Gender", fontsize=12)

ct_sex = pd.crosstab(meta["Study"], meta["Gender"]).reindex(study_order)
ct_sex[["Male", "Female"]].plot.barh(
    stacked=True, ax=axes[1],
    color=[sex_colors["Male"], sex_colors["Female"]],
    edgecolor="white", linewidth=0.8)
axes[1].invert_yaxis()
axes[1].set_xlabel("Δείγματα")
axes[1].set_title("Gender ανά Study", fontsize=12)
axes[1].legend(title="Gender", loc="lower right")
axes[1].spines["top"].set_visible(False)
axes[1].spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig3_gender_distribution.png"), dpi=300)
plt.close()

# ════════════════════════════════════════════════════════════════
# Fig 4 — Age distribution per study (boxplot)
# ════════════════════════════════════════════════════════════════
print("fig4 — age distribution")
fig, ax = plt.subplots(figsize=(10, 5))
data_age = [meta.loc[meta["Study"] == s, "Age"].dropna().values for s in study_order]
bp = ax.boxplot(data_age, vert=False, labels=study_order, patch_artist=True, widths=0.6)
for patch, color in zip(bp["boxes"], study_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.75)
ax.set_xlabel("Ηλικία (έτη)", fontsize=11)
ax.set_title("Κατανομή Ηλικίας ανά Μελέτη", fontsize=13, pad=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig4_age_distribution.png"), dpi=300)
plt.close()

# ════════════════════════════════════════════════════════════════
# Fig 5 — Condition × Gender (grouped bar)
# ════════════════════════════════════════════════════════════════
print("fig5 — condition × gender")
fig, ax = plt.subplots(figsize=(7, 4))
ct_cg = pd.crosstab(meta["Condition"], meta["Gender"])
ct_cg[["Male", "Female"]].plot.bar(
    ax=ax, color=[sex_colors["Male"], sex_colors["Female"]],
    edgecolor="white", linewidth=1.2)
ax.set_ylabel("Δείγματα")
ax.set_title("Condition × Gender", fontsize=12)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
ax.legend(title="Gender")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig5_condition_gender.png"), dpi=300)
plt.close()

# ════════════════════════════════════════════════════════════════
# Fig 6 — Filtering funnel
# ════════════════════════════════════════════════════════════════
print("fig6 — filtering funnel")
fig, ax = plt.subplots(figsize=(8, 5))
steps = [
    ("Raw species (curatedMetagenomicData)", 934),
    ("Prevalence filter (≥10% in ≥6 studies)", 107),
    ("Post Gupta removal (10 studies)", 107),
]
y_positions = [3, 2, 1]
max_w = steps[0][1]
bar_height = 0.5

for (label, val), y in zip(steps, y_positions):
    width = val / max_w * 0.9
    left = (1 - width) / 2
    color = "#2ecc71" if y == 1 else ("#3498db" if y == 2 else "#95a5a6")
    ax.barh(y, width, left=left, height=bar_height, color=color, edgecolor="white", linewidth=2)
    ax.text(0.5, y, f"{label}\n({val})", ha="center", va="center", fontsize=10, fontweight="bold")

ax.set_xlim(0, 1)
ax.set_ylim(0.3, 3.7)
ax.axis("off")
ax.set_title("Species Filtering Funnel", fontsize=13, pad=15)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig6_filtering_funnel.png"), dpi=300)
plt.close()

# ════════════════════════════════════════════════════════════════
# Fig 7 — BMI distribution (histogram)
# ════════════════════════════════════════════════════════════════
print("fig7 — BMI distribution")
fig, ax = plt.subplots(figsize=(7, 4))
bmi = meta["BMI"].dropna()
ax.hist(bmi, bins=30, color="#9b59b6", edgecolor="white", linewidth=0.8, alpha=0.85)
ax.axvline(bmi.median(), color="#e74c3c", ls="--", lw=1.5, label=f"Median = {bmi.median():.1f}")
ax.set_xlabel("BMI")
ax.set_ylabel("Πλήθος")
ax.set_title(f"Κατανομή BMI  (n = {len(bmi)}, missing = {meta['BMI'].isna().sum()})", fontsize=12)
ax.legend()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig7_bmi_distribution.png"), dpi=300)
plt.close()

# ════════════════════════════════════════════════════════════════
# Fig 8 — Sparsity per species
# ════════════════════════════════════════════════════════════════
print("fig8 — sparsity")
feature_cols = [c for c in species.columns if c != "Sample"]
zero_frac = (species[feature_cols] == 0).mean(axis=0).sort_values(ascending=False)

total_cells = species[feature_cols].size
total_zeros = (species[feature_cols] == 0).sum().sum()
sparsity = total_zeros / total_cells * 100

fig, ax = plt.subplots(figsize=(10, 4))
ax.bar(range(len(zero_frac)), zero_frac.values * 100, color="#e67e22", width=1.0, edgecolor="none")
ax.axhline(50, color="#e74c3c", ls="--", lw=1, alpha=0.7, label="50% threshold")
ax.set_xlabel("Species (ταξινομημένα κατά sparsity)", fontsize=11)
ax.set_ylabel("% Μηδενικών", fontsize=11)
ax.set_title(f"Sparsity ανά Species  (overall = {sparsity:.1f}%,  cells = {total_cells:,},  zeros = {total_zeros:,})", fontsize=12)
ax.legend()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig8_sparsity.png"), dpi=300)
plt.close()

print(f"\n✓ All 8 figures regenerated in {FIG_DIR}")
