import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

INPUT_DIR = "OUTPUTS/phase_4/4.1/4.1_results"
AGG_DIR = "OUTPUTS/phase_4/4.2_aggregation"
FIG_DIR = os.path.join(AGG_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

df = pd.read_csv(os.path.join(INPUT_DIR, "advanced_metrics.csv"))
summary = pd.read_csv(os.path.join(AGG_DIR, "4.2a_per_model_summary.csv"))

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.size': 12})

model_labels = {'Model_A': 'A: Universal', 'Model_B': 'B: Male', 'Model_C': 'C: Female', 'Model_D': 'D: Adenoma'}
df['Model_Label'] = df['Model'].map(model_labels)

# ── 1. Boxplot: AUC per fold, per model × algorithm ──
fig, ax = plt.subplots(figsize=(12, 6))
sns.boxplot(x='Model_Label', y='AUC', hue='Algorithm', data=df, ax=ax, palette='Set2')
sns.stripplot(x='Model_Label', y='AUC', hue='Algorithm', data=df, ax=ax,
              dodge=True, color='black', alpha=0.4, size=5, legend=False)
ax.axhline(0.5, ls='--', color='red', alpha=0.5, label='Random (0.5)')
ax.set_title('4.2a — Κατανομή AUC ανά LODO Fold (Median ± IQR)', fontsize=14, fontweight='bold')
ax.set_ylabel('AUC')
ax.set_xlabel('')
ax.set_ylim(0.3, 1.05)
ax.legend(loc='lower right')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '4.2a_auc_boxplot.png'), dpi=300)
plt.close()

# ── 2. Per-fold heatmap (RandomForest only) ──
rf = df[df['Algorithm'] == 'RandomForest'].pivot(index='Study', columns='Model_Label', values='AUC')
rf = rf.sort_index()

fig, ax = plt.subplots(figsize=(8, 7))
sns.heatmap(rf, annot=True, fmt='.2f', cmap='RdYlGn', vmin=0.4, vmax=1.0, 
            linewidths=0.5, ax=ax, cbar_kws={'label': 'AUC'})
ax.set_title('4.2a — AUC ανά Study × Model (Random Forest)', fontsize=14, fontweight='bold')
ax.set_ylabel('Study (LODO Fold)')
ax.set_xlabel('')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '4.2a_perfold_heatmap_RF.png'), dpi=300)
plt.close()

# ── 3. Per-fold heatmap (LASSO) ──
la = df[df['Algorithm'] == 'LASSO'].pivot(index='Study', columns='Model_Label', values='AUC')
la = la.sort_index()

fig, ax = plt.subplots(figsize=(8, 7))
sns.heatmap(la, annot=True, fmt='.2f', cmap='RdYlGn', vmin=0.4, vmax=1.0,
            linewidths=0.5, ax=ax, cbar_kws={'label': 'AUC'})
ax.set_title('4.2a — AUC ανά Study × Model (LASSO)', fontsize=14, fontweight='bold')
ax.set_ylabel('Study (LODO Fold)')
ax.set_xlabel('')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '4.2a_perfold_heatmap_LASSO.png'), dpi=300)
plt.close()

print(f"Τα διαγράμματα αποθηκεύτηκαν στο {FIG_DIR}")
