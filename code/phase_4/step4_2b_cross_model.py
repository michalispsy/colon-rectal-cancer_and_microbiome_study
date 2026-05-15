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

# ── Cross-Model Comparison Table (AUC only, best algo per model) ──
auc_summary = summary[summary['Metric'] == 'AUC'].copy()

# For each model, pick the best algorithm by Median AUC
best_per_model = auc_summary.loc[auc_summary.groupby('Model')['Median'].idxmax()].copy()
best_per_model['IQR_str'] = best_per_model.apply(lambda r: f"[{r['IQR_25']:.3f} – {r['IQR_75']:.3f}]", axis=1)
best_per_model['Median_str'] = best_per_model['Median'].apply(lambda x: f"{x:.3f}")

model_labels = {'Model_A': 'A: Universal (CRC vs Ctrl)', 
                'Model_B': 'B: Male-Only (CRC vs Ctrl)',
                'Model_C': 'C: Female-Only (CRC vs Ctrl)', 
                'Model_D': 'D: Adenoma vs Control'}
best_per_model['Model_Label'] = best_per_model['Model'].map(model_labels)

cross_table = best_per_model[['Model_Label', 'Algorithm', 'Median_str', 'IQR_str', 'N_folds']].copy()
cross_table.columns = ['Model', 'Best Algorithm', 'Median AUC', 'IQR', 'N folds']
cross_table.to_csv(os.path.join(AGG_DIR, "4.2b_cross_model_comparison.csv"), index=False)

print("\n" + "="*70)
print("4.2b  Cross-Model Comparison")
print("="*70)
print(cross_table.to_string(index=False))

# ── Visualization 1: Dot plot with error bars (Median ± IQR) ──
fig, ax = plt.subplots(figsize=(10, 5))

models_order = ['Model_A', 'Model_B', 'Model_C', 'Model_D']
colors_algo = {'RandomForest': '#2196F3', 'LASSO': '#FF9800'}
labels_model = {'Model_A': 'A: Universal', 'Model_B': 'B: Male', 'Model_C': 'C: Female', 'Model_D': 'D: Adenoma'}

y_pos = np.arange(len(models_order))
offset = 0.15

for i, algo in enumerate(['RandomForest', 'LASSO']):
    medians, lows, highs = [], [], []
    for model in models_order:
        row = auc_summary[(auc_summary['Model'] == model) & (auc_summary['Algorithm'] == algo)]
        if len(row) == 0:
            medians.append(np.nan); lows.append(0); highs.append(0)
        else:
            r = row.iloc[0]
            medians.append(r['Median'])
            lows.append(r['Median'] - r['IQR_25'])
            highs.append(r['IQR_75'] - r['Median'])
    
    ax.errorbar(medians, y_pos + (i-0.5)*offset*2, xerr=[lows, highs],
                fmt='o', markersize=10, capsize=5, capthick=2, linewidth=2,
                color=colors_algo[algo], label=algo)

ax.set_yticks(y_pos)
ax.set_yticklabels([labels_model[m] for m in models_order])
ax.axvline(0.5, ls='--', color='red', alpha=0.4, label='Random (0.5)')
ax.set_xlabel('Median AUC (LODO)', fontsize=12)
ax.set_title('4.2b — Cross-Model Comparison (Median AUC ± IQR)', fontsize=14, fontweight='bold')
ax.set_xlim(0.35, 1.0)
ax.legend(loc='lower right')
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '4.2b_cross_model_dotplot.png'), dpi=300)
plt.close()

# ── Visualization 2: Paired bar chart RF vs LASSO ──
fig, ax = plt.subplots(figsize=(10, 5))

x = np.arange(len(models_order))
width = 0.35

rf_medians = [auc_summary[(auc_summary['Model']==m) & (auc_summary['Algorithm']=='RandomForest')]['Median'].values[0] for m in models_order]
la_medians = [auc_summary[(auc_summary['Model']==m) & (auc_summary['Algorithm']=='LASSO')]['Median'].values[0] for m in models_order]

bars1 = ax.bar(x - width/2, rf_medians, width, label='Random Forest', color='#2196F3', alpha=0.85)
bars2 = ax.bar(x + width/2, la_medians, width, label='LASSO', color='#FF9800', alpha=0.85)

ax.bar_label(bars1, fmt='%.3f', padding=3)
ax.bar_label(bars2, fmt='%.3f', padding=3)

ax.set_xticks(x)
ax.set_xticklabels([labels_model[m] for m in models_order])
ax.axhline(0.5, ls='--', color='red', alpha=0.4)
ax.set_ylabel('Median AUC')
ax.set_title('4.2b — RF vs LASSO: Median AUC per Model', fontsize=14, fontweight='bold')
ax.set_ylim(0.35, 0.95)
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '4.2b_rf_vs_lasso_bars.png'), dpi=300)
plt.close()

print(f"\nΔιαγράμματα αποθηκεύτηκαν στο {FIG_DIR}")
