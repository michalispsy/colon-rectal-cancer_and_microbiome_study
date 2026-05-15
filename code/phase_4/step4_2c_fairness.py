import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

INPUT_DIR = "OUTPUTS/phase_4/4.1/4.1_results"
AGG_DIR = "OUTPUTS/phase_4/4.2_aggregation"
FIG_DIR = os.path.join(AGG_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

df = pd.read_csv(os.path.join(INPUT_DIR, "fairness_metrics.csv"))

# We only focus on RandomForest for the main fairness report, or calculate for both?
# Let's group by Algorithm and calculate medians
df['Predictive_Parity_Gap'] = (df['PPV_M'] - df['PPV_F']).abs()
df['Demographic_Parity_Gap'] = (df['DP_M'] - df['DP_F']).abs()
# Equal_Opportunity is already calculated as |TPR_M - TPR_F|

metrics = ['AUC_F', 'AUC_M', 'Equal_Opportunity', 'Predictive_Parity_Gap', 'Demographic_Parity_Gap']

summary_rows = []
for algo in ['RandomForest', 'LASSO']:
    sub = df[df['Algorithm'] == algo]
    if len(sub) == 0: continue
    
    row = {'Algorithm': algo}
    for m in metrics:
        row[f'{m}_Median'] = sub[m].median()
        row[f'{m}_IQR_25'] = sub[m].quantile(0.25)
        row[f'{m}_IQR_75'] = sub[m].quantile(0.75)
    summary_rows.append(row)

summary = pd.DataFrame(summary_rows)
summary.to_csv(os.path.join(AGG_DIR, "4.2c_fairness_summary.csv"), index=False)

# ── Visualization: Fairness Gaps (RF) ──
rf_data = df[df['Algorithm'] == 'RandomForest'].copy()

# Plot 1: AUC_M vs AUC_F Paired Boxplot
plt.figure(figsize=(6, 5))
plot_df = pd.melt(rf_data[['AUC_M', 'AUC_F']], var_name='Sex', value_name='AUC')
plot_df['Sex'] = plot_df['Sex'].map({'AUC_M': 'Άνδρες', 'AUC_F': 'Γυναίκες'})
sns.boxplot(x='Sex', y='AUC', data=plot_df, palette=['#66b3ff', '#ff9999'])
sns.stripplot(x='Sex', y='AUC', data=plot_df, color='black', alpha=0.5, jitter=True)
plt.title('4.2c — Sex-Stratified AUC (Model A - RF)')
plt.ylim(0.4, 1.05)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '4.2c_fairness_auc.png'), dpi=300)
plt.close()

# Plot 2: Fairness Gaps Boxplot
plt.figure(figsize=(8, 5))
gaps_df = pd.melt(rf_data[['Demographic_Parity_Gap', 'Equal_Opportunity', 'Predictive_Parity_Gap']], 
                  var_name='Metric', value_name='Gap')
gaps_df['Metric'] = gaps_df['Metric'].map({
    'Demographic_Parity_Gap': 'Demographic Parity\n(|DP_M - DP_F|)',
    'Equal_Opportunity': 'Equal Opportunity\n(|TPR_M - TPR_F|)',
    'Predictive_Parity_Gap': 'Predictive Parity\n(|PPV_M - PPV_F|)'
})
sns.boxplot(x='Metric', y='Gap', data=gaps_df, palette='Set3')
sns.stripplot(x='Metric', y='Gap', data=gaps_df, color='black', alpha=0.5, jitter=True)
plt.axhline(0.1, ls='--', color='red', alpha=0.5, label='Ιδανικό όριο (10%)')
plt.title('4.2c — Fairness Gaps (Model A - RF)')
plt.ylabel('Απόλυτη Διαφορά (Gap)')
plt.xlabel('')
plt.ylim(0, 0.5)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '4.2c_fairness_gaps.png'), dpi=300)
plt.close()

print(f"Αποτελέσματα Fairness (4.2c) αποθηκεύτηκαν στο {AGG_DIR}")
