import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.ensemble import RandomForestClassifier

PREP_DIR = "OUTPUTS/phase_4/4.1/4.1_prep"
INPUT_DIR = "OUTPUTS/phase_4/4.1/4.1_results"
AGG_DIR = "OUTPUTS/phase_4/4.2_aggregation"
FIG_DIR = os.path.join(AGG_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.size': 12})

# ==============================================================
# 4.2d: Adenoma Positioning Summary (Model A)
# ==============================================================
print("Εκτέλεση 4.2d: Adenoma Positioning...")

# We only focus on RandomForest for this part
model_dir = os.path.join(PREP_DIR, 'Model_A')
studies = [d for d in os.listdir(model_dir) if os.path.isdir(os.path.join(model_dir, d))]

clf = RandomForestClassifier(n_estimators=500, random_state=42, class_weight='balanced')
test_predictions = []

for study in studies:
    fold_dir = os.path.join(model_dir, study)
    try:
        X_train = pd.read_csv(os.path.join(fold_dir, "X_train.csv"), index_col=0)
        y_train = pd.read_csv(os.path.join(fold_dir, "y_train.csv"), index_col=0)['Label']
        X_test = pd.read_csv(os.path.join(fold_dir, "X_test.csv"), index_col=0)
        y_test = pd.read_csv(os.path.join(fold_dir, "y_test.csv"), index_col=0)['Label']
    except FileNotFoundError:
        continue
        
    clf.fit(X_train, y_train)
    y_pred_proba = clf.predict_proba(X_test)[:, 1]
    
    for yt, prob in zip(y_test, y_pred_proba):
        test_predictions.append({
            'Study': study,
            'Condition': 'CRC' if yt == 1 else 'Control',
            'P_CRC': prob
        })

df_test = pd.DataFrame(test_predictions)

# Load Adenoma predictions
ade_preds = pd.read_csv(os.path.join(INPUT_DIR, "adenoma_predictions.csv"))
df_ade = ade_preds[ade_preds['Algorithm'] == 'RandomForest'].copy()
df_ade['Condition'] = 'Adenoma'
df_ade = df_ade[['Study', 'Condition', 'P_CRC']]

df_all = pd.concat([df_test, df_ade], ignore_index=True)

# Order for plots
order = ['Control', 'Adenoma', 'CRC']
colors = {'Control': '#2ecc71', 'Adenoma': '#9b59b6', 'CRC': '#e74c3c'}

# Plot 4.2d: KDE Distribution
plt.figure(figsize=(8, 5))
for cond in order:
    sns.kdeplot(data=df_all[df_all['Condition'] == cond], x='P_CRC', 
                fill=True, label=cond, color=colors[cond], alpha=0.5)
plt.axvline(0.5, ls='--', color='black', alpha=0.8, label='Threshold (0.5)')
plt.title('4.2d — Κατανομή P(CRC): Control vs Adenoma vs CRC')
plt.xlabel('Predicted P(CRC) (Model A - RF)')
plt.ylabel('Πυκνότητα (Density)')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '4.2d_adenoma_distribution.png'), dpi=300)
plt.close()

# Plot 4.2d: Boxplot
plt.figure(figsize=(6, 5))
sns.boxplot(x='Condition', y='P_CRC', data=df_all, order=order, palette=colors)
sns.stripplot(x='Condition', y='P_CRC', data=df_all, order=order, color='black', alpha=0.1, jitter=True)
plt.axhline(0.5, ls='--', color='black', alpha=0.8)
plt.title('4.2d — Risk Scores ανά Κατάσταση')
plt.ylabel('Predicted P(CRC)')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '4.2d_adenoma_boxplot.png'), dpi=300)
plt.close()

# Stats for report
pct_crc_like = (df_ade['P_CRC'] > 0.5).mean() * 100
mean_p_ade = df_ade['P_CRC'].mean()
mean_p_ctrl = df_test[df_test['Condition'] == 'Control']['P_CRC'].mean()
mean_p_crc = df_test[df_test['Condition'] == 'CRC']['P_CRC'].mean()


# ==============================================================
# 4.2e: Sex-Specific Comparison
# ==============================================================
print("Εκτέλεση 4.2e: Sex-Specific Comparison...")

adv_df = pd.read_csv(os.path.join(INPUT_DIR, "advanced_metrics.csv"))
fair_df = pd.read_csv(os.path.join(INPUT_DIR, "fairness_metrics.csv"))

# Filter for RF
adv_rf = adv_df[adv_df['Algorithm'] == 'RandomForest'].copy()
fair_rf = fair_df[fair_df['Algorithm'] == 'RandomForest'].copy()

# Model B (Males only)
df_b = adv_rf[adv_rf['Model'] == 'Model_B'][['Study', 'AUC']].rename(columns={'AUC': 'AUC_ModelB'})
# Model A on Males
df_a_m = fair_rf[['Study', 'AUC_M']].rename(columns={'AUC_M': 'AUC_ModelA_on_Males'})
comp_males = pd.merge(df_b, df_a_m, on='Study').dropna()

# Model C (Females only)
df_c = adv_rf[adv_rf['Model'] == 'Model_C'][['Study', 'AUC']].rename(columns={'AUC': 'AUC_ModelC'})
# Model A on Females
df_a_f = fair_rf[['Study', 'AUC_F']].rename(columns={'AUC_F': 'AUC_ModelA_on_Females'})
comp_females = pd.merge(df_c, df_a_f, on='Study').dropna()

# Plot 4.2e: Males (B vs A_males)
plt.figure(figsize=(6, 5))
plot_m = pd.melt(comp_males, id_vars=['Study'], value_vars=['AUC_ModelA_on_Males', 'AUC_ModelB'],
                 var_name='Model_Type', value_name='AUC')
plot_m['Model_Type'] = plot_m['Model_Type'].map({'AUC_ModelA_on_Males': 'Universal (A)', 'AUC_ModelB': 'Sex-Specific (B)'})

sns.boxplot(x='Model_Type', y='AUC', data=plot_m, palette=['#cccccc', '#66b3ff'])
for study in comp_males['Study']:
    row = comp_males[comp_males['Study'] == study].iloc[0]
    plt.plot([0, 1], [row['AUC_ModelA_on_Males'], row['AUC_ModelB']], color='gray', alpha=0.5, marker='o')

plt.title('4.2e — Άνδρες: Universal vs Sex-Specific Model')
plt.ylabel('AUC')
plt.xlabel('')
plt.ylim(0.4, 1.05)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '4.2e_males_comparison.png'), dpi=300)
plt.close()

# Plot 4.2e: Females (C vs A_females)
plt.figure(figsize=(6, 5))
plot_f = pd.melt(comp_females, id_vars=['Study'], value_vars=['AUC_ModelA_on_Females', 'AUC_ModelC'],
                 var_name='Model_Type', value_name='AUC')
plot_f['Model_Type'] = plot_f['Model_Type'].map({'AUC_ModelA_on_Females': 'Universal (A)', 'AUC_ModelC': 'Sex-Specific (C)'})

sns.boxplot(x='Model_Type', y='AUC', data=plot_f, palette=['#cccccc', '#ff9999'])
for study in comp_females['Study']:
    row = comp_females[comp_females['Study'] == study].iloc[0]
    plt.plot([0, 1], [row['AUC_ModelA_on_Females'], row['AUC_ModelC']], color='gray', alpha=0.5, marker='o')

plt.title('4.2e — Γυναίκες: Universal vs Sex-Specific Model')
plt.ylabel('AUC')
plt.xlabel('')
plt.ylim(0.4, 1.05)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, '4.2e_females_comparison.png'), dpi=300)
plt.close()

# Medians for report
med_A_m = comp_males['AUC_ModelA_on_Males'].median()
med_B = comp_males['AUC_ModelB'].median()

med_A_f = comp_females['AUC_ModelA_on_Females'].median()
med_C = comp_females['AUC_ModelC'].median()

# Write stats to text file to easily read them in python
with open(os.path.join(AGG_DIR, "4.2de_stats.txt"), "w") as f:
    f.write(f"mean_p_ctrl={mean_p_ctrl}\n")
    f.write(f"mean_p_ade={mean_p_ade}\n")
    f.write(f"mean_p_crc={mean_p_crc}\n")
    f.write(f"pct_crc_like={pct_crc_like}\n")
    f.write(f"med_A_m={med_A_m}\n")
    f.write(f"med_B={med_B}\n")
    f.write(f"med_A_f={med_A_f}\n")
    f.write(f"med_C={med_C}\n")

print("Ολοκληρώθηκε!")
