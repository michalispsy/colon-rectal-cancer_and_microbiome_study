import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import roc_auc_score

OUTPUT_DIR = "OUTPUTS/phase_5/5.5_subgroup_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Loading data for Subgroup Analysis...")
X = pd.read_csv("OUTPUTS/phase_4/4.0/data_Model_A.csv")
meta = pd.read_csv("OUTPUTS/phase_4/4.0/meta_Model_A.csv")

def clean_features(X):
    if 'Sample' in X.columns:
        X = X.drop(columns=['Sample'])
    return X.apply(pd.to_numeric, errors='coerce').fillna(0)

X_clean = clean_features(X)
if 'Age' in meta.columns:
    X_clean['Age'] = meta['Age'].fillna(meta['Age'].median())
    # Create Age Groups
    meta['Age_Group'] = pd.cut(meta['Age'], bins=[0, 50, 65, 100], labels=['<50', '50-65', '>65'])
else:
    meta['Age_Group'] = 'Unknown'

y = (meta['Condition'] == 'CRC').astype(int)
groups = meta['Study']

# Replicate LODO CV to get unbiased predictions
print("Running Leave-One-Dataset-Out (LODO) Cross-Validation...")
logo = LeaveOneGroupOut()
meta['P_CRC_CV'] = np.nan

for train_idx, test_idx in logo.split(X_clean, y, groups):
    X_train, X_test = X_clean.iloc[train_idx], X_clean.iloc[test_idx]
    y_train = y.iloc[train_idx]
    
    # Normalize using train
    X_mean = X_train.mean()
    X_std = X_train.std().replace(0, 1)
    
    X_train_norm = (X_train - X_mean) / X_std
    X_test_norm = (X_test - X_mean) / X_std
    
    clf = LogisticRegression(penalty='l1', solver='liblinear', class_weight='balanced', max_iter=1000, random_state=42)
    clf.fit(X_train_norm, y_train)
    
    meta.loc[test_idx, 'P_CRC_CV'] = clf.predict_proba(X_test_norm)[:, 1]

meta['y_true'] = y

def safe_auc(y_true, y_score):
    if len(np.unique(y_true)) < 2:
        return np.nan
    return roc_auc_score(y_true, y_score)

# Evaluate Subgroups
print("\nEvaluating Subgroups...")
subgroups = ['Age_Group', 'Gender', 'Study']
results = []

for sg in subgroups:
    if sg not in meta.columns:
        continue
    for val in meta[sg].dropna().unique():
        subset = meta[meta[sg] == val]
        auc = safe_auc(subset['y_true'], subset['P_CRC_CV'])
        n_samples = len(subset)
        results.append({'Subgroup': sg, 'Value': val, 'N_Samples': n_samples, 'AUC': auc})

df_res = pd.DataFrame(results).dropna(subset=['AUC']).sort_values(['Subgroup', 'AUC'], ascending=[True, False])
df_res.to_csv(f"{OUTPUT_DIR}/subgroup_auc_results.csv", index=False)

# Plot Results
for sg in df_res['Subgroup'].unique():
    subset = df_res[df_res['Subgroup'] == sg]
    plt.figure(figsize=(8, 5))
    plt.bar(subset['Value'].astype(str), subset['AUC'], color='teal')
    plt.axhline(0.5, color='red', linestyle='--')
    plt.ylabel('AUC')
    plt.title(f'AUC by {sg}')
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/auc_by_{sg.lower()}.png", dpi=300)
    plt.close()

print(f"Done. Subgroup analysis saved in {OUTPUT_DIR}")
