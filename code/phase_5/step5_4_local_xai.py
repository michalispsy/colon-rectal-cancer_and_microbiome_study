import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
import os
from sklearn.ensemble import RandomForestClassifier

OUTPUT_DIR = "OUTPUTS/phase_5/5.4_local_xai"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Loading Universal Model A data...")
X_A = pd.read_csv("OUTPUTS/phase_4/4.0/data_Model_A.csv")
meta_A = pd.read_csv("OUTPUTS/phase_4/4.0/meta_Model_A.csv")

def clean_features(X):
    if 'Sample' in X.columns:
        X = X.drop(columns=['Sample'])
    return X.apply(pd.to_numeric, errors='coerce').fillna(0)

X_A_clean = clean_features(X_A)
if 'Age' in meta_A.columns:
    X_A_clean['Age'] = meta_A['Age'].fillna(meta_A['Age'].median())

y_A = (meta_A['Condition'] == 'CRC').astype(int)

# Pick a test study for out-of-sample evaluation (LODO simulation)
test_study = "WirbelJ_2018"
print(f"Setting up Out-Of-Sample evaluation. Test study: {test_study}")

train_idx = meta_A['Study'] != test_study
test_idx = meta_A['Study'] == test_study

X_train = X_A_clean[train_idx]
y_train = y_A[train_idx]

X_test = X_A_clean[test_idx].reset_index(drop=True)
y_test = y_A[test_idx].reset_index(drop=True)
meta_test = meta_A[test_idx].reset_index(drop=True)

# Train the Universal Model
print("Training Universal Model (RandomForest) on remaining studies...")
X_train_mean = X_train.mean()
X_train_std = X_train.std().replace(0, 1)

X_train_norm = (X_train - X_train_mean) / X_train_std
X_test_norm = (X_test - X_train_mean) / X_train_std

clf = RandomForestClassifier(n_estimators=500, class_weight='balanced', n_jobs=-1, random_state=42)
clf.fit(X_train_norm, y_train)

# Evaluate on TEST set
print("Predicting Probabilities on Test Set...")
probs_all = clf.predict_proba(X_test_norm)[:, 1]
meta_test['P_CRC'] = probs_all

# Separate into CRC and Control
crc_mask = meta_test['Condition'] == 'CRC'
ctrl_mask = meta_test['Condition'] == 'Control'

# Find the 4 specific clinical cases
crc_indices = np.where(crc_mask)[0]
probs_crc = probs_all[crc_indices]
idx_crc_high = crc_indices[np.argmax(probs_crc)] # True Positive
idx_crc_low = crc_indices[np.argmin(probs_crc)]  # False Negative

ctrl_indices = np.where(ctrl_mask)[0]
probs_ctrl = probs_all[ctrl_indices]
idx_ctrl_high = ctrl_indices[np.argmax(probs_ctrl)] # False Positive
idx_ctrl_low = ctrl_indices[np.argmin(probs_ctrl)]  # True Negative

cases = [
    (idx_crc_high, "true_positive", "True Positive (CRC detected)"),
    (idx_crc_low, "false_negative", "False Negative (CRC missed)"),
    (idx_ctrl_high, "false_positive", "False Positive (Healthy misclassified as CRC)"),
    (idx_ctrl_low, "true_negative", "True Negative (Healthy confirmed)")
]

print("Calculating SHAP values for the Test Set...")
explainer = shap.TreeExplainer(clf)
shap_values_raw = explainer(X_test_norm)

if len(shap_values_raw.shape) == 3:
    shap_values_test = shap_values_raw[:, :, 1]
else:
    shap_values_test = shap_values_raw

def plot_waterfall(idx, suffix, label):
    sample_info = meta_test.iloc[idx]
    plt.figure(figsize=(26, 12))
    shap.waterfall_plot(shap_values_test[idx], max_display=15, show=False)
    
    title = (f"Patient: {sample_info['Sample']} | Study: {sample_info['Study']} | Actual: {sample_info['Condition']}\n"
             f"RandomForest P(CRC): {sample_info['P_CRC']:.2f} -> {label}")
    plt.title(title, pad=20, fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    filepath = f"{OUTPUT_DIR}/patient_{suffix}_waterfall.png"
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filepath}")

print("Generating waterfall plots for the 4 key cases...")
for idx, suffix, label in cases:
    plot_waterfall(idx, suffix, label)

print(f"Done. Plots saved in {OUTPUT_DIR}")
