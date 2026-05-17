import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
import os
from sklearn.ensemble import RandomForestClassifier

OUTPUT_DIR = "OUTPUTS/phase_5/5.3_shap_visualizations"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Loading data...")
X_A = pd.read_csv("OUTPUTS/phase_4/4.0/data_Model_A.csv")
meta_A = pd.read_csv("OUTPUTS/phase_4/4.0/meta_Model_A.csv")

def clean_features(X):
    if 'Sample' in X.columns:
        X = X.drop(columns=['Sample'])
    return X.apply(pd.to_numeric, errors='coerce').fillna(0)

print("\nProcessing All with RandomForest...")
X_clean = clean_features(X_A)

if 'Age' in meta_A.columns:
    X_clean['Age'] = meta_A['Age'].fillna(meta_A['Age'].median())

y = (meta_A['Condition'] == 'CRC').astype(int)

# Normalization (Z-score)
X_mean = X_clean.mean()
X_std = X_clean.std().replace(0, 1)
X_norm = (X_clean - X_mean) / X_std

clf = RandomForestClassifier(n_estimators=500, class_weight='balanced', n_jobs=-1, random_state=42)
clf.fit(X_norm, y)

# TreeExplainer
explainer = shap.TreeExplainer(clf)
shap_values_raw = explainer.shap_values(X_norm, check_additivity=False)

if isinstance(shap_values_raw, np.ndarray) and shap_values_raw.ndim == 3:
    shap_values = shap_values_raw[:, :, 1]
elif isinstance(shap_values_raw, list):
    shap_values = shap_values_raw[1]
else:
    shap_values = shap_values_raw

# Generate Summary Plot (Beeswarm)
plt.figure(figsize=(10, 8))
shap.summary_plot(shap_values, X_norm, max_display=20, show=False)
plt.title("SHAP Summary Plot (Beeswarm) - All (RandomForest)")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/All_RF_shap_beeswarm.png", dpi=300)
plt.close()

print(f"Saved All_RF_shap_beeswarm.png to {OUTPUT_DIR}")
