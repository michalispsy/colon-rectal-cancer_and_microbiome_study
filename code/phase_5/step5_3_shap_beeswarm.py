import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
import os
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

OUTPUT_DIR = "OUTPUTS/phase_5/5.3_shap_visualizations"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load raw Phase 4 data
print("Loading data...")
X_A = pd.read_csv("OUTPUTS/phase_4/4.0/data_Model_A.csv")
meta_A = pd.read_csv("OUTPUTS/phase_4/4.0/meta_Model_A.csv")

X_B = pd.read_csv("OUTPUTS/phase_4/4.0/data_Model_B.csv")
meta_B = pd.read_csv("OUTPUTS/phase_4/4.0/meta_Model_B.csv")

X_C = pd.read_csv("OUTPUTS/phase_4/4.0/data_Model_C.csv")
meta_C = pd.read_csv("OUTPUTS/phase_4/4.0/meta_Model_C.csv")

models = {
    "All": (X_A, meta_A, "LASSO"),
    "OnlyMale": (X_B, meta_B, "RandomForest"),
    "OnlyFemale": (X_C, meta_C, "RandomForest")
}

def clean_features(X):
    if 'Sample' in X.columns:
        X = X.drop(columns=['Sample'])
    return X.apply(pd.to_numeric, errors='coerce').fillna(0)

for model_name, (X, meta, algo) in models.items():
    print(f"\nProcessing {model_name} ({algo})...")
    X_clean = clean_features(X)
    
    # Add age as covariate if possible, just like Phase 4
    if 'Age' in meta.columns:
        X_clean['Age'] = meta['Age'].fillna(meta['Age'].median())
    
    y = (meta['Condition'] == 'CRC').astype(int)
    
    # Normalization (Z-score)
    X_mean = X_clean.mean()
    X_std = X_clean.std().replace(0, 1)
    X_norm = (X_clean - X_mean) / X_std

    if algo == "LASSO":
        clf = LogisticRegression(penalty='l1', solver='liblinear', class_weight='balanced', max_iter=1000, random_state=42)
        clf.fit(X_norm, y)
        # LinearExplainer
        explainer = shap.LinearExplainer(clf, shap.sample(X_norm, 100))
        shap_values = explainer.shap_values(X_norm)
    else:
        clf = RandomForestClassifier(n_estimators=500, class_weight='balanced', n_jobs=-1, random_state=42)
        clf.fit(X_norm, y)
        # TreeExplainer
        explainer = shap.TreeExplainer(clf)
        shap_values_raw = explainer.shap_values(X_norm, check_additivity=False)
        shap_values = shap_values_raw[:, :, 1] if isinstance(shap_values_raw, np.ndarray) and shap_values_raw.ndim == 3 else shap_values_raw[1] if isinstance(shap_values_raw, list) else shap_values_raw

    # Generate Summary Plot (Beeswarm)
    plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, X_norm, max_display=20, show=False)
    plt.title(f"SHAP Summary Plot (Beeswarm) - {model_name}")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/{model_name}_shap_beeswarm.png", dpi=300)
    plt.close()
    
    print(f"Saved {model_name}_shap_beeswarm.png")

print("\nDone. Beeswarm plots saved in:", OUTPUT_DIR)
