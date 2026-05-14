import pandas as pd
import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score, confusion_matrix

INPUT_DIR = "OUTPUTS/phase_4/4.1/4.1_prep"
OUTPUT_DIR = "OUTPUTS/phase_4/4.3_advanced"

models = ['Model_A', 'Model_B', 'Model_C', 'Model_D']
algorithms = {
    'RandomForest': RandomForestClassifier(n_estimators=500, random_state=42, class_weight='balanced'),
    'LASSO': LogisticRegression(penalty='l1', solver='liblinear', random_state=42, class_weight='balanced', max_iter=1000)
}

def evaluate_model(y_true, y_pred_proba, y_pred_class):
    auc = roc_auc_score(y_true, y_pred_proba)
    auprc = average_precision_score(y_true, y_pred_proba)
    f1 = f1_score(y_true, y_pred_class)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred_class).ravel()
    return {'AUC': auc, 'AUPRC': auprc, 'F1': f1, 'TP': tp, 'FP': fp, 'TN': tn, 'FN': fn}

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    all_results = []
    fairness_results = []
    adenoma_predictions = []
    
    for model_name in models:
        model_dir = os.path.join(INPUT_DIR, model_name)
        if not os.path.exists(model_dir):
            continue
            
        studies = [d for d in os.listdir(model_dir) if os.path.isdir(os.path.join(model_dir, d))]
        
        print(f"\n{'='*40}\nAdvanced Analysis: {model_name}\n{'='*40}")
        
        for algo_name, clf in algorithms.items():
            print(f"--- Αλγόριθμος: {algo_name} ---")
            fold_metrics = []
            
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
                y_pred_class = clf.predict(X_test)
                
                metrics = evaluate_model(y_test, y_pred_proba, y_pred_class)
                metrics.update({'Study': study, 'Model': model_name, 'Algorithm': algo_name})
                fold_metrics.append(metrics)
                
                # Fairness Metrics & Sex-stratified AUC for Model A
                if model_name == 'Model_A' and 'Gender' in X_test.columns:
                    # Gender is standard scaled. values > 0 means 1 (Female) initially.
                    is_female = X_test['Gender'] > 0
                    
                    if is_female.sum() > 0 and (~is_female).sum() > 0:
                        y_test_f = y_test[is_female]
                        y_pred_proba_f = y_pred_proba[is_female]
                        y_pred_class_f = y_pred_class[is_female]
                        
                        y_test_m = y_test[~is_female]
                        y_pred_proba_m = y_pred_proba[~is_female]
                        y_pred_class_m = y_pred_class[~is_female]
                        
                        # Only calc AUC if both classes present
                        auc_f = roc_auc_score(y_test_f, y_pred_proba_f) if len(np.unique(y_test_f)) > 1 else np.nan
                        auc_m = roc_auc_score(y_test_m, y_pred_proba_m) if len(np.unique(y_test_m)) > 1 else np.nan
                        
                        # TPR (Recall)
                        tpr_f = np.sum((y_pred_class_f == 1) & (y_test_f == 1)) / np.sum(y_test_f == 1) if np.sum(y_test_f == 1) > 0 else np.nan
                        tpr_m = np.sum((y_pred_class_m == 1) & (y_test_m == 1)) / np.sum(y_test_m == 1) if np.sum(y_test_m == 1) > 0 else np.nan
                        equal_opportunity = abs(tpr_m - tpr_f) if not np.isnan(tpr_m) and not np.isnan(tpr_f) else np.nan
                        
                        # PPV (Precision)
                        ppv_f = np.sum((y_pred_class_f == 1) & (y_test_f == 1)) / np.sum(y_pred_class_f == 1) if np.sum(y_pred_class_f == 1) > 0 else np.nan
                        ppv_m = np.sum((y_pred_class_m == 1) & (y_test_m == 1)) / np.sum(y_pred_class_m == 1) if np.sum(y_pred_class_m == 1) > 0 else np.nan
                        
                        # Demographic Parity: P(y_hat=1 | group)
                        dp_f = np.mean(y_pred_class_f == 1)
                        dp_m = np.mean(y_pred_class_m == 1)
                        
                        fairness_results.append({
                            'Study': study, 'Algorithm': algo_name,
                            'AUC_F': auc_f, 'AUC_M': auc_m,
                            'TPR_F': tpr_f, 'TPR_M': tpr_m, 'Equal_Opportunity': equal_opportunity,
                            'PPV_F': ppv_f, 'PPV_M': ppv_m,
                            'DP_F': dp_f, 'DP_M': dp_m
                        })
                
                # Adenoma Scoring for Model A
                if model_name == 'Model_A':
                    ade_file = os.path.join(fold_dir, "X_adenoma.csv")
                    if os.path.exists(ade_file):
                        X_ade = pd.read_csv(ade_file, index_col=0)
                        # Reorder columns to match train exactly
                        X_ade = X_ade[X_train.columns]
                        ade_pred_proba = clf.predict_proba(X_ade)[:, 1]
                        for sample_id, prob in zip(X_ade.index, ade_pred_proba):
                            adenoma_predictions.append({
                                'Sample': sample_id,
                                'Study': study,
                                'Algorithm': algo_name,
                                'P_CRC': prob,
                                'Is_CRC_Like': 1 if prob > 0.5 else 0
                            })
                            
            all_results.extend(fold_metrics)
            
    # Αποθήκευση
    if len(all_results) > 0:
        pd.DataFrame(all_results).to_csv(os.path.join(OUTPUT_DIR, "advanced_metrics.csv"), index=False)
    if len(fairness_results) > 0:
        pd.DataFrame(fairness_results).to_csv(os.path.join(OUTPUT_DIR, "fairness_metrics.csv"), index=False)
    if len(adenoma_predictions) > 0:
        pd.DataFrame(adenoma_predictions).to_csv(os.path.join(OUTPUT_DIR, "adenoma_predictions.csv"), index=False)
        
    print("\nΟλοκληρώθηκε! Τα αποτελέσματα σώθηκαν στο", OUTPUT_DIR)

if __name__ == "__main__":
    main()
