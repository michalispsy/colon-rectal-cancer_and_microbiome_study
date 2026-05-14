import pandas as pd
import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score, precision_score, recall_score

INPUT_DIR = "OUTPUTS/phase_4/4.1/4.1_prep"
OUTPUT_DIR = "OUTPUTS/phase_4/4.2_training"

models = ['Model_A', 'Model_B', 'Model_C', 'Model_D']
algorithms = {
    'RandomForest': RandomForestClassifier(n_estimators=500, random_state=42, class_weight='balanced'),
    'LASSO': LogisticRegression(penalty='l1', solver='liblinear', random_state=42, class_weight='balanced', max_iter=1000)
}

def evaluate_model(y_true, y_pred_proba, y_pred_class):
    auc = roc_auc_score(y_true, y_pred_proba)
    auprc = average_precision_score(y_true, y_pred_proba)
    f1 = f1_score(y_true, y_pred_class)
    return {'AUC': auc, 'AUPRC': auprc, 'F1': f1}

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    all_results = []
    
    for model_name in models:
        model_dir = os.path.join(INPUT_DIR, model_name)
        if not os.path.exists(model_dir):
            continue
            
        studies = [d for d in os.listdir(model_dir) if os.path.isdir(os.path.join(model_dir, d))]
        
        print(f"\n{'='*40}\nΕκπαίδευση {model_name}\n{'='*40}")
        
        for algo_name, clf in algorithms.items():
            print(f"--- Αλγόριθμος: {algo_name} ---")
            
            fold_metrics = []
            
            for study in studies:
                fold_dir = os.path.join(model_dir, study)
                
                # Φόρτωση δεδομένων
                try:
                    X_train = pd.read_csv(os.path.join(fold_dir, "X_train.csv"), index_col=0)
                    y_train = pd.read_csv(os.path.join(fold_dir, "y_train.csv"), index_col=0)['Label']
                    X_test = pd.read_csv(os.path.join(fold_dir, "X_test.csv"), index_col=0)
                    y_test = pd.read_csv(os.path.join(fold_dir, "y_test.csv"), index_col=0)['Label']
                except FileNotFoundError:
                    continue
                    
                # Εκπαίδευση
                clf.fit(X_train, y_train)
                
                # Πρόβλεψη
                y_pred_proba = clf.predict_proba(X_test)[:, 1]
                y_pred_class = clf.predict(X_test)
                
                # Αξιολόγηση
                metrics = evaluate_model(y_test, y_pred_proba, y_pred_class)
                metrics['Study'] = study
                metrics['Model'] = model_name
                metrics['Algorithm'] = algo_name
                metrics['Train_Size'] = len(X_train)
                metrics['Test_Size'] = len(X_test)
                
                fold_metrics.append(metrics)
                
            if len(fold_metrics) > 0:
                df_metrics = pd.DataFrame(fold_metrics)
                mean_auc = df_metrics['AUC'].mean()
                mean_f1 = df_metrics['F1'].mean()
                print(f"Μέσο AUC: {mean_auc:.3f} | Μέσο F1: {mean_f1:.3f}")
                
                all_results.extend(fold_metrics)

    # Αποθήκευση συνολικών αποτελεσμάτων
    if len(all_results) > 0:
        final_df = pd.DataFrame(all_results)
        final_df.to_csv(os.path.join(OUTPUT_DIR, "all_models_cv_results.csv"), index=False)
        
        # Summary table (Mέσοι όροι)
        summary = final_df.groupby(['Model', 'Algorithm'])[['AUC', 'AUPRC', 'F1']].mean().reset_index()
        summary.to_csv(os.path.join(OUTPUT_DIR, "models_summary.csv"), index=False)
        print("\nΤα αποτελέσματα αποθηκεύτηκαν στο OUTPUTS/phase_4/4.2_training/")

if __name__ == "__main__":
    main()
