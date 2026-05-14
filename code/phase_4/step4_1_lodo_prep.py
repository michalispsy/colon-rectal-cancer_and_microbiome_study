import pandas as pd
import numpy as np
import os
from combat.pycombat import pycombat
from sklearn.preprocessing import StandardScaler

# Ορισμός διαδρομών
INPUT_DIR = "OUTPUTS/phase_4/4.0"
OUTPUT_DIR = "OUTPUTS/phase_4/4.1_prep"

models = {
    'Model_A': {'covariates': ['Condition', 'Gender', 'Age']},
    'Model_B': {'covariates': ['Condition', 'Age']},
    'Model_C': {'covariates': ['Condition', 'Age']},
    'Model_D': {'covariates': ['Condition', 'Gender', 'Age']}
}

def process_model(model_name, covariates):
    print(f"\n{'='*50}\nΠροετοιμασία: {model_name}\n{'='*50}")
    
    # Φόρτωση δεδομένων
    data_path = os.path.join(INPUT_DIR, f"data_{model_name}.csv")
    meta_path = os.path.join(INPUT_DIR, f"meta_{model_name}.csv")
    
    if not os.path.exists(data_path) or not os.path.exists(meta_path):
        print(f"Σφάλμα: Δεν βρέθηκαν τα αρχεία για το {model_name}.")
        return
        
    X = pd.read_csv(data_path)
    meta = pd.read_csv(meta_path)
    
    X_ade, meta_ade = None, None
    if model_name == 'Model_A':
        ade_data_path = os.path.join(INPUT_DIR, "data_Adenoma.csv")
        ade_meta_path = os.path.join(INPUT_DIR, "meta_Adenoma.csv")
        if os.path.exists(ade_data_path) and os.path.exists(ade_meta_path):
            X_ade = pd.read_csv(ade_data_path)
            meta_ade = pd.read_csv(ade_meta_path)
    
    # Διαχωρισμός στηλών
    pa_cols = [c for c in X.columns if c.startswith('PA_')]
    rclr_cols = [c for c in X.columns if not c.startswith('PA_') and c != 'Sample']
    
    studies = sorted(meta['Study'].unique())
    
    model_out_dir = os.path.join(OUTPUT_DIR, model_name)
    os.makedirs(model_out_dir, exist_ok=True)
    
    # LODO Loop
    for test_study in studies:
        print(f"--- Fold: {test_study} ---")
        
        # Εύρεση δεικτών
        test_idx = meta['Study'] == test_study
        train_idx = ~test_idx
        
        meta_train = meta[train_idx]
        meta_test = meta[test_idx]
        
        # Αν το test set δεν έχει αρκετά δείγματα για την κλάση στόχο, το κάνουμε skip
        target_counts = meta_test['Condition'].value_counts()
        if len(target_counts) < 2 or target_counts.min() < 5:
            print(f"  -> Παράλειψη (πολύ μικρό test set: {target_counts.to_dict()})")
            continue
            
        X_train = X[train_idx]
        X_test = X[test_idx]
        
        fold_out_dir = os.path.join(model_out_dir, test_study)
        os.makedirs(fold_out_dir, exist_ok=True)
        
        # BATCH CORRECTION με pyComBat στο Training Set (μόνο rCLR)
        X_train_rclr = X_train[['Sample'] + rclr_cols].set_index('Sample')
        meta_train_indexed = meta_train.set_index('Sample')
        
        # pycombat expects: rows=features, cols=samples
        data_for_combat = X_train_rclr.T
        batch = meta_train_indexed['Study'].tolist()
        
        # Covariates encoding
        mod = []
        if 'Condition' in covariates:
            cond_encoded = (meta_train_indexed['Condition'] != 'Control').astype(int).tolist()
            mod.append(cond_encoded)
        if 'Gender' in covariates:
            gender_encoded = (meta_train_indexed['Gender'] == 'Female').astype(int).tolist()
            mod.append(gender_encoded)
        if 'Age' in covariates:
            age_encoded = meta_train_indexed['Age'].fillna(meta_train_indexed['Age'].mean()).tolist()
            mod.append(age_encoded)
            
        try:
            print(f"  -> pyComBat σε {len(rclr_cols)} features...")
            # Εκτέλεση pyComBat
            corrected_rclr_T = pycombat(data_for_combat, batch, mod=mod)
            X_train_rclr_corr = corrected_rclr_T.T
        except Exception as e:
            print(f"  -> Σφάλμα pyComBat: {e}")
            continue
            
        # Συνένωση PA και Corrected rCLR
        X_train_pa = X_train[['Sample'] + pa_cols].set_index('Sample')
        
        # Συνένωση και προσθήκη metadata (Age, Gender)
        X_train_final = pd.concat([X_train_pa, X_train_rclr_corr], axis=1)
        
        if 'Gender' in covariates:
            X_train_final['Gender'] = (meta_train_indexed['Gender'] == 'Female').astype(int)
        if 'Age' in covariates:
            X_train_final['Age'] = meta_train_indexed['Age'].fillna(meta_train_indexed['Age'].mean())
            
        # Για το Test Set (ΧΩΡΙΣ Batch Correction, χρησιμοποιούμε τα αρχικά rCLR)
        X_test_pa = X_test[['Sample'] + pa_cols].set_index('Sample')
        X_test_rclr = X_test[['Sample'] + rclr_cols].set_index('Sample')
        meta_test_indexed = meta_test.set_index('Sample')
        
        X_test_final = pd.concat([X_test_pa, X_test_rclr], axis=1)
        if 'Gender' in covariates:
            X_test_final['Gender'] = (meta_test_indexed['Gender'] == 'Female').astype(int)
        if 'Age' in covariates:
            X_test_final['Age'] = meta_test_indexed['Age'].fillna(meta_test_indexed['Age'].mean())
            
        # Εξαγωγή Labels (0 = Control, 1 = Disease)
        y_train = (meta_train_indexed['Condition'] != 'Control').astype(int)
        y_test = (meta_test_indexed['Condition'] != 'Control').astype(int)
        
        # SCALING (Z-score)
        # Κάνουμε fit ΜΟΝΟ στο train set, transform σε train και test
        scaler = StandardScaler()
        
        X_train_scaled = pd.DataFrame(scaler.fit_transform(X_train_final), 
                                      index=X_train_final.index, 
                                      columns=X_train_final.columns)
        X_test_scaled = pd.DataFrame(scaler.transform(X_test_final), 
                                     index=X_test_final.index, 
                                     columns=X_test_final.columns)
        
        # Αποθήκευση τελικών έτοιμων matrices για το Fold
        X_train_scaled.to_csv(os.path.join(fold_out_dir, "X_train.csv"))
        y_train.to_frame(name="Label").to_csv(os.path.join(fold_out_dir, "y_train.csv"))
        X_test_scaled.to_csv(os.path.join(fold_out_dir, "X_test.csv"))
        y_test.to_frame(name="Label").to_csv(os.path.join(fold_out_dir, "y_test.csv"))
        
        # Ειδικά για Model_A: Προετοιμασία Adenoma samples της ίδιας μελέτης (αν υπάρχουν)
        if model_name == 'Model_A' and X_ade is not None:
            ade_idx = meta_ade['Study'] == test_study
            if ade_idx.sum() > 0:
                X_ade_study = X_ade[ade_idx]
                meta_ade_study = meta_ade[ade_idx]
                
                X_ade_pa = X_ade_study[['Sample'] + pa_cols].set_index('Sample')
                X_ade_rclr = X_ade_study[['Sample'] + rclr_cols].set_index('Sample')
                meta_ade_indexed = meta_ade_study.set_index('Sample')
                
                X_ade_final = pd.concat([X_ade_pa, X_ade_rclr], axis=1)
                X_ade_final['Gender'] = (meta_ade_indexed['Gender'] == 'Female').astype(int)
                X_ade_final['Age'] = meta_ade_indexed['Age'].fillna(meta_train_indexed['Age'].mean())
                
                X_ade_scaled = pd.DataFrame(scaler.transform(X_ade_final),
                                            index=X_ade_final.index,
                                            columns=X_ade_final.columns)
                y_ade = (meta_ade_indexed['Condition'] != 'Control').astype(int) # Όλα θα είναι 1 ουσιαστικά, αλλά δε μας νοιάζει το label
                
                X_ade_scaled.to_csv(os.path.join(fold_out_dir, "X_adenoma.csv"))
                y_ade.to_frame(name="Label").to_csv(os.path.join(fold_out_dir, "y_adenoma.csv"))
                print(f"  -> Προετοιμάστηκαν {len(X_ade_scaled)} δείγματα Adenoma")
        
        print(f"  -> Ολοκληρώθηκε: Train={len(X_train_scaled)}, Test={len(X_test_scaled)}")

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for model_name, config in models.items():
        process_model(model_name, config['covariates'])
        
    print("\nΗ προετοιμασία των LODO datasets ολοκληρώθηκε επιτυχώς με pyComBat!")

if __name__ == "__main__":
    main()
