import pandas as pd
import os

# Define paths
X_EXPLORE_PATH = "data/crc_study_final_data/species_level/X_explore.csv"
METADATA_PATH = "data/crc_study_final_data/species_level/metadata.csv"
OUTPUT_DIR = "OUTPUTS/phase_4/4.0"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Φόρτωση δεδομένων...")
X = pd.read_csv(X_EXPLORE_PATH)
meta = pd.read_csv(METADATA_PATH)

# Ensure sample alignment
if not X['Sample'].equals(meta['Sample']):
    print("ΠΡΟΕΙΔΟΠΟΙΗΣΗ: Η στήλη Sample δεν είναι απόλυτα ευθυγραμμισμένη. Γίνεται merge...")
    merged = pd.merge(meta, X, on='Sample', how='inner')
    meta = merged[meta.columns]
    X = merged[X.columns]

print(f"Αρχικό dataset: {len(X)} δείγματα.")

# --- Model A: CRC vs Control (All genders) ---
mask_A = meta['Condition'].isin(['CRC', 'Control'])
meta_A = meta[mask_A]
X_A = X[mask_A]

# --- Model B: CRC vs Control (Male only) ---
mask_B = mask_A & (meta['Gender'] == 'Male')
meta_B = meta[mask_B]
X_B = X[mask_B]

# --- Model C: CRC vs Control (Female only) ---
mask_C = mask_A & (meta['Gender'] == 'Female')
meta_C = meta[mask_C]
X_C = X[mask_C]

# --- Model D: Adenoma vs Control (All genders) ---
mask_D = meta['Condition'].isin(['Adenoma', 'Control'])
meta_D = meta[mask_D]
X_D = X[mask_D]

# --- Adenoma Scoring Set (for Model A later) ---
mask_Adenoma = meta['Condition'] == 'Adenoma'
meta_Adenoma = meta[mask_Adenoma]
X_Adenoma = X[mask_Adenoma]

# --- Save outputs ---
print("\nΑποθήκευση συνόλων δεδομένων...")

# Model A
X_A.to_csv(os.path.join(OUTPUT_DIR, "data_Model_A.csv"), index=False)
meta_A.to_csv(os.path.join(OUTPUT_DIR, "meta_Model_A.csv"), index=False)
print(f"Model A (Όλα CRC+Control): {len(X_A)} δείγματα")

# Model B
X_B.to_csv(os.path.join(OUTPUT_DIR, "data_Model_B.csv"), index=False)
meta_B.to_csv(os.path.join(OUTPUT_DIR, "meta_Model_B.csv"), index=False)
print(f"Model B (Άνδρες CRC+Control): {len(X_B)} δείγματα")

# Model C
X_C.to_csv(os.path.join(OUTPUT_DIR, "data_Model_C.csv"), index=False)
meta_C.to_csv(os.path.join(OUTPUT_DIR, "meta_Model_C.csv"), index=False)
print(f"Model C (Γυναίκες CRC+Control): {len(X_C)} δείγματα")

# Model D
X_D.to_csv(os.path.join(OUTPUT_DIR, "data_Model_D.csv"), index=False)
meta_D.to_csv(os.path.join(OUTPUT_DIR, "meta_Model_D.csv"), index=False)
print(f"Model D (Όλα Adenoma+Control): {len(X_D)} δείγματα")

# Adenoma
X_Adenoma.to_csv(os.path.join(OUTPUT_DIR, "data_Adenoma.csv"), index=False)
meta_Adenoma.to_csv(os.path.join(OUTPUT_DIR, "meta_Adenoma.csv"), index=False)
print(f"Adenoma (Για scoring): {len(X_Adenoma)} δείγματα")

print("\nΌλα τα αρχεία δημιουργήθηκαν επιτυχώς στο:", OUTPUT_DIR)
