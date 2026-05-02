import os
import pandas as pd
import shutil

ROOT_DIR = "/home/michalis/Documents/ece_ntua/8th/health"
DATA_DIR = os.path.join(ROOT_DIR, "data/crc_study_final_data/species_level")

METADATA_PATH = os.path.join(DATA_DIR, "metadata.csv")
X_EXPLORE_PATH = os.path.join(DATA_DIR, "X_explore.csv")
SPECIES_FILTERED_PATH = os.path.join(DATA_DIR, "species_filtered_crossstudy.csv")

# Backup the original files just in case
print("Backing up original files...")
for f in [METADATA_PATH, X_EXPLORE_PATH, SPECIES_FILTERED_PATH]:
    if os.path.exists(f) and not os.path.exists(f + ".bak"):
        shutil.copy(f, f + ".bak")

print("Loading metadata...")
meta = pd.read_csv(METADATA_PATH)
initial_samples = len(meta)

# Find samples from GuptaA_2019
gupta_samples = meta[meta['Study'] == 'GuptaA_2019']['Sample'].tolist()
print(f"Found {len(gupta_samples)} samples from GuptaA_2019.")

# Filter metadata
meta_filtered = meta[meta['Study'] != 'GuptaA_2019']
print(f"Metadata samples: {initial_samples} -> {len(meta_filtered)}")
meta_filtered.to_csv(METADATA_PATH, index=False)

# Filter X_explore
if os.path.exists(X_EXPLORE_PATH):
    print("Loading X_explore.csv...")
    x_explore = pd.read_csv(X_EXPLORE_PATH)
    x_filtered = x_explore[~x_explore['Sample'].isin(gupta_samples)]
    print(f"X_explore samples: {len(x_explore)} -> {len(x_filtered)}")
    x_filtered.to_csv(X_EXPLORE_PATH, index=False)

# Filter species_filtered_crossstudy
if os.path.exists(SPECIES_FILTERED_PATH):
    print("Loading species_filtered_crossstudy.csv...")
    species = pd.read_csv(SPECIES_FILTERED_PATH)
    species_filtered = species[~species['Sample'].isin(gupta_samples)]
    print(f"species_filtered samples: {len(species)} -> {len(species_filtered)}")
    species_filtered.to_csv(SPECIES_FILTERED_PATH, index=False)

print("Successfully removed GuptaA_2019 from the dataset.")
