import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

# Configure paths
ROOT_DIR = "/home/michalis/Documents/ece_ntua/8th/health"
X_EXPLORE_PATH = os.path.join(ROOT_DIR, "data/crc_study_final_data/species_level/X_explore.csv")
METADATA_PATH = os.path.join(ROOT_DIR, "data/crc_study_final_data/species_level/metadata.csv")
OUT_DIR = os.path.join(ROOT_DIR, "OUTPUTS/phase_2/2.1")
os.makedirs(OUT_DIR, exist_ok=True)

print("Loading data...")
X_df = pd.read_csv(X_EXPLORE_PATH)
meta_df = pd.read_csv(METADATA_PATH)

# Merge
df = pd.merge(X_df, meta_df, on="Sample", how="inner")
print(f"Merged dataset shape: {df.shape}")

# Features
feature_cols = [c for c in X_df.columns if c.startswith("PA_") or c.startswith("rCLR_")]
print(f"Number of features: {len(feature_cols)}")

X = df[feature_cols].values
X_scaled = StandardScaler().fit_transform(X)

print("Running PCA...")
pca = PCA(n_components=2, random_state=42)
pca_res = pca.fit_transform(X_scaled)
df['PCA1'] = pca_res[:, 0]
df['PCA2'] = pca_res[:, 1]
print(f"PCA Explained Variance Ratio: {pca.explained_variance_ratio_}")

print("Running t-SNE...")
tsne = TSNE(n_components=2, random_state=42, perplexity=30)
tsne_res = tsne.fit_transform(X_scaled)
df['tSNE1'] = tsne_res[:, 0]
df['tSNE2'] = tsne_res[:, 1]

# Visualization function
def plot_scatter(data, x_col, y_col, hue_col, title, filename):
    plt.figure(figsize=(10, 8))
    
    # Get unique groups
    # Fill NaN to avoid errors
    data_plot = data.copy()
    data_plot[hue_col] = data_plot[hue_col].fillna("Unknown")
    groups = data_plot[hue_col].unique()
    
    # Select colormap based on number of groups
    if len(groups) <= 10:
        cmap = plt.get_cmap("tab10")
    else:
        cmap = plt.get_cmap("tab20")
        
    for i, group in enumerate(groups):
        subset = data_plot[data_plot[hue_col] == group]
        plt.scatter(
            subset[x_col], subset[y_col],
            label=group,
            alpha=0.7,
            s=30,
            linewidth=0,
            color=cmap(i % cmap.N)
        )
        
    plt.title(title, fontsize=14, pad=15)
    plt.legend(title=hue_col, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, filename), dpi=300, bbox_inches='tight')
    plt.close()

print("Generating PCA plots...")
plot_scatter(df, "PCA1", "PCA2", "Study", "PCA - Colored by Study (Batch Effect)", "pca_study.png")
plot_scatter(df, "PCA1", "PCA2", "Condition", "PCA - Colored by Condition (CRC Signal)", "pca_condition.png")
plot_scatter(df, "PCA1", "PCA2", "Gender", "PCA - Colored by Gender (Sex Variation)", "pca_gender.png")

print("Generating t-SNE plots...")
plot_scatter(df, "tSNE1", "tSNE2", "Study", "t-SNE - Colored by Study (Batch Effect)", "tsne_study.png")
plot_scatter(df, "tSNE1", "tSNE2", "Condition", "t-SNE - Colored by Condition (CRC Signal)", "tsne_condition.png")
plot_scatter(df, "tSNE1", "tSNE2", "Gender", "t-SNE - Colored by Gender (Sex Variation)", "tsne_gender.png")

# Save coordinates for later if needed
coords_df = df[['Sample', 'Study', 'Condition', 'Gender', 'PCA1', 'PCA2', 'tSNE1', 'tSNE2']]
coords_df.to_csv(os.path.join(OUT_DIR, "embeddings_2d.csv"), index=False)
print("Saved embeddings to embeddings_2d.csv")
print("Step 2.1 complete.")
