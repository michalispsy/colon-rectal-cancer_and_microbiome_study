import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

OUTPUT_DIR = "OUTPUTS/phase_5/5.6_biomarker_shifts"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define paths to consensus files
paths = {
    "All": "OUTPUTS/phase_5/5.1_shap/All/All_RandomForest_species_level_consensus.csv",
    "Male": "OUTPUTS/phase_5/5.1_shap/OnlyMale/OnlyMale_RandomForest_species_level_consensus.csv",
    "Female": "OUTPUTS/phase_5/5.1_shap/OnlyFemale/OnlyFemale_RandomForest_species_level_consensus.csv"
}

dfs = []
for model, path in paths.items():
    if os.path.exists(path):
        df = pd.read_csv(path)
        rank_col = 'SpeciesRank' if 'SpeciesRank' in df.columns else 'ConsensusRank'
        df_sub = df[['Species', rank_col]].copy()
        df_sub.columns = ['Species', 'Rank']
        df_sub['Model'] = model
        dfs.append(df_sub)

if dfs:
    all_ranks = pd.concat(dfs, ignore_index=True)
    
    def create_bump_chart(starting_model, title_suffix):
        # Determine column order based on starting model
        if starting_model == "All":
            cols = ["All", "Male", "Female"]
        elif starting_model == "Male":
            cols = ["Male", "All", "Female"]
        elif starting_model == "Female":
            cols = ["Female", "All", "Male"]
            
        top_microbes = all_ranks[(all_ranks['Model'] == starting_model) & (all_ranks['Rank'] <= 15)]['Species'].tolist()
        
        track_df = all_ranks[all_ranks['Species'].isin(top_microbes)].copy()
        pivot_df = track_df.pivot(index='Species', columns='Model', values='Rank')
        
        # Ensure all columns exist, fill missing with 50
        for c in cols:
            if c not in pivot_df.columns:
                pivot_df[c] = 50
                
        pivot_df = pivot_df[cols].fillna(50)
        
        # Sort index so colors match ranking order at the starting column
        # Rank values for starting model
        start_ranks = pivot_df[starting_model].sort_values()
        pivot_df = pivot_df.reindex(start_ranks.index)
        
        fig, ax = plt.subplots(figsize=(12, 8))
        cmap = plt.get_cmap('tab20')
        colors = [cmap(i) for i in np.linspace(0, 1, len(top_microbes))]
        
        for i, species in enumerate(pivot_df.index):
            ranks = pivot_df.loc[species].values
            plt.plot(cols, ranks, marker='o', linewidth=3, markersize=10, 
                     markeredgewidth=1.5, markeredgecolor='white', label=species, color=colors[i], alpha=0.85)

        plt.gca().invert_yaxis()
        plt.xticks(fontsize=12, fontweight='bold')
        plt.yticks([1, 5, 10, 15, 20, 30, 40, 50], labels=['1', '5', '10', '15', '20', '30', '40', 'Not in Top 50'], fontsize=10)
        plt.ylabel("Consensus Rank (Based on SHAP Stability & Impact)", fontsize=12, fontweight='bold')
        plt.title(f"Shift in Microbiome Biomarker Importance\n(Tracking Top 15 {title_suffix} Microbes)", fontsize=15, fontweight='bold', pad=20)
        
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0., frameon=False, fontsize=10, title=f"Top 15 {title_suffix}", title_fontsize=11)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.5)

        plt.tight_layout()
        filepath = f"{OUTPUT_DIR}/biomarker_rank_shifts_starting_{starting_model}.png"
        plt.savefig(filepath, dpi=300)
        plt.close()
        print(f"Bump chart saved to {filepath}")

    # Generate all three plots
    create_bump_chart("All", "Universal")
    create_bump_chart("Male", "Male-Specific")
    create_bump_chart("Female", "Female-Specific")
else:
    print("Consensus files not found.")
