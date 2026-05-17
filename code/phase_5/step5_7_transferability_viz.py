#!/usr/bin/env python3
"""
STEP 5.7 — TRANSFERABILITY VISUALIZATION (Pure Matplotlib)
===========================================================

Reads LODO fold breakdowns (RandomForest) for Model A (Universal), 
Model B (Male-Only), and Model C (Female-Only) and generates 
high-quality, publication-ready visualizations using ONLY matplotlib.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define paths
AGG_DIR = "OUTPUTS/phase_4/4.2_aggregation"
OUT_DIR = "OUTPUTS/phase_5/5.7_transferability"
os.makedirs(OUT_DIR, exist_ok=True)

# Load LODO RF results
def load_results():
    files = {
        'Universal (All)': '4.2a_fold_breakdown_Model_A_RandomForest.csv',
        'Male-Only': '4.2a_fold_breakdown_Model_B_RandomForest.csv',
        'Female-Only': '4.2a_fold_breakdown_Model_C_RandomForest.csv'
    }
    
    data = {}
    for label, fname in files.items():
        path = os.path.join(AGG_DIR, fname)
        if os.path.exists(path):
            df = pd.read_csv(path)
            data[label] = df.set_index('Study')['AUC'].to_dict()
        else:
            print(f"Warning: {path} not found.")
            
    if not data:
        raise FileNotFoundError("No fold breakdown CSV files found.")
        
    return data

def main():
    # Set premium aesthetic styles
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'figure.titlesize': 16
    })
    
    data = load_results()
    
    # Sort studies by Universal model performance for a clean visual gradient
    universal_auc = data['Universal (All)']
    study_order = sorted(universal_auc.keys(), key=lambda k: universal_auc[k], reverse=True)
    
    # 1. Grouped Bar Chart of AUC per Study
    fig, ax = plt.subplots(figsize=(14, 7), dpi=300)
    
    # Elegant color palette
    colors = {
        'Universal (All)': '#4A5568',  # Slate Gray
        'Male-Only': '#2B6CB0',       # Deep Blue
        'Female-Only': '#D53F8C'      # Vibrant Pink/Magenta
    }
    
    x = np.arange(len(study_order))
    width = 0.25
    
    rects1 = ax.bar(x - width, [data['Universal (All)'].get(s, 0) for s in study_order], width, label='Universal (All)', color=colors['Universal (All)'], edgecolor='black', linewidth=0.8, alpha=0.9)
    rects2 = ax.bar(x, [data['Male-Only'].get(s, 0) for s in study_order], width, label='Male-Only', color=colors['Male-Only'], edgecolor='black', linewidth=0.8, alpha=0.9)
    rects3 = ax.bar(x + width, [data['Female-Only'].get(s, 0) for s in study_order], width, label='Female-Only', color=colors['Female-Only'], edgecolor='black', linewidth=0.8, alpha=0.9)
    
    # Add a horizontal line at AUC = 0.50 (random baseline) and AUC = 0.70 (transferability threshold)
    ax.axhline(0.50, color='#E53E3E', linestyle='--', linewidth=1.2, label='Random Guessing (AUC = 0.50)')
    ax.axhline(0.70, color='#319795', linestyle=':', linewidth=1.5, label='Transferability Threshold (AUC = 0.70)')
    
    ax.set_title('Cross-Cohort Generalizability: LODO AUC per Study (Random Forest)', pad=15, fontweight='bold')
    ax.set_xlabel('Held-out Study (Unseen during training)', labelpad=12)
    ax.set_ylabel('Area Under ROC Curve (AUC)', labelpad=12)
    ax.set_ylim(0.40, 1.05)
    ax.set_xticks(x)
    ax.set_xticklabels(study_order, rotation=25, ha='right')
    ax.grid(axis='y', linestyle=':', alpha=0.6)
    
    # Add values on top of bars
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            if height > 0.05:
                ax.annotate(f'{height:.2f}',
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 3),  # 3 points vertical offset
                            textcoords="offset points",
                            ha='center', va='bottom', fontsize=8, color='#2D3748')
                
    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)
    
    ax.legend(frameon=True, facecolor='white', edgecolor='#E2E8F0', loc='lower left')
    plt.tight_layout()
    bar_path = os.path.join(OUT_DIR, 'lodo_auc_comparison_bar.png')
    plt.savefig(bar_path, bbox_inches='tight')
    plt.close()
    print(f"Saved grouped bar chart to {bar_path}")

    # 2. Boxplot of AUC Distribution (Dispersion / Stability)
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    
    cohorts = ['Universal (All)', 'Male-Only', 'Female-Only']
    box_data = [[data[c][s] for s in study_order if s in data[c]] for c in cohorts]
    
    # Draw boxplot
    bp = ax.boxplot(
        box_data,
        patch_artist=True,
        widths=0.4,
        boxprops=dict(linewidth=1.5, alpha=0.7),
        whiskerprops=dict(linewidth=1.5, color='#4A5568'),
        capprops=dict(linewidth=1.5, color='#4A5568'),
        medianprops=dict(linewidth=2.0, color='#1A202C'),
        flierprops=dict(marker='o', markerfacecolor='#E53E3E', markersize=6)
    )
    
    # Color boxplots
    for patch, color in zip(bp['boxes'], [colors[c] for c in cohorts]):
        patch.set_facecolor(color)
        
    # Overlay individual points for completeness
    for i, c in enumerate(cohorts):
        vals = box_data[i]
        x_jitter = np.random.normal(i + 1, 0.04, size=len(vals))
        ax.scatter(x_jitter, vals, color='#1A202C', s=35, alpha=0.8, zorder=3)
        
        # Annotate median value
        median_val = np.median(vals)
        ax.text(i + 1, median_val + 0.02, f'Median: {median_val:.3f}', 
                ha='center', va='bottom', fontweight='bold', color='#1A202C', fontsize=10)
        
    ax.axhline(0.70, color='#319795', linestyle=':', linewidth=1.2)
    ax.set_title('Robustness Comparison: LODO AUC Dispersion across 10 Cohorts', pad=15, fontweight='bold')
    ax.set_xlabel('Cohort Model', labelpad=12)
    ax.set_ylabel('LODO AUC Distribution', labelpad=12)
    ax.set_xticklabels(cohorts)
    ax.set_ylim(0.45, 1.05)
    ax.grid(axis='y', linestyle=':', alpha=0.6)
    
    plt.tight_layout()
    box_path = os.path.join(OUT_DIR, 'lodo_auc_dispersion_box.png')
    plt.savefig(box_path, bbox_inches='tight')
    plt.close()
    print(f"Saved boxplot chart to {box_path}")

if __name__ == '__main__':
    main()
