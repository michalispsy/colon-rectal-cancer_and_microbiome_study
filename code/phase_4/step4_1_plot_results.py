import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

INPUT_DIR = "OUTPUTS/phase_4/4.1/4.1_training"
OUTPUT_DIR = os.path.join(INPUT_DIR, "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Φόρτωση δεδομένων
try:
    df_summary = pd.read_csv(os.path.join(INPUT_DIR, "models_summary.csv"))
    df_cv = pd.read_csv(os.path.join(INPUT_DIR, "all_models_cv_results.csv"))
except FileNotFoundError as e:
    print(f"Σφάλμα: Δεν βρέθηκαν τα αρχεία. {e}")
    exit(1)

# Set style
sns.set_theme(style="whitegrid", palette="muted")
plt.rcParams.update({'font.size': 12})

# 1. Barplot για το AUC
plt.figure(figsize=(10, 6))
ax = sns.barplot(x='Model', y='AUC', hue='Algorithm', data=df_summary)
plt.title('Μέσο AUC ανά Μοντέλο και Αλγόριθμο (LODO Cross-Validation)', pad=20, fontsize=14, fontweight='bold')
plt.ylim(0.4, 0.9)
plt.axhline(y=0.5, color='r', linestyle='--', alpha=0.5, label='Random Guess (0.5)')
plt.ylabel('Area Under ROC Curve (AUC)')

# Προσθήκη ετικετών (νούμερα) πάνω από τις μπάρες
for container in ax.containers:
    ax.bar_label(container, fmt='%.3f', padding=3)

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'auc_comparison_bar.png'), dpi=300, bbox_inches='tight')
plt.close()

# 2. Boxplot για τη διασπορά του AUC στα Folds (από το df_cv)
plt.figure(figsize=(12, 6))
sns.boxplot(x='Model', y='AUC', hue='Algorithm', data=df_cv)
# Προσθήκη των πραγματικών σημείων (folds) πάνω από τα boxes
sns.stripplot(x='Model', y='AUC', hue='Algorithm', data=df_cv, 
              dodge=True, color='black', alpha=0.5, size=4, legend=False)
plt.title('Διασπορά του AUC ανά LODO Fold', pad=20, fontsize=14, fontweight='bold')
plt.axhline(y=0.5, color='r', linestyle='--', alpha=0.5)
plt.ylabel('AUC σε κάθε ξεχωριστό Test Fold')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'auc_dispersion_box.png'), dpi=300, bbox_inches='tight')
plt.close()

# 3. Barplot για F1 Score
plt.figure(figsize=(10, 6))
ax = sns.barplot(x='Model', y='F1', hue='Algorithm', data=df_summary)
plt.title('Μέσο F1-Score ανά Μοντέλο και Αλγόριθμο', pad=20, fontsize=14, fontweight='bold')
plt.ylim(0, 0.9)
plt.ylabel('F1 Score')

for container in ax.containers:
    ax.bar_label(container, fmt='%.3f', padding=3)

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'f1_comparison_bar.png'), dpi=300, bbox_inches='tight')
plt.close()

print(f"Τα διαγράμματα αποθηκεύτηκαν επιτυχώς στο {OUTPUT_DIR}")
