import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

INPUT_DIR = "OUTPUTS/phase_4/4.3_advanced"
OUTPUT_DIR = os.path.join(INPUT_DIR, "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Set style
sns.set_theme(style="whitegrid", palette="muted")
plt.rcParams.update({'font.size': 12})

def load_data():
    adv = pd.read_csv(os.path.join(INPUT_DIR, "advanced_metrics.csv"))
    fair = pd.read_csv(os.path.join(INPUT_DIR, "fairness_metrics.csv"))
    ade = pd.read_csv(os.path.join(INPUT_DIR, "adenoma_predictions.csv"))
    return adv, fair, ade

def plot_confusion_matrix(adv):
    # Sum CM for Model A, RandomForest
    ma_rf = adv[(adv['Model'] == 'Model_A') & (adv['Algorithm'] == 'RandomForest')]
    if len(ma_rf) == 0: return
    
    tp = ma_rf['TP'].sum()
    fp = ma_rf['FP'].sum()
    tn = ma_rf['TN'].sum()
    fn = ma_rf['FN'].sum()
    
    cm = np.array([[tn, fp], [fn, tp]])
    
    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=['Control', 'CRC'], 
                yticklabels=['Control', 'CRC'])
    plt.title('Συνολικό Confusion Matrix\n(Model A - Random Forest)')
    plt.ylabel('Πραγματική Κατάσταση')
    plt.xlabel('Πρόβλεψη Μοντέλου')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'confusion_matrix_ModelA_RF.png'), dpi=300)
    plt.close()

def plot_fairness(fair):
    # Mean across folds for RandomForest
    rf_fair = fair[fair['Algorithm'] == 'RandomForest'].mean(numeric_only=True)
    
    # Sex-Stratified AUC
    plt.figure(figsize=(6, 5))
    sns.barplot(x=['Γυναίκες', 'Άνδρες'], y=[rf_fair['AUC_F'], rf_fair['AUC_M']], palette=['#ff9999', '#66b3ff'])
    plt.title('Sex-Stratified AUC (Model A - Random Forest)')
    plt.ylim(0.5, 0.9)
    plt.ylabel('AUC')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'fairness_auc_stratified.png'), dpi=300)
    plt.close()
    
    # Parity Metrics
    metrics = ['Demographic Parity\n(P(y=1))', 'Equal Opportunity\n(TPR/Recall)', 'Predictive Parity\n(PPV/Precision)']
    f_vals = [rf_fair['DP_F'], rf_fair['TPR_F'], rf_fair['PPV_F']]
    m_vals = [rf_fair['DP_M'], rf_fair['TPR_M'], rf_fair['PPV_M']]
    
    x = np.arange(len(metrics))
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(8, 6))
    rects1 = ax.bar(x - width/2, f_vals, width, label='Γυναίκες', color='#ff9999')
    rects2 = ax.bar(x + width/2, m_vals, width, label='Άνδρες', color='#66b3ff')
    
    ax.set_ylabel('Τιμή Μετρικής')
    ax.set_title('Fairness Metrics ανά Φύλο (Model A - RF)')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.legend()
    plt.ylim(0, 1)
    
    for container in ax.containers:
        ax.bar_label(container, fmt='%.3f', padding=3)
        
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'fairness_metrics_comparison.png'), dpi=300)
    plt.close()

def plot_adenoma_scores(ade):
    rf_ade = ade[ade['Algorithm'] == 'RandomForest']
    if len(rf_ade) == 0: return
    
    plt.figure(figsize=(8, 5))
    sns.histplot(rf_ade['P_CRC'], bins=20, kde=True, color='purple')
    plt.axvline(x=0.5, color='red', linestyle='--', label='Threshold (0.5)')
    
    # Calculate % CRC-like
    pct_crc = (rf_ade['P_CRC'] > 0.5).mean() * 100
    
    plt.title(f'Κατανομή Κινδύνου CRC στα Αδενώματα\n(Ποσοστό CRC-like: {pct_crc:.1f}%)')
    plt.xlabel('Πιθανότητα CRC (Predicted P(CRC))')
    plt.ylabel('Αριθμός Δειγμάτων')
    plt.xlim(0, 1)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'adenoma_risk_distribution.png'), dpi=300)
    plt.close()

if __name__ == "__main__":
    try:
        adv, fair, ade = load_data()
        plot_confusion_matrix(adv)
        plot_fairness(fair)
        plot_adenoma_scores(ade)
        print(f"Ολοκληρώθηκε! Τα διαγράμματα σώθηκαν στο {OUTPUT_DIR}")
    except Exception as e:
        print("Σφάλμα:", e)
