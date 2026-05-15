import pandas as pd
import numpy as np
import os

INPUT_DIR = "OUTPUTS/phase_4/4.1/4.1_results"
OUTPUT_DIR = "OUTPUTS/phase_4/4.2_aggregation"
os.makedirs(OUTPUT_DIR, exist_ok=True)

df = pd.read_csv(os.path.join(INPUT_DIR, "advanced_metrics.csv"))

# ─────────────────────────────────────────────────
# 4.2a  Per-Model Summary  (Median AUC ± IQR)
#       + Per-fold breakdown table
# ─────────────────────────────────────────────────

models = ['Model_A', 'Model_B', 'Model_C', 'Model_D']
algos  = ['RandomForest', 'LASSO']

# ── Summary table ──
rows = []
for model in models:
    for algo in algos:
        sub = df[(df['Model'] == model) & (df['Algorithm'] == algo)]
        if len(sub) == 0:
            continue
        for metric in ['AUC', 'AUPRC', 'F1']:
            vals = sub[metric]
            rows.append({
                'Model': model,
                'Algorithm': algo,
                'Metric': metric,
                'Median': vals.median(),
                'IQR_25': vals.quantile(0.25),
                'IQR_75': vals.quantile(0.75),
                'Mean': vals.mean(),
                'Std': vals.std(),
                'N_folds': len(vals)
            })

summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(OUTPUT_DIR, "4.2a_per_model_summary.csv"), index=False)

# ── Per-fold breakdown tables (one per model × algo) ──
for model in models:
    for algo in algos:
        sub = df[(df['Model'] == model) & (df['Algorithm'] == algo)].copy()
        if len(sub) == 0:
            continue
        sub = sub.sort_values('Study')
        cols = ['Study', 'AUC', 'AUPRC', 'F1', 'TP', 'FP', 'TN', 'FN']
        fname = f"4.2a_fold_breakdown_{model}_{algo}.csv"
        sub[cols].to_csv(os.path.join(OUTPUT_DIR, fname), index=False)

# ── Print summary ──
print("\n" + "="*70)
print("4.2a  Per-Model Summary  (Median ± IQR)")
print("="*70)
for model in models:
    print(f"\n── {model} ──")
    for algo in algos:
        s = summary[(summary['Model'] == model) & (summary['Algorithm'] == algo) & (summary['Metric'] == 'AUC')]
        if len(s) == 0:
            continue
        r = s.iloc[0]
        print(f"  {algo:15s}  AUC = {r['Median']:.3f}  [{r['IQR_25']:.3f} – {r['IQR_75']:.3f}]  (n={int(r['N_folds'])} folds)")

print(f"\nΑρχεία αποθηκεύτηκαν στο {OUTPUT_DIR}")
