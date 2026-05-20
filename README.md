# 🦠 CRC Microbiome Study — Sex-Stratified Machine Learning for Colorectal Cancer Detection

<div align="center">

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-4.x-276DC3?logo=r&logoColor=white)
![scikit-learn](https://img.shields.io/badge/scikit--learn-ML-F7931E?logo=scikit-learn&logoColor=white)
![SHAP](https://img.shields.io/badge/SHAP-Explainability-ff6b6b)
![License](https://img.shields.io/badge/License-Academic-green)
![Status](https://img.shields.io/badge/Status-Complete-brightgreen)

**A multi-cohort metagenomics study investigating sex-specific gut microbiome dysbiosis in Colorectal Cancer using LODO cross-validation, Hurdle encoding, and SHAP-based interpretability.**

*NTUA — School of Electrical & Computer Engineering | Course: Mobile & Electronic Health Technologies*

</div>

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Key Results](#-key-results)
- [Dataset](#-dataset)
- [Methodology](#-methodology)
- [Repository Structure](#-repository-structure)
- [Pipeline Phases](#-pipeline-phases)
- [Setup & Reproducibility](#-setup--reproducibility)
- [Key Findings](#-key-findings)
- [Authors](#-authors)

---

## 🔬 Overview

This project presents a rigorous, multi-cohort machine learning pipeline to detect **Colorectal Cancer (CRC)** from fecal metagenomics data, with a primary focus on uncovering **sex-specific microbial signatures**. The analysis spans **10 independent international cohorts** (1,544 samples, 8 countries, 3 continents) and employs **Leave-One-Dataset-Out (LODO) cross-validation** — the gold standard for evaluating cross-study generalizability in microbiome research.

**Core Research Questions:**
1. Can gut microbiome data non-invasively detect CRC with clinically useful accuracy?
2. Do men and women exhibit fundamentally different microbial dysbiosis pathways toward CRC?
3. Does a universal model treat both sexes fairly, or is sex-specific screening warranted?
4. At what stage of disease progression (adenoma → carcinoma) do microbial signatures emerge?

---

## 🏆 Key Results

| Model | Target | Median LODO AUC | IQR |
|:---|:---|:---:|:---|
| **A — Universal** | CRC vs Control (all) | **0.748** | [0.691 – 0.824] |
| **B — Male-Only** | CRC vs Control (♂) | **0.774** | [0.722 – 0.833] |
| **C — Female-Only** | CRC vs Control (♀) | **0.686** | [0.614 – 0.771] |
| **D — Early Detection** | Adenoma vs Control | TBD | — |
| **E — Stage Discrimination** | CRC vs Adenoma | TBD | — |
| **F — Pan-Neoplasia** | (CRC+Adenoma) vs Control | TBD | — |

**Fairness (Model A):** AUC gap = 0.008 · Equal Opportunity gap = 9.5% · Demographic Parity gap = 9.7% → **Clinically fair model**

**Adenoma Progression:** 42.1% of adenoma samples scored as "CRC-like" (P > 0.5), confirming measurable microbial dysbiosis onset at the pre-cancerous stage.

---

## 📊 Dataset

| Property | Value |
|:---|:---|
| **Source** | `curatedMetagenomicData` (Bioconductor) |
| **Processing pipeline** | MetaPhlAn3 + HUMAnN3 (runDate: 2021-03-31) |
| **Taxonomic level** | Species-level profiles |
| **Total samples** | 1,544 |
| **Cohorts** | 10 (post GuptaA_2019 exclusion) |
| **Countries / Continents** | 8 countries · Europe, Asia, North America |
| **CRC / Control / Adenoma** | 671 / 664 / 209 |
| **Male / Female** | 929 (60.2%) / 615 (39.8%) |
| **Species (post-filtering)** | 107 (from 934 original; >10% prevalence threshold) |
| **Feature space (Hurdle)** | 214 (107 PA + 107 rCLR) + sex + age = **216** |
| **Matrix sparsity** | 46.66% |

<details>
<summary><b>Cohort breakdown</b></summary>

| Study | n | Country | Continent | Has Adenoma |
|:---|:---:|:---|:---|:---:|
| YachidaS_2019 | 576 | 🇯🇵 Japan | Asia | ✅ |
| ZellerG_2014 | 156 | 🇫🇷 France | Europe | ✅ |
| FengQ_2015 | 154 | 🇦🇹 Austria | Europe | ✅ |
| YuJ_2015 | 128 | 🇨🇳 China | Asia | ❌ |
| WirbelJ_2018 | 125 | 🇩🇪 Germany | Europe | ❌ |
| VogtmannE_2016 | 104 | 🇺🇸 USA | N. America | ❌ |
| HanniganGD_2017 | 81 | 🇺🇸 USA / 🇨🇦 Canada | N. America | ✅ |
| ThomasAM_2018a | 80 | 🇮🇹 Italy | Europe | ✅ |
| ThomasAM_2019_c | 80 | 🇮🇹 Italy | Europe | ❌ |
| ThomasAM_2018b | 60 | 🇮🇹 Italy | Europe | ❌ |

</details>

---

## 🧬 Methodology

### Hurdle Encoding Strategy

To handle the zero-inflated, compositional nature of metagenomic data, each species generates **two features** per sample — with **no pseudocount**:

```
PA_i   = 1  if x_i > 0, else 0          (presence/absence)

rCLR_i = log(x_i / geomean_nonzero(sample))  if x_i > 0
       = 0                                    if x_i = 0
```

This disentangles the **detection signal** (PA) from the **relative abundance signal** (rCLR), allowing the model to separately weight "is this bacterium present?" vs "how abundant is it?".

### LODO Cross-Validation

```
For each of the 10 cohorts as held-out test set:
  1. Train on 9 remaining cohorts
  2. Batch-correct rCLR features on train only (MMUPHin R / ComBat Python)
     → avoids data leakage into test set
  3. Frozen Z-score normalization (fit on train → transform test)
  4. Train: Random Forest + LASSO Logistic Regression
     (inner 5-fold CV for hyperparameter tuning)
  5. Evaluate: AUC, AUPRC, F1, sex-stratified AUC, fairness metrics
```

### Statistical Validation (Phase 3)

Before ML, three complementary statistical tests establish biologically valid biomarkers:
- **Blocked Wilcoxon** (rCLR features) — controls for study as a blocking factor
- **Cochran-Mantel-Haenszel** (PA features) — stratified odds ratios across cohorts
- **Random Effects Meta-Analysis** — pooled effect sizes with I² heterogeneity

These serve as ground truth for SHAP cross-validation in Phase 5.

---

## 📁 Repository Structure

```
health/
├── code/
│   ├── phase_1/            # Data preparation & Hurdle encoding
│   │   └── regen_phase1_figures.py
│   ├── phase_2/            # Exploratory Data Analysis
│   │   ├── step2_1_dimensionality_reduction.py   # PCA, t-SNE, UMAP
│   │   ├── step2_2_1_permanova.py                # Global + per-study PERMANOVA
│   │   ├── step2_2_2_remove_gupta.py             # GuptaA_2019 exclusion
│   │   ├── step2_3_centroid_distances.py         # PCA centroid analysis
│   │   └── step2_4_unsupervised_clustering.py    # K-Means, DBSCAN, ARI
│   ├── phase_3/            # Differential Abundance Analysis
│   │   ├── step3_1_blocked_wilcoxon_rclr.py      # rCLR Wilcoxon tests
│   │   ├── step3_2_cmh_pa.py                     # CMH PA tests
│   │   ├── step3_3_random_effects_meta_analysis.py
│   │   ├── step3_4_visualization.py              # Volcano, forest, heatmap
│   │   └── step3_5_adenoma_carcinoma_sequence.py # EARLY/LATE/PROGRESSION timeline
│   ├── phase_4/            # ML Pipeline — LODO + Fairness
│   │   ├── step4_0_prepare_data.py
│   │   ├── step4_1_lodo_prep.py
│   │   ├── step4_1_training.py
│   │   ├── step4_2a_aggregation.py               # Per-model summary
│   │   ├── step4_2b_cross_model.py               # Cross-model comparison
│   │   ├── step4_2c_fairness.py                  # Sex fairness metrics
│   │   ├── step4_2de_aggregation.py              # Adenoma scoring
│   │   ├── step4_3_advanced_training.py
│   │   └── mmuphin_wrapper.R                     # Batch correction (R)
│   └── phase_5/            # Interpretability & Validation
│       ├── step5_1_shap_analysis.py              # SHAP per model (PA vs rCLR)
│       ├── step5_2_cross_validate_shap_phase3.py # SHAP × Phase 3 overlap
│       ├── step5_2b_lasso_interpretability.py
│       ├── step5_3_shap_beeswarm.py              # Beeswarm plots (sex split)
│       ├── step5_3b_shap_beeswarm_all_rf.py
│       ├── step5_4_local_xai.py                  # SHAP waterfall (per patient)
│       ├── step5_5_subgroup_analysis.py           # Age, country, study breakdown
│       ├── step5_6_biomarker_shifts.py            # Bump charts (rank shifts)
│       └── step5_7_transferability_viz.py         # LODO generalizability viz
├── data/
│   └── crc_study_final_data/    # Processed species-level profiles (tracked)
│       └── species_level/
├── OUTPUTS/
│   ├── phase_1/ … phase_5/      # All figures, CSV results, reports
└── pipeline.txt                  # Full methodological specification
```

---

## 🔄 Pipeline Phases

```
Phase 1: Data Preparation
  └─ MetaPhlAn3 profiles → Prevalence filtering (107 species)
     → Hurdle encoding (PA + rCLR) → 214-feature matrix (1,544 × 214)

Phase 2: Exploratory Analysis  ──┐
  └─ PCA/t-SNE/UMAP              │   (parallel, independent of ML)
     PERMANOVA (R² per factor)    │
     Centroid distances           │
     Unsupervised clustering      │

Phase 3: Differential Abundance ─┘
  └─ Blocked Wilcoxon + CMH + Meta-Analysis
     → 9 cross-study reproducible biomarkers

Phase 4: ML Pipeline (LODO)
  └─ 6 models (Universal / Male / Female / Adenoma / Stage / Pan-Neoplasia)
     → Random Forest + LASSO, inner CV tuning
     → AUC, AUPRC, Fairness metrics, Adenoma scoring

Phase 5: Interpretability & Validation
  └─ SHAP TreeExplainer → Consensus biomarker rankings
     Cross-validation vs Phase 3 statistics
     Bump charts (rank shifts across sex splits)
     Local XAI (waterfall per patient: TP/FN/TN/FP)
     Subgroup analysis, Transferability analysis
```

---

## ⚙️ Setup & Reproducibility

### Requirements

**Python (≥ 3.10)**
```bash
pip install numpy pandas scikit-learn shap matplotlib seaborn scipy statsmodels
pip install umap-learn pycombat pingouin
```

**R (≥ 4.x)** — for batch correction only
```r
install.packages("MMUPHin")
```

### Data

The processed `crc_study_final_data/` directory is tracked in git (raw downloads are gitignored). Data originates from the [`curatedMetagenomicData`](https://bioconductor.org/packages/curatedMetagenomicData/) Bioconductor package, MetaPhlAn3 pipeline, runDate `2021-03-31`.

### Running the Pipeline

```bash
# Phase 1 — Regenerate figures
python code/phase_1/regen_phase1_figures.py

# Phase 2 — EDA
python code/phase_2/step2_1_dimensionality_reduction.py
python code/phase_2/step2_2_1_permanova.py
python code/phase_2/step2_3_centroid_distances.py
python code/phase_2/step2_4_unsupervised_clustering.py

# Phase 3 — Differential Abundance
python code/phase_3/step3_1_blocked_wilcoxon_rclr.py
python code/phase_3/step3_2_cmh_pa.py
python code/phase_3/step3_3_random_effects_meta_analysis.py
python code/phase_3/step3_5_adenoma_carcinoma_sequence.py

# Phase 4 — ML (computationally intensive)
python code/phase_4/step4_0_prepare_data.py
python code/phase_4/step4_1_training.py      # runs LODO for all 6 models
python code/phase_4/step4_2a_aggregation.py
python code/phase_4/step4_2c_fairness.py

# Phase 5 — Interpretability
python code/phase_5/step5_1_shap_analysis.py
python code/phase_5/step5_2_cross_validate_shap_phase3.py
python code/phase_5/step5_3b_shap_beeswarm_all_rf.py
python code/phase_5/step5_4_local_xai.py
python code/phase_5/step5_7_transferability_viz.py
```

---

## 🔍 Key Findings

### 1. Two Fundamentally Different Dysbiosis Pathways

SHAP beeswarm analysis (Models B & C) reveals that men and women do not simply differ in *degree* of CRC-associated dysbiosis — they differ in *mechanism*:

| | Men (Model B) | Women (Model C) |
|:---|:---|:---|
| **Primary signal** | Collapse of protective microbiome | Dysregulation of *Bacteroides* |
| **Key protective species** | *Roseburia intestinalis* (SHAP −0.07), *R. faecis*, *Lachnospira eligens*, *F. prausnitzii* | *Dorea formicigenerans* (#3, SHAP −0.05) |
| **Key risk species** | *Akkermansia muciniphila* (SHAP +0.08) | *B. fragilis* ETBF strains (SHAP +0.05), *B. thetaiotaomicron* |
| **Universal #1 biomarker** | *Parvimonas micra* (SHAP +0.13) | *Parvimonas micra* (SHAP +0.08, weaker) |

> The male signal is characterised by Lachnospiraceae loss; the female signal by Bacteroides gain — pointing to different molecular oncogenesis pathways.

### 2. Biomarker Stability (Bump Charts)

Only ***Parvimonas micra*** and ***Ruthenibacterium lactatiformans*** maintain stable top rankings across all three model variants (All/Male/Female), confirming them as **universal pan-sex CRC biomarkers**.

Sex-specific biomarkers show dramatic rank shifts: *Roseburia intestinalis* drops from #2 (men) to #17 (women); *Akkermansia muciniphila* from #7 to #37.

### 3. Fairness & Clinical Recommendation

The Universal model (Model A) achieves excellent demographic fairness (AUC gap: 0.008). Counterintuitively, the Female-only model (Model C, AUC=0.686) *underperforms* the Universal model on female patients (AUC=0.749), because it lacks the benefit of the stronger male signal. **Clinical recommendation: deploy the Universal model**, not sex-specific models.

### 4. Local XAI — Interpretable Errors

SHAP waterfall analysis on 4 archetype patients from WirbelJ_2018 (held-out):

| Case | P(CRC) | Explanation |
|:---|:---:|:---|
| True Positive | 0.78 | *P. micra* PA + rCLR together contribute +0.11 |
| False Negative | 0.30 | Absent *P. micra*; rich *Roseburia* — likely early-stage CRC |
| True Negative | 0.23 | All 15 top features push toward Control |
| False Positive | 0.65 | No dominant driver; diffuse microbial imbalance — possibly subclinical gut inflammation |

### 5. Transferability Constraints

Cross-cohort generalizability is bounded by three domain shifts:
- **Clinical shift**: Stage III/IV cohorts (ThomasAM_2019_c: AUC 0.996–1.000) far outperform early-detection screening cohorts (YachidaS_2019: AUC ~0.72)
- **Technical shift**: FOBT card preservation (VogtmannE_2016) and virome-optimised protocols (HanniganGD_2017: AUC ≈ 0.50) degrade anaerobic bacterial signals
- **Biogeographical shift**: Japanese dietary baseline (*Bacteroides plebeius* / seaweed-processing enzymes) introduces confounders invisible to Western-trained models

---

## 👥 Authors

- Παπακώστας Αχιλλέας
- Τσάκαλος Θεόδωρος
- Τσουρούφλη Ζωή
- Ψυχής Μιχάλης

**Supervisors:** Prof. K. Nikita, Prof. G. Matsopoulos, Prof. P. Tsanakas, Prof. I. Delis  
**Academic Advisor:** Theodόrou Glykeria  
**Institution:** National Technical University of Athens (NTUA)  
**Course:** Mobile & Electronic Health Technologies  
**Submission:** May 2026

---

<div align="center">

*For questions or collaboration, feel free to open an issue or reach out.*

</div>
