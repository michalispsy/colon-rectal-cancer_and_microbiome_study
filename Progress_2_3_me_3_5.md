# Progress

## Dataset used

- Input matrix: `X_explore.csv`
- Metadata: `metadata.csv`
- Samples used: **1544**
- Features used: **214** total = **107 PA features + 107 rCLR features**
- Conditions:
  - Control: **664**
  - Adenoma: **209**
  - CRC: **671**
- Studies: **10**
- Sex distribution:
  - Male: **929**
  - Female: **615**

---

# Phase 2 — Exploratory analysis

## Step 2.3 — PCA centroid distances

**What we did:**  
Projected the hurdle-encoded data into PCA space and computed Euclidean distances between condition centroids: Control, Adenoma, CRC.

**Main result:**  
In the PCA space explaining 80% variance:

| Distance             |     Value |
| -------------------- | --------: |
| d(Adenoma, CRC)      | **2.136** |
| d(Adenoma, Control)  | **1.687** |
| d(CRC, Control)      | **1.612** |
| Adenoma CRC-likeness | **0.441** |

**Interpretation:**  
Adenoma was closer to **Control** than to **CRC** in global PCA centroid space. However, this is exploratory and Adenoma exists only in 5 studies, so we should not merge Adenoma with either group based only on this.

image:
	OUTPUTS/phase_2/step2_3_pca_centroids_PC12.png

---

## Step 2.4 — Unsupervised clustering

**What we did:**  
Ran K-Means with `k=2`, `k=3`, and `k=number_of_studies`, plus DBSCAN. Evaluated clustering using silhouette and ARI against study, condition, and sex.

**Main result:**

| Best agreement target | Best model | ARI |
|---|---|---:|
| Study | KMeans k=10 | **0.1445** |
| Condition | KMeans k=10 | **0.0058** |
| Sex | KMeans k=10 | **0.0019** |

**Interpretation:**  
Unsupervised clusters align much more with **study/batch** than with disease condition or sex. This supports using LODO later.

**Adenoma in k=2:**

| k=2 cluster label | Adenoma samples | Fraction |
|---|---:|---:|
| Weakly CRC-like | 125 / 209 | **59.8%** |
| Weakly Control-like | 84 / 209 | **40.2%** |

This is weak evidence only, because the CRC/control composition of the two clusters was almost balanced.

image:
	OUTPUTS/phase_2/step2_4_adenoma_k2_assignment.png

---

# Phase 3 — Differential abundance analysis

## Step 3.1 — rCLR features, blocked Wilcoxon / van Elteren test

**What we did:**  
Tested rCLR abundance differences using a study-blocked Wilcoxon-style test, equivalent to a van Elteren stratified rank test. FDR correction was applied within each comparison.

**Significant species counts, FDR < 0.05:**

| Comparison | Significant rCLR species |
|---|---:|
| CRC vs Control, all | **15** |
| CRC vs Control, male | **9** |
| CRC vs Control, female | **0** |
| Adenoma vs Control | **0** |
| CRC vs Adenoma | **1** |
| Male vs Female within Control | **5** |
| Male vs Female within Adenoma | **0** |
| Male vs Female within CRC | **4** |

**Important signal:**  
The strongest rCLR signal is **CRC vs Control**, especially in males. Adenoma vs Control was not significant after FDR in the blocked test.

Top CRC vs Control rCLR species included: `[Ruminococcus] torques`, `Parvimonas micra`, `Collinsella aerofaciens`, `Ruthenibacterium lactatiformans`, `Faecalibacterium prausnitzii`, `Roseburia faecis`.

image:
	OUTPUTS/phase_3/step3_1_blocked_wilcoxon_summary.png

---

## Step 3.2 — PA features, Cochran-Mantel-Haenszel test

**What we did:**  
Tested presence/absence differences using Cochran-Mantel-Haenszel tests, stratified by study. FDR correction was applied within each comparison.

**Significant species counts, FDR < 0.05:**

| Comparison                    | Significant PA species |
| ----------------------------- | ---------------------: |
| CRC vs Control, all           |                 **19** |
| CRC vs Control, male          |                 **18** |
| CRC vs Control, female        |                  **1** |
| Adenoma vs Control            |                  **2** |
| CRC vs Adenoma                |                  **1** |
| Male vs Female within Control |                  **5** |
| Male vs Female within Adenoma |                  **0** |
| Male vs Female within CRC     |                 **13** |

**Important signal:**  
Presence/absence also shows the strongest disease signal for **CRC vs Control**, again stronger in males than females.

Top PA species included `Parvimonas micra`, which was much more prevalent in CRC than Control.

image:
	OUTPUTS/phase_3/step3_2_cmh_summary.png

---

## Step 3.3 — Random-effects meta-analysis

**What we did:**  
Computed per-study standardized mean differences for each rCLR species and combined them with REML random-effects meta-analysis. We also computed I² and Cochran's Q.

**Significant species counts, FDR < 0.05:**

| Comparison | Stratum | Significant species |
|---|---|---:|
| CRC vs Control | All | **9** |
| CRC vs Control | Male | **7** |
| CRC vs Control | Female | **0** |
| Adenoma vs Control | All / Male / Female | **0** |
| CRC vs Adenoma | All / Male / Female | **0** |

**Top significant CRC vs Control species:**

| Species | Direction |
|---|---|
| `Ruthenibacterium lactatiformans` | Higher in CRC |
| `Lachnospira pectinoschiza` | Higher in Control |
| `Anaerostipes hadrus` | Higher in Control |
| `Flavonifractor plautii` | Higher in CRC |
| `Odoribacter splanchnicus` | Higher in CRC |
| `Roseburia faecis` | Higher in Control |
| `Lachnospira eligens` | Higher in Control |
| `Parabacteroides distasonis` | Higher in CRC |
| `Bacteroides fragilis` | Higher in CRC |

**Interpretation:**  
The most reproducible meta-analytic signal is **CRC vs Control**, particularly in males. Adenoma-related contrasts were not significant after FDR in the random-effects analysis.

OUTPUTS/phase_3/step3_3_meta_analysis_summary.png

---

## Step 3.4 — Visualization summary

**What we did:**  
Created visualization-oriented summaries: volcano plots, a forest plot, a heatmap by condition/sex, and a CSV of top differentially abundant species for later SHAP comparison.

**Important note:**  
Step 3.4 uses simple visualization statistics and should not be treated as the main statistical evidence. Formal evidence should come from Steps 3.1, 3.2, and 3.3.

**Visualization-level significant counts, FDR < 0.05:**

| Comparison | rCLR species | PA species |
|---|---:|---:|
| CRC vs Control | **26** | **18** |
| Adenoma vs Control | **4** | **9** |
| CRC vs Adenoma | **13** | **15** |

**Use:**  
This step is mainly useful for communicating patterns and for selecting candidate species to compare with SHAP features in Phase 5.

image:
	OUTPUTS/phase_3/step3_4_heatmap_top_species_condition_by_sex.png

---

## Step 3.5 — Adenoma-carcinoma sequence analysis

**What we did:**  
Tested whether species follow an ordered trend across:

```text
Control → Adenoma → CRC
```

using Jonckheere-Terpstra trend tests, plus pairwise checks and biomarker categorization.

**Main counts:**

| Category | Number of species |
|---|---:|
| Ordered increasing markers | **13** |
| Ordered decreasing markers | **11** |
| Early markers | **2** |
| Late markers | **24** |
| Progression markers | **13** |

**Interpretation:**  
There is some evidence for ordered microbiome shifts across Control → Adenoma → CRC, but most robust formal evidence from Steps 3.1–3.3 still points to CRC vs Control as the clearest signal. Adenoma findings should remain exploratory because Adenoma samples come from fewer studies.

images:

	OUTPUTS/phase_3/step3_5_biomarker_overlap.png
	OUTPUTS/phase_3/step3_5_effect_size_gradient_top_species.png

---

# Overall interpretation so far

1. **Do not merge Adenoma with CRC or Control yet.**  
   Adenoma is clinically distinct and the exploratory analyses do not justify forcing it into another group.

2. **The strongest and most reproducible signal is CRC vs Control.**  
   This appears in rCLR, PA, and meta-analysis results.

3. **The signal is stronger in males than females.**  
   Male-stratified CRC vs Control repeatedly produced more significant species than female-stratified analysis.

4. **Study/batch effects are present.**  
   Unsupervised clustering aligned more with study than condition, supporting Leave-One-Dataset-Out evaluation later.

5. **Adenoma results are interesting but less stable.**  
   Some visualization and trend analyses suggest possible sequence patterns, but formal blocked/meta-analytic evidence is weaker.

---

# Recommended next steps

## Immediate next step: Phase 4 ML pipeline

Proceed with multiple ML tasks instead of one forced Adenoma decision:

| Model | Task | Adenoma handling |
|---|---|---|
| Model A | CRC vs Control | Exclude Adenoma from training; score Adenoma afterward |
| Model B | Male-only CRC vs Control | Exclude Adenoma |
| Model C | Female-only CRC vs Control | Exclude Adenoma |
| Model D | Adenoma vs Control | Use Adenoma as positive class |
| Model E | CRC vs Adenoma | Directly test progression distinction |
| Model F | Any neoplasia vs Control | Merge CRC + Adenoma as positive class |

## Important ML design points

- Use **LODO** evaluation because study/batch effects are visible.
- Do not batch-correct test folds.
- Keep PA and rCLR features separate/interpretable.
- Compare SHAP features later with the top species from Phase 3.
- Report sex-stratified performance and fairness metrics for the universal CRC model.

