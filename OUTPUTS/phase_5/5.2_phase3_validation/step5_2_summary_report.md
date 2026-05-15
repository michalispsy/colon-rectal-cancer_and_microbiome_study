# Step 5.2 — Cross-validation of SHAP Features with Phase 3

Models included: **All**, **OnlyMale**, **OnlyFemale**. Adenoma / Models D,E,F were intentionally excluded.

## Inputs used

Formal Phase 3 validation sources:

- 3.1 blocked Wilcoxon / Van Elteren significant rCLR species

- 3.2 Cochran-Mantel-Haenszel significant PA species

- 3.3 random-effects meta-analysis significant rCLR species

- 3.5 adenoma-carcinoma biomarker categories were added as optional biological annotation



## Overlap summary

| DisplayModel   |   Top_N |   N_SHAP_species |   N_overlap_any_formal_Phase3 |   Percent_overlap_any_formal_Phase3 |   N_overlap_Phase3_1_BlockedWilcoxon |   N_overlap_Phase3_2_CMH |   N_overlap_Phase3_3_MetaAnalysis |   N_supported_by_2plus_formal_tests | RedFlag_No_Overlap   |
|:---------------|--------:|-----------------:|------------------------------:|------------------------------------:|-------------------------------------:|-------------------------:|----------------------------------:|------------------------------------:|:---------------------|
| All            |      10 |               10 |                             4 |                                  40 |                                    2 |                        2 |                                 1 |                                   1 | False                |
| All            |      20 |               20 |                             9 |                                  45 |                                    4 |                        5 |                                 2 |                                   2 | False                |
| All            |      50 |               50 |                            20 |                                  40 |                                   12 |                       12 |                                 5 |                                   8 | False                |
| OnlyMale       |      10 |               10 |                            10 |                                 100 |                                    8 |                        3 |                                 7 |                                   8 | False                |
| OnlyMale       |      20 |               20 |                            14 |                                  70 |                                    8 |                        7 |                                 7 |                                   8 | False                |
| OnlyMale       |      50 |               50 |                            23 |                                  46 |                                    9 |                       15 |                                 7 |                                   8 | False                |
| OnlyFemale     |      10 |               10 |                             1 |                                  10 |                                    0 |                        1 |                                 0 |                                   0 | False                |
| OnlyFemale     |      20 |               20 |                             1 |                                   5 |                                    0 |                        1 |                                 0 |                                   0 | False                |
| OnlyFemale     |      50 |               50 |                             1 |                                   2 |                                    0 |                        1 |                                 0 |                                   0 | False                |


## All

Top 20 SHAP species with Phase 3 support:

| Species_clean                |   SpeciesRank |   Species_MeanAbsSHAP | Dominant_Signal   | In_Phase3_1_BlockedWilcoxon_rCLR   | In_Phase3_2_CMH_PA   | In_Phase3_3_MetaAnalysis_rCLR   |   N_Formal_Phase3_Supports | Phase3_Support_Level                                   | categories              | is_early   | is_late   | is_progression   |
|:-----------------------------|--------------:|----------------------:|:------------------|:-----------------------------------|:---------------------|:--------------------------------|---------------------------:|:-------------------------------------------------------|:------------------------|:-----------|:----------|:-----------------|
| Parvimonas micra             |             1 |              1.32678  | PA                | True                               | True                 | False                           |                          2 | Strong: supported by at least two formal Phase 3 tests | LATE                    | False      | True      | False            |
| Streptococcus salivarius     |             2 |              0.980825 | PA                | False                              | True                 | False                           |                          1 | Moderate: supported by one formal Phase 3 test         | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Veillonella parvula          |             3 |              0.61825  | PA                | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | ORDERED_INCREASING      | False      | False     | False            |
| Firmicutes bacterium CAG:110 |             4 |              0.502097 | PA                | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Streptococcus thermophilus   |             5 |              0.415792 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | PROGRESSION             | False      | False     | True             |
| Klebsiella variicola         |             6 |              0.386114 | PA                | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Dorea formicigenerans        |             7 |              0.358611 | rCLR              | True                               | False                | False                           |                          1 | Moderate: supported by one formal Phase 3 test         | LATE;ORDERED_INCREASING | False      | True      | False            |
| Eggerthella lenta            |             8 |              0.357058 | PA                | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Flavonifractor plautii       |             9 |              0.326901 | rCLR              | False                              | False                | True                            |                          1 | Strong: reproduced by random-effects meta-analysis     | LATE;PROGRESSION        | False      | True      | True             |
| Bacteroides intestinalis     |            10 |              0.324041 | PA                | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Anaerostipes hadrus          |            11 |              0.320847 | rCLR              | True                               | False                | True                            |                          2 | Strong: reproduced by random-effects meta-analysis     | LATE;ORDERED_DECREASING | False      | True      | False            |
| Parabacteroides merdae       |            12 |              0.319956 | PA                | False                              | True                 | False                           |                          1 | Moderate: supported by one formal Phase 3 test         | PROGRESSION             | False      | False     | True             |
| Butyricimonas virosa         |            13 |              0.311247 | rCLR              | False                              | True                 | False                           |                          1 | Moderate: supported by one formal Phase 3 test         | LATE;ORDERED_INCREASING | False      | True      | False            |
| Klebsiella pneumoniae        |            14 |              0.307455 | PA                | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Streptococcus parasanguinis  |            15 |              0.283843 | PA                | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Coprococcus comes            |            16 |              0.281853 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | LATE;ORDERED_DECREASING | False      | True      | False            |
| Methanobrevibacter smithii   |            17 |              0.269413 | PA                | False                              | True                 | False                           |                          1 | Moderate: supported by one formal Phase 3 test         | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| [Ruminococcus] torques       |            18 |              0.266773 | rCLR              | True                               | False                | False                           |                          1 | Moderate: supported by one formal Phase 3 test         | LATE                    | False      | True      | False            |
| [Eubacterium] siraeum        |            19 |              0.26426  | PA                | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Gemmiger formicilis          |            20 |              0.264156 | PA                | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |


> Interpretation: 9/20 top SHAP species overlap with formal Phase 3 results, including 2 supported by meta-analysis.


## OnlyMale

Top 20 SHAP species with Phase 3 support:

| Species_clean                     |   SpeciesRank |   Species_MeanAbsSHAP | Dominant_Signal   | In_Phase3_1_BlockedWilcoxon_rCLR   | In_Phase3_2_CMH_PA   | In_Phase3_3_MetaAnalysis_rCLR   |   N_Formal_Phase3_Supports | Phase3_Support_Level                                   | categories              | is_early   | is_late   | is_progression   |
|:----------------------------------|--------------:|----------------------:|:------------------|:-----------------------------------|:---------------------|:--------------------------------|---------------------------:|:-------------------------------------------------------|:------------------------|:-----------|:----------|:-----------------|
| Parvimonas micra                  |             1 |            0.0508442  | PA                | True                               | True                 | False                           |                          2 | Strong: supported by at least two formal Phase 3 tests | LATE                    | False      | True      | False            |
| Roseburia intestinalis            |             2 |            0.0120392  | rCLR              | True                               | True                 | False                           |                          2 | Strong: supported by at least two formal Phase 3 tests | EARLY                   | True       | False     | False            |
| Roseburia faecis                  |             3 |            0.0119566  | rCLR              | True                               | False                | True                            |                          2 | Strong: reproduced by random-effects meta-analysis     | LATE;ORDERED_DECREASING | False      | True      | False            |
| Ruthenibacterium lactatiformans   |             4 |            0.0118855  | rCLR              | False                              | False                | True                            |                          1 | Strong: reproduced by random-effects meta-analysis     | LATE;ORDERED_INCREASING | False      | True      | False            |
| Odoribacter splanchnicus          |             5 |            0.011228   | rCLR              | True                               | False                | True                            |                          2 | Strong: reproduced by random-effects meta-analysis     | LATE;PROGRESSION        | False      | True      | True             |
| Akkermansia muciniphila           |             6 |            0.0112236  | rCLR              | False                              | True                 | False                           |                          1 | Moderate: supported by one formal Phase 3 test         | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Fusicatenibacter saccharivorans   |             7 |            0.0112116  | rCLR              | True                               | False                | True                            |                          2 | Strong: reproduced by random-effects meta-analysis     | LATE;ORDERED_DECREASING | False      | True      | False            |
| Lachnospira eligens               |             8 |            0.0107274  | rCLR              | True                               | False                | True                            |                          2 | Strong: reproduced by random-effects meta-analysis     | EARLY                   | True       | False     | False            |
| Eubacterium ventriosum            |             9 |            0.00893974 | rCLR              | True                               | False                | True                            |                          2 | Strong: reproduced by random-effects meta-analysis     | LATE;ORDERED_DECREASING | False      | True      | False            |
| Lachnospira pectinoschiza         |            10 |            0.00818904 | rCLR              | True                               | False                | True                            |                          2 | Strong: reproduced by random-effects meta-analysis     | LATE                    | False      | True      | False            |
| Firmicutes bacterium CAG:145      |            11 |            0.00806158 | rCLR              | False                              | True                 | False                           |                          1 | Moderate: supported by one formal Phase 3 test         | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Bacteroides caccae                |            12 |            0.00748336 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | LATE                    | False      | True      | False            |
| Methanobrevibacter smithii        |            13 |            0.00736163 | rCLR              | False                              | True                 | False                           |                          1 | Moderate: supported by one formal Phase 3 test         | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Bifidobacterium pseudocatenulatum |            14 |            0.00703513 | rCLR              | False                              | True                 | False                           |                          1 | Moderate: supported by one formal Phase 3 test         | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Eubacterium sp. CAG:251           |            15 |            0.00684197 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | LATE;ORDERED_DECREASING | False      | True      | False            |
| Paraprevotella xylaniphila        |            16 |            0.00666658 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | LATE                    | False      | True      | False            |
| Coprococcus catus                 |            17 |            0.00611124 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Collinsella aerofaciens           |            18 |            0.00580356 | rCLR              | False                              | True                 | False                           |                          1 | Moderate: supported by one formal Phase 3 test         | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Phocaeicola coprocola             |            19 |            0.00554209 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Eubacterium sp. CAG:38            |            20 |            0.00529324 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap                     | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |


> Interpretation: 14/20 top SHAP species overlap with formal Phase 3 results, including 7 supported by meta-analysis.


## OnlyFemale

Top 20 SHAP species with Phase 3 support:

| Species_clean                   |   SpeciesRank |   Species_MeanAbsSHAP | Dominant_Signal   | In_Phase3_1_BlockedWilcoxon_rCLR   | In_Phase3_2_CMH_PA   | In_Phase3_3_MetaAnalysis_rCLR   |   N_Formal_Phase3_Supports | Phase3_Support_Level                           | categories              | is_early   | is_late   | is_progression   |
|:--------------------------------|--------------:|----------------------:|:------------------|:-----------------------------------|:---------------------|:--------------------------------|---------------------------:|:-----------------------------------------------|:------------------------|:-----------|:----------|:-----------------|
| Ruthenibacterium lactatiformans |             1 |            0.0229693  | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | LATE;ORDERED_INCREASING | False      | True      | False            |
| Parvimonas micra                |             2 |            0.0223346  | PA                | False                              | True                 | False                           |                          1 | Moderate: supported by one formal Phase 3 test | LATE                    | False      | True      | False            |
| Anaerostipes hadrus             |             3 |            0.0134342  | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | LATE;ORDERED_DECREASING | False      | True      | False            |
| Bacteroides fragilis            |             4 |            0.0126423  | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | LATE;ORDERED_INCREASING | False      | True      | False            |
| Dorea formicigenerans           |             5 |            0.0106012  | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | LATE;ORDERED_INCREASING | False      | True      | False            |
| Bacteroides thetaiotaomicron    |             6 |            0.0100292  | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | LATE;PROGRESSION        | False      | True      | True             |
| Flavonifractor plautii          |             7 |            0.00881873 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | LATE;PROGRESSION        | False      | True      | True             |
| Anaerobutyricum hallii          |             8 |            0.00858661 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | PROGRESSION             | False      | False     | True             |
| Parabacteroides distasonis      |             9 |            0.00808864 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | LATE;PROGRESSION        | False      | True      | True             |
| [Ruminococcus] torques          |            10 |            0.00730533 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | LATE                    | False      | True      | False            |
| Coprococcus comes               |            11 |            0.00720315 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | LATE;ORDERED_DECREASING | False      | True      | False            |
| Phocaeicola plebeius            |            12 |            0.00623292 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Bacteroides caccae              |            13 |            0.00595898 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | LATE                    | False      | True      | False            |
| Alistipes finegoldii            |            14 |            0.00592141 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | LATE;ORDERED_INCREASING | False      | True      | False            |
| [Ruminococcus] gnavus           |            15 |            0.00560753 | PA                | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Adlercreutzia equolifaciens     |            16 |            0.00548962 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |
| Eubacterium sp. CAG:251         |            17 |            0.00498512 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | LATE;ORDERED_DECREASING | False      | True      | False            |
| Parabacteroides merdae          |            18 |            0.00496844 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | PROGRESSION             | False      | False     | True             |
| Eubacterium sp. CAG:180         |            19 |            0.00496682 | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | ORDERED_INCREASING      | False      | False     | False            |
| Bacteroides sp. CAG:144         |            20 |            0.004918   | rCLR              | False                              | False                | False                           |                          0 | ML-only: no formal Phase 3 overlap             | NOT_SIGNIFICANT_BY_FDR  | False      | False     | False            |


> Interpretation: 1/20 top SHAP species overlap with formal Phase 3 results, but none with the meta-analysis core list.


## Phase 3 reference set sizes

| DisplayModel   |   Phase3Step |   N_species |
|:---------------|-------------:|------------:|
| All            |          3.1 |          15 |
| All            |          3.2 |          19 |
| All            |          3.3 |           9 |
| OnlyFemale     |          3.2 |           1 |
| OnlyMale       |          3.1 |           9 |
| OnlyMale       |          3.2 |          18 |
| OnlyMale       |          3.3 |           7 |

