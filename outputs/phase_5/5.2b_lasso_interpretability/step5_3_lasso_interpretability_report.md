# Step 5.3 — LASSO Logistic Regression Interpretability


## Purpose

This analysis explores interpretability using an L1-penalized Logistic Regression model. The L1 penalty induces sparsity by shrinking some feature coefficients exactly to zero. Therefore, selected non-zero features form a compact and directly interpretable microbial signature.


## Arrangements and model

- Arrangements: **All**, **OnlyMale**, **OnlyFemale**.

- Model: **LASSO Logistic Regression**.

- Validation: same LODO folds as Phase 4.

- C tuning inside train folds: **no**.


## Performance and sparsity summary

| Arrangement   |   N_folds |   Median_AUC |   IQR_AUC_low |   IQR_AUC_high |   Median_AUPRC |   Median_F1 |   Median_N_nonzero_features |   IQR_N_nonzero_low |   IQR_N_nonzero_high |   Median_Percent_nonzero_features |   Median_Selected_C |
|:--------------|----------:|-------------:|--------------:|---------------:|---------------:|------------:|----------------------------:|--------------------:|---------------------:|----------------------------------:|--------------------:|
| All           |        10 |     0.734635 |      0.668503 |       0.797041 |       0.745024 |    0.646586 |                        31   |               28    |                   34 |                          14.3519  |                0.03 |
| OnlyFemale    |        10 |     0.722525 |      0.608392 |       0.747837 |       0.727324 |    0.516667 |                         5.5 |                5    |                    6 |                           2.55814 |                0.03 |
| OnlyMale      |        10 |     0.747209 |      0.685343 |       0.81964  |       0.803915 |    0.642857 |                        14.5 |               12.25 |                   15 |                           6.74419 |                0.03 |



## All

Top consensus coefficients:

|   ConsensusRank | Feature                              | Feature_Type   | Species_clean                   | ConsensusDirection   |   SelectionFrequency |   N_Folds_NonZero |   MeanNonZeroCoefficient |   MeanAbsNonZeroCoefficient |   ConsensusScore |
|----------------:|:-------------------------------------|:---------------|:--------------------------------|:---------------------|---------------------:|------------------:|-------------------------:|----------------------------:|-----------------:|
|               1 | PA_Parvimonas micra                  | PA             | Parvimonas micra                | CRC-associated       |                  1   |                10 |                0.524638  |                   0.524638  |         1        |
|               2 | PA_Streptococcus salivarius          | PA             | Streptococcus salivarius        | Control-associated   |                  1   |                10 |               -0.155672  |                   0.155672  |         0.824181 |
|               3 | rCLR_Ruthenibacterium lactatiformans | rCLR           | Ruthenibacterium lactatiformans | CRC-associated       |                  1   |                10 |                0.148245  |                   0.148245  |         0.820642 |
|               4 | rCLR_Anaerostipes hadrus             | rCLR           | Anaerostipes hadrus             | Control-associated   |                  1   |                10 |               -0.129264  |                   0.129264  |         0.811597 |
|               5 | Age                                  | covariate      | Age                             | CRC-associated       |                  1   |                10 |                0.102443  |                   0.102443  |         0.798816 |
|               6 | rCLR_Dorea formicigenerans           | rCLR           | Dorea formicigenerans           | CRC-associated       |                  1   |                10 |                0.0792405 |                   0.0792405 |         0.78776  |
|               7 | rCLR_Odoribacter splanchnicus        | rCLR           | Odoribacter splanchnicus        | CRC-associated       |                  1   |                10 |                0.0642428 |                   0.0642428 |         0.780613 |
|               8 | PA_Roseburia intestinalis            | PA             | Roseburia intestinalis          | Control-associated   |                  1   |                10 |               -0.0445104 |                   0.0445104 |         0.77121  |
|               9 | rCLR_[Ruminococcus] torques          | rCLR           | [Ruminococcus] torques          | CRC-associated       |                  0.9 |                 9 |                0.0998441 |                   0.0998441 |         0.742578 |
|              10 | PA_Phocaeicola plebeius              | PA             | Phocaeicola plebeius            | CRC-associated       |                  0.9 |                 9 |                0.0905485 |                   0.0905485 |         0.738148 |
|              11 | rCLR_Streptococcus thermophilus      | rCLR           | Streptococcus thermophilus      | Control-associated   |                  0.9 |                 9 |               -0.0808424 |                   0.0808424 |         0.733523 |
|              12 | rCLR_Eubacterium ventriosum          | rCLR           | Eubacterium ventriosum          | Control-associated   |                  0.9 |                 9 |               -0.0730187 |                   0.0730187 |         0.729795 |
|              13 | rCLR_Flavonifractor plautii          | rCLR           | Flavonifractor plautii          | CRC-associated       |                  0.9 |                 9 |                0.0492789 |                   0.0492789 |         0.718482 |
|              14 | PA_Anaerotignum lactatifermentans    | PA             | Anaerotignum lactatifermentans  | CRC-associated       |                  0.9 |                 9 |                0.0454713 |                   0.0454713 |         0.716668 |
|              15 | PA_Methanobrevibacter smithii        | PA             | Methanobrevibacter smithii      | CRC-associated       |                  0.9 |                 9 |                0.0415347 |                   0.0415347 |         0.714792 |
|              16 | rCLR_Lachnospira eligens             | rCLR           | Lachnospira eligens             | Control-associated   |                  0.9 |                 9 |               -0.0342009 |                   0.0342009 |         0.711297 |
|              17 | PA_Eubacterium sp. CAG:38            | PA             | Eubacterium sp. CAG:38          | Control-associated   |                  0.9 |                 9 |               -0.0225927 |                   0.0225927 |         0.705766 |
|              18 | rCLR_Paraprevotella xylaniphila      | rCLR           | Paraprevotella xylaniphila      | Control-associated   |                  0.8 |                 8 |               -0.0532147 |                   0.0532147 |         0.665358 |
|              19 | rCLR_Lachnospira pectinoschiza       | rCLR           | Lachnospira pectinoschiza       | Control-associated   |                  0.8 |                 8 |               -0.0333872 |                   0.0333872 |         0.65591  |
|              20 | rCLR_Bacteroides thetaiotaomicron    | rCLR           | Bacteroides thetaiotaomicron    | CRC-associated       |                  0.8 |                 8 |                0.0271206 |                   0.0271206 |         0.652923 |


Feature-type split among top coefficients:

| Feature_Type   |   N_Features |   MeanAbsNonZeroCoef_Total |   Mean_SelectionFrequency |   Percent_TopN |   Percent_Coefficient_TopN |
|:---------------|-------------:|---------------------------:|--------------------------:|---------------:|---------------------------:|
| PA             |            7 |                   0.924967 |                  0.942857 |             35 |                    48.7002 |
| rCLR           |           12 |                   0.8719   |                  0.908333 |             60 |                    45.9061 |
| covariate      |            1 |                   0.102443 |                  1        |              5 |                     5.3937 |


Top species-level LASSO consensus:

|   SpeciesRank | Species_clean                   | Dominant_Signal   |   Species_MaxSelectionFrequency |   Species_SumMeanAbsNonZeroCoef | Any_CRC_Associated   | Any_Control_Associated   |
|--------------:|:--------------------------------|:------------------|--------------------------------:|--------------------------------:|:---------------------|:-------------------------|
|             1 | Parvimonas micra                | PA                |                             1   |                       0.538573  | True                 | False                    |
|             2 | Streptococcus salivarius        | PA                |                             1   |                       0.171854  | False                | True                     |
|             3 | Ruthenibacterium lactatiformans | rCLR              |                             1   |                       0.158066  | True                 | False                    |
|             4 | Anaerostipes hadrus             | rCLR              |                             1   |                       0.129264  | False                | True                     |
|             5 | Dorea formicigenerans           | rCLR              |                             1   |                       0.0792405 | True                 | False                    |
|             6 | Odoribacter splanchnicus        | rCLR              |                             1   |                       0.0642428 | True                 | False                    |
|             7 | Roseburia intestinalis          | PA                |                             1   |                       0.0624305 | False                | True                     |
|             8 | [Ruminococcus] torques          | rCLR              |                             0.9 |                       0.0998441 | True                 | False                    |
|             9 | Streptococcus thermophilus      | rCLR              |                             0.9 |                       0.0956802 | False                | True                     |
|            10 | Phocaeicola plebeius            | PA                |                             0.9 |                       0.0905485 | True                 | False                    |
|            11 | Flavonifractor plautii          | rCLR              |                             0.9 |                       0.0811154 | True                 | False                    |
|            12 | Eubacterium ventriosum          | rCLR              |                             0.9 |                       0.0730187 | False                | True                     |
|            13 | Lachnospira eligens             | rCLR              |                             0.9 |                       0.0526162 | False                | True                     |
|            14 | Anaerotignum lactatifermentans  | PA                |                             0.9 |                       0.0454713 | True                 | False                    |
|            15 | Methanobrevibacter smithii      | PA                |                             0.9 |                       0.0415347 | True                 | False                    |
|            16 | Eubacterium sp. CAG:38          | PA                |                             0.9 |                       0.0225927 | False                | True                     |
|            17 | Paraprevotella xylaniphila      | rCLR              |                             0.8 |                       0.0532147 | False                | True                     |
|            18 | Lachnospira pectinoschiza       | rCLR              |                             0.8 |                       0.0333872 | False                | True                     |
|            19 | Bacteroides thetaiotaomicron    | rCLR              |                             0.8 |                       0.0271206 | True                 | False                    |
|            20 | Parabacteroides distasonis      | rCLR              |                             0.8 |                       0.0221499 | True                 | False                    |



LASSO vs Random Forest SHAP: 13/20 top LASSO species are also top RF SHAP species.

LASSO vs Phase 3: 17/20 top LASSO species overlap with formal Phase 3 results; 7/20 overlap with random-effects meta-analysis.


## OnlyMale

Top consensus coefficients:

|   ConsensusRank | Feature                              | Feature_Type   | Species_clean                     | ConsensusDirection   |   SelectionFrequency |   N_Folds_NonZero |   MeanNonZeroCoefficient |   MeanAbsNonZeroCoefficient |   ConsensusScore |
|----------------:|:-------------------------------------|:---------------|:----------------------------------|:---------------------|---------------------:|------------------:|-------------------------:|----------------------------:|-----------------:|
|               1 | PA_Parvimonas micra                  | PA             | Parvimonas micra                  | CRC-associated       |                  1   |                10 |               0.502642   |                  0.502642   |         1        |
|               2 | PA_Streptococcus salivarius          | PA             | Streptococcus salivarius          | Control-associated   |                  1   |                10 |              -0.116122   |                  0.116122   |         0.807756 |
|               3 | rCLR_Ruthenibacterium lactatiformans | rCLR           | Ruthenibacterium lactatiformans   | CRC-associated       |                  1   |                10 |               0.115845   |                  0.115845   |         0.807618 |
|               4 | rCLR_Eubacterium ventriosum          | rCLR           | Eubacterium ventriosum            | Control-associated   |                  1   |                10 |              -0.0942184  |                  0.0942184  |         0.796862 |
|               5 | rCLR_Lachnospira eligens             | rCLR           | Lachnospira eligens               | Control-associated   |                  1   |                10 |              -0.0573831  |                  0.0573831  |         0.778541 |
|               6 | rCLR_Odoribacter splanchnicus        | rCLR           | Odoribacter splanchnicus          | CRC-associated       |                  0.9 |                 9 |               0.060418   |                  0.060418   |         0.72505  |
|               7 | rCLR_Roseburia intestinalis          | rCLR           | Roseburia intestinalis            | Control-associated   |                  0.9 |                 9 |              -0.0601584  |                  0.0601584  |         0.724921 |
|               8 | Age                                  | covariate      | Age                               | CRC-associated       |                  0.9 |                 9 |               0.0471962  |                  0.0471962  |         0.718474 |
|               9 | PA_Bifidobacterium bifidum           | PA             | Bifidobacterium bifidum           | Control-associated   |                  0.8 |                 8 |              -0.0124031  |                  0.0124031  |         0.646169 |
|              10 | PA_Collinsella aerofaciens           | PA             | Collinsella aerofaciens           | CRC-associated       |                  0.7 |                 7 |               0.034187   |                  0.034187   |         0.602004 |
|              11 | rCLR_Paraprevotella xylaniphila      | rCLR           | Paraprevotella xylaniphila        | Control-associated   |                  0.6 |                 6 |              -0.0252189  |                  0.0252189  |         0.542543 |
|              12 | rCLR_Fusicatenibacter saccharivorans | rCLR           | Fusicatenibacter saccharivorans   | Control-associated   |                  0.6 |                 6 |              -0.0199359  |                  0.0199359  |         0.539916 |
|              13 | PA_Anaerotignum lactatifermentans    | PA             | Anaerotignum lactatifermentans    | CRC-associated       |                  0.6 |                 6 |               0.0163788  |                  0.0163788  |         0.538146 |
|              14 | PA_Bifidobacterium pseudocatenulatum | PA             | Bifidobacterium pseudocatenulatum | Control-associated   |                  0.4 |                 4 |              -0.0158674  |                  0.0158674  |         0.427892 |
|              15 | rCLR_Lachnospira pectinoschiza       | rCLR           | Lachnospira pectinoschiza         | Control-associated   |                  0.4 |                 4 |              -0.0153973  |                  0.0153973  |         0.427658 |
|              16 | PA_Methanobrevibacter smithii        | PA             | Methanobrevibacter smithii        | CRC-associated       |                  0.3 |                 3 |               0.0263453  |                  0.0263453  |         0.378103 |
|              17 | PA_Bacteroides salyersiae            | PA             | Bacteroides salyersiae            | CRC-associated       |                  0.3 |                 3 |               0.0095758  |                  0.0095758  |         0.369763 |
|              18 | rCLR_Bacteroides caccae              | rCLR           | Bacteroides caccae                | CRC-associated       |                  0.3 |                 3 |               0.00723459 |                  0.00723459 |         0.368598 |
|              19 | rCLR_Roseburia faecis                | rCLR           | Roseburia faecis                  | Control-associated   |                  0.2 |                 2 |              -0.00647039 |                  0.00647039 |         0.313218 |
|              20 | PA_Akkermansia muciniphila           | PA             | Akkermansia muciniphila           | CRC-associated       |                  0.1 |                 1 |               0.0460007  |                  0.0460007  |         0.277879 |


Feature-type split among top coefficients:

| Feature_Type   |   N_Features |   MeanAbsNonZeroCoef_Total |   Mean_SelectionFrequency |   Percent_TopN |   Percent_Coefficient_TopN |
|:---------------|-------------:|---------------------------:|--------------------------:|---------------:|---------------------------:|
| PA             |            9 |                  0.779523  |                  0.577778 |             45 |                   60.4751  |
| rCLR           |           10 |                  0.46228   |                  0.69     |             50 |                   35.8635  |
| covariate      |            1 |                  0.0471962 |                  0.9      |              5 |                    3.66146 |


Top species-level LASSO consensus:

|   SpeciesRank | Species_clean                     | Dominant_Signal   |   Species_MaxSelectionFrequency |   Species_SumMeanAbsNonZeroCoef | Any_CRC_Associated   | Any_Control_Associated   |
|--------------:|:----------------------------------|:------------------|--------------------------------:|--------------------------------:|:---------------------|:-------------------------|
|             1 | Parvimonas micra                  | PA                |                             1   |                      0.502642   | True                 | False                    |
|             2 | Streptococcus salivarius          | PA                |                             1   |                      0.116122   | False                | True                     |
|             3 | Ruthenibacterium lactatiformans   | rCLR              |                             1   |                      0.115845   | True                 | False                    |
|             4 | Eubacterium ventriosum            | rCLR              |                             1   |                      0.0942184  | False                | True                     |
|             5 | Lachnospira eligens               | rCLR              |                             1   |                      0.0573831  | False                | True                     |
|             6 | Odoribacter splanchnicus          | rCLR              |                             0.9 |                      0.060418   | True                 | False                    |
|             7 | Roseburia intestinalis            | rCLR              |                             0.9 |                      0.0601584  | False                | True                     |
|             8 | Bifidobacterium bifidum           | PA                |                             0.8 |                      0.0124031  | False                | True                     |
|             9 | Collinsella aerofaciens           | PA                |                             0.7 |                      0.034187   | True                 | False                    |
|            10 | Paraprevotella xylaniphila        | rCLR              |                             0.6 |                      0.0252189  | False                | True                     |
|            11 | Fusicatenibacter saccharivorans   | rCLR              |                             0.6 |                      0.0199359  | False                | True                     |
|            12 | Anaerotignum lactatifermentans    | PA                |                             0.6 |                      0.0163788  | True                 | False                    |
|            13 | Bifidobacterium pseudocatenulatum | PA                |                             0.4 |                      0.0158674  | False                | True                     |
|            14 | Lachnospira pectinoschiza         | rCLR              |                             0.4 |                      0.0153973  | False                | True                     |
|            15 | Methanobrevibacter smithii        | PA                |                             0.3 |                      0.0263453  | True                 | False                    |
|            16 | Bacteroides salyersiae            | PA                |                             0.3 |                      0.0095758  | True                 | False                    |
|            17 | Bacteroides caccae                | rCLR              |                             0.3 |                      0.00723459 | True                 | False                    |
|            18 | Roseburia faecis                  | rCLR              |                             0.2 |                      0.00647039 | False                | True                     |
|            19 | Akkermansia muciniphila           | PA                |                             0.1 |                      0.0460007  | True                 | False                    |
|            20 | Anaerostipes hadrus               | rCLR              |                             0.1 |                      0.0118053  | False                | True                     |



LASSO vs Random Forest SHAP: 15/20 top LASSO species are also top RF SHAP species.

LASSO vs Phase 3: 17/20 top LASSO species overlap with formal Phase 3 results; 7/20 overlap with random-effects meta-analysis.


## OnlyFemale

Top consensus coefficients:

|   ConsensusRank | Feature                              | Feature_Type   | Species_clean                   | ConsensusDirection   |   SelectionFrequency |   N_Folds_NonZero |   MeanNonZeroCoefficient |   MeanAbsNonZeroCoefficient |   ConsensusScore |
|----------------:|:-------------------------------------|:---------------|:--------------------------------|:---------------------|---------------------:|------------------:|-------------------------:|----------------------------:|-----------------:|
|               1 | PA_Parvimonas micra                  | PA             | Parvimonas micra                | CRC-associated       |                  1   |                10 |               0.256146   |                  0.256146   |         1        |
|               2 | rCLR_Ruthenibacterium lactatiformans | rCLR           | Ruthenibacterium lactatiformans | CRC-associated       |                  0.9 |                 9 |               0.0329679  |                  0.0329679  |         0.727177 |
|               3 | rCLR_Dorea formicigenerans           | rCLR           | Dorea formicigenerans           | CRC-associated       |                  0.8 |                 8 |               0.0573599  |                  0.0573599  |         0.695984 |
|               4 | rCLR_Anaerostipes hadrus             | rCLR           | Anaerostipes hadrus             | Control-associated   |                  0.8 |                 8 |              -0.0535673  |                  0.0535673  |         0.692282 |
|               5 | PA_Flavonifractor plautii            | PA             | Flavonifractor plautii          | CRC-associated       |                  0.8 |                 8 |               0.0480677  |                  0.0480677  |         0.686914 |
|               6 | rCLR_Bacteroides fragilis            | rCLR           | Bacteroides fragilis            | CRC-associated       |                  0.6 |                 6 |               0.0212347  |                  0.0212347  |         0.550725 |
|               7 | rCLR_Anaerobutyricum hallii          | rCLR           | Anaerobutyricum hallii          | Control-associated   |                  0.3 |                 3 |              -0.0166958  |                  0.0166958  |         0.381295 |
|               8 | PA_Phocaeicola plebeius              | PA             | Phocaeicola plebeius            | CRC-associated       |                  0.2 |                 2 |               0.00883994 |                  0.00883994 |         0.318628 |
|               9 | PA_[Ruminococcus] gnavus             | PA             | [Ruminococcus] gnavus           | CRC-associated       |                  0.1 |                 1 |               0.0288142  |                  0.0288142  |         0.283123 |
|              10 | rCLR_Bacteroides thetaiotaomicron    | rCLR           | Bacteroides thetaiotaomicron    | CRC-associated       |                  0.1 |                 1 |               0.0163057  |                  0.0163057  |         0.270914 |
|              11 | Age                                  | covariate      | Age                             | mixed/zero           |                  0   |                 0 |               0          |                  0          |         0        |
|              12 | PA_Acidaminococcus intestini         | PA             | Acidaminococcus intestini       | mixed/zero           |                  0   |                 0 |               0          |                  0          |         0        |
|              13 | PA_Adlercreutzia equolifaciens       | PA             | Adlercreutzia equolifaciens     | mixed/zero           |                  0   |                 0 |               0          |                  0          |         0        |
|              14 | PA_Adlercreutzia equolifaciens_1     | PA             | Adlercreutzia equolifaciens 1   | mixed/zero           |                  0   |                 0 |               0          |                  0          |         0        |
|              15 | PA_Agathobaculum butyriciproducens   | PA             | Agathobaculum butyriciproducens | mixed/zero           |                  0   |                 0 |               0          |                  0          |         0        |
|              16 | PA_Akkermansia muciniphila           | PA             | Akkermansia muciniphila         | mixed/zero           |                  0   |                 0 |               0          |                  0          |         0        |
|              17 | PA_Alistipes finegoldii              | PA             | Alistipes finegoldii            | mixed/zero           |                  0   |                 0 |               0          |                  0          |         0        |
|              18 | PA_Alistipes inops                   | PA             | Alistipes inops                 | mixed/zero           |                  0   |                 0 |               0          |                  0          |         0        |
|              19 | PA_Alistipes putredinis              | PA             | Alistipes putredinis            | mixed/zero           |                  0   |                 0 |               0          |                  0          |         0        |
|              20 | PA_Alistipes shahii                  | PA             | Alistipes shahii                | mixed/zero           |                  0   |                 0 |               0          |                  0          |         0        |


Feature-type split among top coefficients:

| Feature_Type   |   N_Features |   MeanAbsNonZeroCoef_Total |   Mean_SelectionFrequency |   Percent_TopN |   Percent_Coefficient_TopN |
|:---------------|-------------:|---------------------------:|--------------------------:|---------------:|---------------------------:|
| PA             |           13 |                   0.341868 |                  0.161538 |             65 |                     63.309 |
| rCLR           |            6 |                   0.198131 |                  0.583333 |             30 |                     36.691 |
| covariate      |            1 |                   0        |                  0        |              5 |                      0     |


Top species-level LASSO consensus:

|   SpeciesRank | Species_clean                   | Dominant_Signal   |   Species_MaxSelectionFrequency |   Species_SumMeanAbsNonZeroCoef | Any_CRC_Associated   | Any_Control_Associated   |
|--------------:|:--------------------------------|:------------------|--------------------------------:|--------------------------------:|:---------------------|:-------------------------|
|             1 | Parvimonas micra                | PA                |                             1   |                      0.256146   | True                 | False                    |
|             2 | Ruthenibacterium lactatiformans | rCLR              |                             0.9 |                      0.0329679  | True                 | False                    |
|             3 | Dorea formicigenerans           | rCLR              |                             0.8 |                      0.0573599  | True                 | False                    |
|             4 | Anaerostipes hadrus             | rCLR              |                             0.8 |                      0.0535673  | False                | True                     |
|             5 | Flavonifractor plautii          | PA                |                             0.8 |                      0.0480677  | True                 | False                    |
|             6 | Bacteroides fragilis            | rCLR              |                             0.6 |                      0.0212347  | True                 | False                    |
|             7 | Anaerobutyricum hallii          | rCLR              |                             0.3 |                      0.0166958  | False                | True                     |
|             8 | Phocaeicola plebeius            | PA                |                             0.2 |                      0.00883994 | True                 | False                    |
|             9 | [Ruminococcus] gnavus           | PA                |                             0.1 |                      0.0288142  | True                 | False                    |
|            10 | Bacteroides thetaiotaomicron    | rCLR              |                             0.1 |                      0.0163057  | True                 | False                    |
|            11 | Acidaminococcus intestini       | PA                |                             0   |                      0          | False                | False                    |
|            12 | Adlercreutzia equolifaciens     | PA                |                             0   |                      0          | False                | False                    |
|            13 | Adlercreutzia equolifaciens 1   | PA                |                             0   |                      0          | False                | False                    |
|            14 | Agathobaculum butyriciproducens | PA                |                             0   |                      0          | False                | False                    |
|            15 | Akkermansia muciniphila         | PA                |                             0   |                      0          | False                | False                    |
|            16 | Alistipes finegoldii            | PA                |                             0   |                      0          | False                | False                    |
|            17 | Alistipes inops                 | PA                |                             0   |                      0          | False                | False                    |
|            18 | Alistipes putredinis            | PA                |                             0   |                      0          | False                | False                    |
|            19 | Alistipes shahii                | PA                |                             0   |                      0          | False                | False                    |
|            20 | Anaerotignum lactatifermentans  | PA                |                             0   |                      0          | False                | False                    |



LASSO vs Random Forest SHAP: 12/20 top LASSO species are also top RF SHAP species.

LASSO vs Phase 3: 1/20 top LASSO species overlap with formal Phase 3 results; 0/20 overlap with random-effects meta-analysis.


## Arrangement overlap

|   Top_N |   N_All_Top |   N_OnlyMale_Top |   N_OnlyFemale_Top |   N_Shared_All_Three |   N_Shared_All_OnlyMale |   N_Shared_All_OnlyFemale |   N_Shared_OnlyMale_OnlyFemale |   Jaccard_All_OnlyMale |   Jaccard_All_OnlyFemale |   Jaccard_OnlyMale_OnlyFemale |
|--------:|------------:|-----------------:|-------------------:|---------------------:|------------------------:|--------------------------:|-------------------------------:|-----------------------:|-------------------------:|------------------------------:|
|      20 |          20 |               20 |                 20 |                    4 |                      12 |                         8 |                              5 |               0.428571 |                     0.25 |                      0.142857 |

