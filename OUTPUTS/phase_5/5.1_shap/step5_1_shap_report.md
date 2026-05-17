# Step 5.1 — SHAP Analysis per Model
## What was computed
- Fold-level SHAP values on each held-out LODO test set.
- Top features per fold.
- Consensus ranking across folds using stability and mean absolute SHAP.
- PA vs rCLR contribution split.
- OnlyMale vs OnlyFemale overlap analysis.

## All (RandomForest)
Top 20 consensus features:

|   ConsensusRank | Feature                              | Feature_Type   | Species                         |   StabilityScore |   MeanAbsSHAP_Mean |   ConsensusScore |
|----------------:|:-------------------------------------|:---------------|:--------------------------------|-----------------:|-------------------:|-----------------:|
|               1 | PA_Parvimonas micra                  | PA             | Parvimonas micra                |              1   |         0.035732   |         1        |
|               2 | rCLR_Ruthenibacterium lactatiformans | rCLR           | Ruthenibacterium lactatiformans |              1   |         0.0155982  |         0.802786 |
|               3 | rCLR_Parvimonas micra                | rCLR           | Parvimonas micra                |              1   |         0.0139637  |         0.786776 |
|               4 | rCLR_Anaerostipes hadrus             | rCLR           | Anaerostipes hadrus             |              1   |         0.010678   |         0.754592 |
|               5 | rCLR_Roseburia intestinalis          | rCLR           | Roseburia intestinalis          |              1   |         0.00924486 |         0.740555 |
|               6 | rCLR_[Ruminococcus] torques          | rCLR           | [Ruminococcus] torques          |              1   |         0.00728966 |         0.721403 |
|               7 | rCLR_Lachnospira pectinoschiza       | rCLR           | Lachnospira pectinoschiza       |              0.9 |         0.00932981 |         0.676387 |
|               8 | rCLR_Eubacterium sp. CAG:251         | rCLR           | Eubacterium sp. CAG:251         |              0.9 |         0.00929491 |         0.676045 |
|               9 | rCLR_Bacteroides caccae              | rCLR           | Bacteroides caccae              |              0.9 |         0.00803233 |         0.663678 |
|              10 | rCLR_Dorea formicigenerans           | rCLR           | Dorea formicigenerans           |              0.9 |         0.00749121 |         0.658377 |
|              11 | rCLR_Parabacteroides distasonis      | rCLR           | Parabacteroides distasonis      |              0.9 |         0.00732896 |         0.656788 |
|              12 | rCLR_Lachnospira eligens             | rCLR           | Lachnospira eligens             |              0.9 |         0.00683686 |         0.651968 |
|              13 | rCLR_Eubacterium ventriosum          | rCLR           | Eubacterium ventriosum          |              0.8 |         0.00769537 |         0.595377 |
|              14 | rCLR_Clostridium sp. CAG:58          | rCLR           | Clostridium sp. CAG:58          |              0.8 |         0.00682689 |         0.58687  |
|              15 | rCLR_Bacteroides thetaiotaomicron    | rCLR           | Bacteroides thetaiotaomicron    |              0.7 |         0.00574987 |         0.511321 |
|              16 | rCLR_Odoribacter splanchnicus        | rCLR           | Odoribacter splanchnicus        |              0.6 |         0.00671739 |         0.455798 |
|              17 | rCLR_Butyricimonas virosa            | rCLR           | Butyricimonas virosa            |              0.6 |         0.00602039 |         0.448971 |
|              18 | rCLR_Fusicatenibacter saccharivorans | rCLR           | Fusicatenibacter saccharivorans |              0.6 |         0.0053088  |         0.442    |
|              19 | rCLR_Alistipes finegoldii            | rCLR           | Alistipes finegoldii            |              0.5 |         0.00646132 |         0.38829  |
|              20 | Age                                  | covariate      | Age                             |              0.5 |         0.0056598  |         0.380439 |

Feature-type split among top features:

| Feature_Type   |   N_Features |   MeanAbsSHAP_Total |   MeanAbsSHAP_Mean |   Percent_TopN |   Percent_SHAP_TopN |
|:---------------|-------------:|--------------------:|-------------------:|---------------:|--------------------:|
| rCLR           |           18 |           0.149868  |         0.00832603 |             90 |            78.3584  |
| PA             |            1 |           0.035732  |         0.035732   |              5 |            18.6824  |
| covariate      |            1 |           0.0056598 |         0.0056598  |              5 |             2.95921 |

## OnlyMale (RandomForest)
Top 20 consensus features:

|   ConsensusRank | Feature                                | Feature_Type   | Species                           |   StabilityScore |   MeanAbsSHAP_Mean |   ConsensusScore |
|----------------:|:---------------------------------------|:---------------|:----------------------------------|-----------------:|-------------------:|-----------------:|
|               1 | PA_Parvimonas micra                    | PA             | Parvimonas micra                  |              1   |         0.0345328  |         1        |
|               2 | rCLR_Parvimonas micra                  | rCLR           | Parvimonas micra                  |              1   |         0.0163115  |         0.815322 |
|               3 | rCLR_Ruthenibacterium lactatiformans   | rCLR           | Ruthenibacterium lactatiformans   |              1   |         0.0115913  |         0.767481 |
|               4 | rCLR_Roseburia intestinalis            | rCLR           | Roseburia intestinalis            |              1   |         0.011019   |         0.76168  |
|               5 | rCLR_Odoribacter splanchnicus          | rCLR           | Odoribacter splanchnicus          |              1   |         0.0109757  |         0.761242 |
|               6 | rCLR_Fusicatenibacter saccharivorans   | rCLR           | Fusicatenibacter saccharivorans   |              1   |         0.0108193  |         0.759657 |
|               7 | rCLR_Lachnospira eligens               | rCLR           | Lachnospira eligens               |              1   |         0.0102601  |         0.753989 |
|               8 | rCLR_Akkermansia muciniphila           | rCLR           | Akkermansia muciniphila           |              1   |         0.00978623 |         0.749186 |
|               9 | rCLR_Roseburia faecis                  | rCLR           | Roseburia faecis                  |              0.9 |         0.011425   |         0.700796 |
|              10 | rCLR_Eubacterium ventriosum            | rCLR           | Eubacterium ventriosum            |              0.9 |         0.00870795 |         0.673258 |
|              11 | rCLR_Lachnospira pectinoschiza         | rCLR           | Lachnospira pectinoschiza         |              0.9 |         0.00796732 |         0.665751 |
|              12 | rCLR_Firmicutes bacterium CAG:145      | rCLR           | Firmicutes bacterium CAG:145      |              0.8 |         0.00726665 |         0.59365  |
|              13 | rCLR_Bacteroides caccae                | rCLR           | Bacteroides caccae                |              0.8 |         0.00685695 |         0.589497 |
|              14 | rCLR_Bifidobacterium pseudocatenulatum | rCLR           | Bifidobacterium pseudocatenulatum |              0.6 |         0.00591592 |         0.44996  |
|              15 | rCLR_Eubacterium sp. CAG:251           | rCLR           | Eubacterium sp. CAG:251           |              0.5 |         0.006405   |         0.389917 |
|              16 | Age                                    | covariate      | Age                               |              0.5 |         0.00525731 |         0.378284 |
|              17 | rCLR_Eubacterium sp. CAG:38            | rCLR           | Eubacterium sp. CAG:38            |              0.5 |         0.00481455 |         0.373797 |
|              18 | rCLR_Butyricimonas virosa              | rCLR           | Butyricimonas virosa              |              0.5 |         0.00476648 |         0.37331  |
|              19 | rCLR_Anaerobutyricum hallii            | rCLR           | Anaerobutyricum hallii            |              0.5 |         0.00468661 |         0.3725   |
|              20 | rCLR_Paraprevotella xylaniphila        | rCLR           | Paraprevotella xylaniphila        |              0.4 |         0.00614025 |         0.322233 |

Feature-type split among top features:

| Feature_Type   |   N_Features |   MeanAbsSHAP_Total |   MeanAbsSHAP_Mean |   Percent_TopN |   Percent_SHAP_TopN |
|:---------------|-------------:|--------------------:|-------------------:|---------------:|--------------------:|
| rCLR           |           18 |          0.155716   |         0.00865088 |             90 |            79.6476  |
| PA             |            1 |          0.0345328  |         0.0345328  |              5 |            17.6633  |
| covariate      |            1 |          0.00525731 |         0.00525731 |              5 |             2.68908 |

## OnlyFemale (RandomForest)
Top 20 consensus features:

|   ConsensusRank | Feature                              | Feature_Type   | Species                         |   StabilityScore |   MeanAbsSHAP_Mean |   ConsensusScore |
|----------------:|:-------------------------------------|:---------------|:--------------------------------|-----------------:|-------------------:|-----------------:|
|               1 | rCLR_Ruthenibacterium lactatiformans | rCLR           | Ruthenibacterium lactatiformans |              1   |         0.022793   |         1        |
|               2 | PA_Parvimonas micra                  | PA             | Parvimonas micra                |              1   |         0.0159702  |         0.895232 |
|               3 | rCLR_Anaerostipes hadrus             | rCLR           | Anaerostipes hadrus             |              1   |         0.0133056  |         0.854315 |
|               4 | rCLR_Bacteroides fragilis            | rCLR           | Bacteroides fragilis            |              1   |         0.0110121  |         0.819097 |
|               5 | rCLR_Dorea formicigenerans           | rCLR           | Dorea formicigenerans           |              1   |         0.0103075  |         0.808277 |
|               6 | rCLR_Bacteroides thetaiotaomicron    | rCLR           | Bacteroides thetaiotaomicron    |              1   |         0.00962825 |         0.797848 |
|               7 | rCLR_Anaerobutyricum hallii          | rCLR           | Anaerobutyricum hallii          |              0.9 |         0.00834001 |         0.713066 |
|               8 | rCLR_Parabacteroides distasonis      | rCLR           | Parabacteroides distasonis      |              0.9 |         0.00781643 |         0.705026 |
|               9 | rCLR_Flavonifractor plautii          | rCLR           | Flavonifractor plautii          |              0.9 |         0.00640788 |         0.683397 |
|              10 | rCLR_[Ruminococcus] torques          | rCLR           | [Ruminococcus] torques          |              0.8 |         0.00702607 |         0.62789  |
|              11 | rCLR_Parvimonas micra                | rCLR           | Parvimonas micra                |              0.8 |         0.00636443 |         0.61773  |
|              12 | rCLR_Coprococcus comes               | rCLR           | Coprococcus comes               |              0.7 |         0.00688644 |         0.560745 |
|              13 | rCLR_Alistipes finegoldii            | rCLR           | Alistipes finegoldii            |              0.7 |         0.00564845 |         0.541735 |
|              14 | rCLR_Adlercreutzia equolifaciens     | rCLR           | Adlercreutzia equolifaciens     |              0.7 |         0.00512446 |         0.533689 |
|              15 | rCLR_Bacteroides caccae              | rCLR           | Bacteroides caccae              |              0.6 |         0.00504912 |         0.467532 |
|              16 | rCLR_Eubacterium sp. CAG:251         | rCLR           | Eubacterium sp. CAG:251         |              0.6 |         0.0046836  |         0.46192  |
|              17 | Age                                  | covariate      | Age                             |              0.5 |         0.0049304  |         0.400709 |
|              18 | rCLR_Eubacterium sp. CAG:180         | rCLR           | Eubacterium sp. CAG:180         |              0.5 |         0.00396821 |         0.385934 |
|              19 | rCLR_Parabacteroides merdae          | rCLR           | Parabacteroides merdae          |              0.4 |         0.00462963 |         0.331091 |
|              20 | rCLR_Paraprevotella xylaniphila      | rCLR           | Paraprevotella xylaniphila      |              0.4 |         0.00461509 |         0.330868 |

Feature-type split among top features:

| Feature_Type   |   N_Features |   MeanAbsSHAP_Total |   MeanAbsSHAP_Mean |   Percent_TopN |   Percent_SHAP_TopN |
|:---------------|-------------:|--------------------:|-------------------:|---------------:|--------------------:|
| rCLR           |           18 |           0.143606  |         0.00797812 |             90 |            87.295   |
| PA             |            1 |           0.0159702 |         0.0159702  |              5 |             9.70793 |
| covariate      |            1 |           0.0049304 |         0.0049304  |              5 |             2.99708 |

## OnlyMale vs OnlyFemale overlap
|   Top_N |   N_Male_Top |   N_Female_Top |   N_Shared |   N_Male_Specific |   N_Female_Specific |   Jaccard_Index |   Overlap_Percent_of_Male |   Overlap_Percent_of_Female |
|--------:|-------------:|---------------:|-----------:|------------------:|--------------------:|----------------:|--------------------------:|----------------------------:|
|      20 |           20 |             20 |          8 |                12 |                  12 |            0.25 |                        40 |                          40 |

| Category        | Feature                                | Feature_Type   | Species                           |   Male_Rank |   Female_Rank |   Male_ConsensusScore |   Female_ConsensusScore |   Male_Stability |   Female_Stability |
|:----------------|:---------------------------------------|:---------------|:----------------------------------|------------:|--------------:|----------------------:|------------------------:|-----------------:|-------------------:|
| shared          | Age                                    | covariate      | Age                               |          16 |            17 |              0.378284 |                0.400709 |              0.5 |                0.5 |
| shared          | PA_Parvimonas micra                    | PA             | Parvimonas micra                  |           1 |             2 |              1        |                0.895232 |              1   |                1   |
| shared          | rCLR_Anaerobutyricum hallii            | rCLR           | Anaerobutyricum hallii            |          19 |             7 |              0.3725   |                0.713066 |              0.5 |                0.9 |
| shared          | rCLR_Bacteroides caccae                | rCLR           | Bacteroides caccae                |          13 |            15 |              0.589497 |                0.467532 |              0.8 |                0.6 |
| shared          | rCLR_Eubacterium sp. CAG:251           | rCLR           | Eubacterium sp. CAG:251           |          15 |            16 |              0.389917 |                0.46192  |              0.5 |                0.6 |
| shared          | rCLR_Paraprevotella xylaniphila        | rCLR           | Paraprevotella xylaniphila        |          20 |            20 |              0.322233 |                0.330868 |              0.4 |                0.4 |
| shared          | rCLR_Parvimonas micra                  | rCLR           | Parvimonas micra                  |           2 |            11 |              0.815322 |                0.61773  |              1   |                0.8 |
| shared          | rCLR_Ruthenibacterium lactatiformans   | rCLR           | Ruthenibacterium lactatiformans   |           3 |             1 |              0.767481 |                1        |              1   |                1   |
| male_specific   | rCLR_Akkermansia muciniphila           | rCLR           | Akkermansia muciniphila           |           8 |           nan |              0.749186 |              nan        |              1   |              nan   |
| male_specific   | rCLR_Bifidobacterium pseudocatenulatum | rCLR           | Bifidobacterium pseudocatenulatum |          14 |           nan |              0.44996  |              nan        |              0.6 |              nan   |
| male_specific   | rCLR_Butyricimonas virosa              | rCLR           | Butyricimonas virosa              |          18 |           nan |              0.37331  |              nan        |              0.5 |              nan   |
| male_specific   | rCLR_Eubacterium sp. CAG:38            | rCLR           | Eubacterium sp. CAG:38            |          17 |           nan |              0.373797 |              nan        |              0.5 |              nan   |
| male_specific   | rCLR_Eubacterium ventriosum            | rCLR           | Eubacterium ventriosum            |          10 |           nan |              0.673258 |              nan        |              0.9 |              nan   |
| male_specific   | rCLR_Firmicutes bacterium CAG:145      | rCLR           | Firmicutes bacterium CAG:145      |          12 |           nan |              0.59365  |              nan        |              0.8 |              nan   |
| male_specific   | rCLR_Fusicatenibacter saccharivorans   | rCLR           | Fusicatenibacter saccharivorans   |           6 |           nan |              0.759657 |              nan        |              1   |              nan   |
| male_specific   | rCLR_Lachnospira eligens               | rCLR           | Lachnospira eligens               |           7 |           nan |              0.753989 |              nan        |              1   |              nan   |
| male_specific   | rCLR_Lachnospira pectinoschiza         | rCLR           | Lachnospira pectinoschiza         |          11 |           nan |              0.665751 |              nan        |              0.9 |              nan   |
| male_specific   | rCLR_Odoribacter splanchnicus          | rCLR           | Odoribacter splanchnicus          |           5 |           nan |              0.761242 |              nan        |              1   |              nan   |
| male_specific   | rCLR_Roseburia faecis                  | rCLR           | Roseburia faecis                  |           9 |           nan |              0.700796 |              nan        |              0.9 |              nan   |
| male_specific   | rCLR_Roseburia intestinalis            | rCLR           | Roseburia intestinalis            |           4 |           nan |              0.76168  |              nan        |              1   |              nan   |
| female_specific | rCLR_Adlercreutzia equolifaciens       | rCLR           | Adlercreutzia equolifaciens       |         nan |            14 |            nan        |                0.533689 |            nan   |                0.7 |
| female_specific | rCLR_Alistipes finegoldii              | rCLR           | Alistipes finegoldii              |         nan |            13 |            nan        |                0.541735 |            nan   |                0.7 |
| female_specific | rCLR_Anaerostipes hadrus               | rCLR           | Anaerostipes hadrus               |         nan |             3 |            nan        |                0.854315 |            nan   |                1   |
| female_specific | rCLR_Bacteroides fragilis              | rCLR           | Bacteroides fragilis              |         nan |             4 |            nan        |                0.819097 |            nan   |                1   |
| female_specific | rCLR_Bacteroides thetaiotaomicron      | rCLR           | Bacteroides thetaiotaomicron      |         nan |             6 |            nan        |                0.797848 |            nan   |                1   |
| female_specific | rCLR_Coprococcus comes                 | rCLR           | Coprococcus comes                 |         nan |            12 |            nan        |                0.560745 |            nan   |                0.7 |
| female_specific | rCLR_Dorea formicigenerans             | rCLR           | Dorea formicigenerans             |         nan |             5 |            nan        |                0.808277 |            nan   |                1   |
| female_specific | rCLR_Eubacterium sp. CAG:180           | rCLR           | Eubacterium sp. CAG:180           |         nan |            18 |            nan        |                0.385934 |            nan   |                0.5 |
| female_specific | rCLR_Flavonifractor plautii            | rCLR           | Flavonifractor plautii            |         nan |             9 |            nan        |                0.683397 |            nan   |                0.9 |
| female_specific | rCLR_Parabacteroides distasonis        | rCLR           | Parabacteroides distasonis        |         nan |             8 |            nan        |                0.705026 |            nan   |                0.9 |
| female_specific | rCLR_Parabacteroides merdae            | rCLR           | Parabacteroides merdae            |         nan |            19 |            nan        |                0.331091 |            nan   |                0.4 |
| female_specific | rCLR_[Ruminococcus] torques            | rCLR           | [Ruminococcus] torques            |         nan |            10 |            nan        |                0.62789  |            nan   |                0.8 |
