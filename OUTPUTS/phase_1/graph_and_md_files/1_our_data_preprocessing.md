# Η Δική μας Επεξεργασία Δεδομένων (Data Preprocessing)

> **Αφετηρία:** 1,604 δείγματα × 934 species (relative abundances %)  
> **Τέλος:** 1,604 δείγματα × 214 features (Hurdle Encoding: 107 PA + 107 rCLR)

---

## Βήμα 1.1 — Επαλήθευση Ομοιογένειας Pipeline

Πριν ξεκινήσουμε οποιαδήποτε ανάλυση, επαληθεύσαμε μέσω R ότι **όλες οι 11 CRC μελέτες** που επιλέξαμε από το curatedMetagenomicData πέρασαν από ακριβώς την ίδια έκδοση του bioinformatics pipeline:

- **Κοινό MetaPhlAn3 runDate:** Όλες οι μελέτες φέρουν runDate = **2021-03-31**. Αυτό σημαίνει ότι χρησιμοποιήθηκαν οι ίδιες ακριβώς εκδόσεις εργαλείων και οι ίδιες βάσεις δεδομένων αναφοράς (marker gene databases) για την ταξινομική ανάθεση σε κάθε μελέτη.
- **Κοινό HUMAnN3 pipeline:** Η λειτουργική ανάλυση (pathways, gene families) έγινε επίσης με την ίδια έκδοση — αν και δεν αξιοποιούμε αυτά τα δεδομένα στην παρούσα μελέτη, η ομοιογένεια επιβεβαιώθηκε.

> Αυτό το βήμα εξασφαλίζει ότι τυχόν διαφορές μεταξύ μελετών **δεν** οφείλονται σε διαφορετικές εκδόσεις λογισμικού ή βάσεων δεδομένων — μόνο σε βιολογικούς ή πειραματικούς παράγοντες.

---

## Βήμα 1.2 — Αποκλεισμός Προβληματικών Μελετών (Confounding Checks)

Ελέγξαμε μέσω R αν κάποια από τις 11 μελέτες εισάγει μη-διορθώσιμο batch effect. Συγκεκριμένα:

- **Έλεγχος χωριστής αλληλούχισης cases/controls:** Σε ορισμένες μελέτες στη βιβλιογραφία (π.χ. GuptaA_2019, WirbelJ_2018 σε παλαιότερες αναλύσεις), cases και controls αλληλουχίστηκαν σε διαφορετικά batches, δημιουργώντας πλασματική σύγχυση μεταξύ νόσου και sequencing batch. **Στις δικές μας 11 μελέτες: cases και controls αλληλουχίστηκαν μαζί σε κάθε μελέτη — δεν υπάρχει confounding.**
- **Έλεγχος sequencing platform ανά condition:** Επαληθεύσαμε ότι δεν υπάρχει σύγχυση μεταξύ condition (CRC/Control) και πλατφόρμας αλληλούχισης εντός κάθε μελέτης.
- **Έλεγχος country ανά condition:** Κάθε μελέτη αντιστοιχεί σε μία χώρα (πλην HanniganGD_2017: USA/Canada), με cases και controls από τον ίδιο γεωγραφικό πληθυσμό.

> **Αποτέλεσμα:** Καμία μελέτη δεν αποκλείστηκε. Δεν βρέθηκε confounding μεταξύ condition και platform/country/batch σε καμία από τις 11 μελέτες.

---

## Βήμα 1.3 — Ποιοτικός Έλεγχος Μεταδεδομένων

Πριν αγγίξουμε τα μικροβιακά δεδομένα, ελέγξαμε τα metadata:

- **Δεν υπάρχουν κενά** στα κρίσιμα πεδία: Sample, Study, Condition, Gender (0 missing).
- Πεδία με λίγα κενά: Age (1 missing), BMI (17 missing) — αμελητέα.
- **3 τιμές Condition:** CRC (701), Control (694), Adenoma (209).
- **2 τιμές Gender:** Male (959, 59.8%), Female (645, 40.2%).
- Ελέγξαμε ότι **καμία μελέτη δεν έχει extreme imbalance** σε φύλο (<25%) — η χειρότερη είναι η VogtmannE_2016 (71.2% Male), εντός αποδεκτών ορίων.

---

## Βήμα 1.4 — Prevalence Filtering (1ο Φίλτρο)

**Πρόβλημα:** Από τα 934 species, πολλά ήταν εξαιρετικά σπάνια — εμφανίζονταν σε πολύ λίγα δείγματα. Αυτά εισάγουν στατιστικό θόρυβο λόγω πληθώρας μηδενικών τιμών.

**Τι κάναμε:** Κόψαμε κάθε species που δεν ανιχνεύθηκε σε τουλάχιστον **15% των δειγμάτων** (δηλαδή ≥241 από τα 1,604 δείγματα).

**Αποτέλεσμα:**

| | Πριν | Μετά | Αφαιρέθηκαν |
|---|---|---|---|
| **Species** | 934 | 190 | 744 (79.7%) |
| **Δείγματα** | 1,604 | 1,604 | 0 |

> Κρατήσαμε μόνο species με αρκετή παρουσία ώστε να μπορούμε να βγάλουμε στατιστικά αξιόπιστα συμπεράσματα. Το log αυτού του φίλτρου αποθηκεύτηκε στο αρχείο [`species_prevalence_stats.csv`](file:///home/zoetsou/ntua/iatriki/colon-rectal-cancer_and_microbiome_study/data/crc_study_final_data/species_level/species_prevalence_stats.csv) (934 εγγραφές — μία ανά species, με τo prevalence και αν κρατήθηκε).

---

## Βήμα 1.5 — Cross-Study Abundance Filtering (2ο Φίλτρο)

**Πρόβλημα:** Ακόμα και μετά το prevalence filter, κάποια species μπορεί να εμφανίζονταν κυρίως σε 1–2 μελέτες — αποτέλεσμα τεχνικών batch effects (π.χ. συγκεκριμένο DNA extraction kit "βγάζει" κάποιο μικρόβιο πιο εύκολα). Αν τα κρατούσαμε, θα μάθαιναν τα μοντέλα να αναγνωρίζουν τη **μελέτη**, όχι τη **νόσο**.

**Τι κάναμε:** Κρατήσαμε μόνο species με μέση σχετική αφθονία **≥0.1% (10⁻³)** σε τουλάχιστον **3 διαφορετικές μελέτες**.

**Αποτέλεσμα:**

| | Πριν | Μετά | Αφαιρέθηκαν |
|---|---|---|---|
| **Species** | 190 | 107 | 83 (43.7%) |
| **Δείγματα** | 1,604 | 1,604 | 0 |

> Αυτό το φίλτρο στοχεύει κατευθείαν στο batch effect: αν ένα μικρόβιο εμφανίζεται μόνο σε 1–2 εργαστήρια, δεν μπορούμε να ξέρουμε αν η παρουσία του είναι βιολογική ή τεχνητή. Το log αποθηκεύτηκε στο [`cross_study_filter_log.csv`](file:///home/zoetsou/ntua/iatriki/colon-rectal-cancer_and_microbiome_study/data/crc_study_final_data/species_level/cross_study_filter_log.csv) (190 εγγραφές — σε πόσες μελέτες ξεπερνά το threshold, αν κρατήθηκε, μέση αφθονία).

**Το filtered dataset αποθηκεύτηκε ως:**  
[`species_filtered_crossstudy.csv`](file:///home/zoetsou/ntua/iatriki/colon-rectal-cancer_and_microbiome_study/data/crc_study_final_data/species_level/species_filtered_crossstudy.csv) — 1,604 δείγματα × 108 στήλες (107 species + Sample ID). Οι τιμές εδώ είναι ακόμα **raw relative abundances (%)**.

---

## Βήμα 1.6 — Εξέταση Αραιότητας (Sparsity Check)

Μετά το φιλτράρισμα, ελέγξαμε πόσα μηδενικά παρέμειναν στον πίνακα:

| Μέτρο | Τιμή |
|---|---|
| Κελιά πίνακα | 171,628 |
| Μηδενικά (zeros) | 81,562 |
| **Sparsity** | **47.52%** |

> Σχεδόν τα μισά κελιά είναι μηδενικά. Αυτό σημαίνει ότι **δεν μπορούμε** να εφαρμόσουμε απλό CLR μετασχηματισμό χωρίς ειδική διαχείριση μηδενικών — γι' αυτό προχωρήσαμε σε Hurdle Encoding.

---

## Βήμα 1.7 — Hurdle Encoding (Τελικός Μετασχηματισμός)

**Πρόβλημα:** Τα δεδομένα relative abundance έχουν διπλή φύση:
1. Πολλά μηδενικά (παρουσία/απουσία ενός είδους)
2. Μεγάλη διακύμανση στα μη-μηδενικά (πόσο άφθονο είναι όταν υπάρχει)

Ένας απλός CLR μετασχηματισμός δεν μπορεί να χειριστεί αυτή τη διττή φύση — χρειάζεται pseudocount για τα μηδενικά, που εισάγει τεχνητές τιμές.

**Τι κάναμε:** Για κάθε ένα από τα 107 species, δημιουργήσαμε **2 ξεχωριστά features** ανά δείγμα:

1. **PA (Presence/Absence):**
   - `1` αν το species ανιχνεύθηκε στο δείγμα (abundance > 0)
   - `0` αν δεν ανιχνεύθηκε

2. **rCLR (robust Centered Log-Ratio):**
   - Αν το species υπάρχει: `rCLR = log(αφθονία / γεωμετρικός μέσος μη-μηδενικών species του δείγματος)`
   - Αν δεν υπάρχει: `rCLR = 0`

> **Σημαντικό:** Ο γεωμετρικός μέσος υπολογίστηκε ανά δείγμα, **μόνο** από τα μη-μηδενικά species αυτού του δείγματος — χωρίς pseudocount.

**Αποτέλεσμα:**

| | Πριν | Μετά |
|---|---|---|
| **Δείγματα** | 1,604 | 1,604 |
| **Features** | 107 (raw abundances) | **214** (107 PA + 107 rCLR) |

> Ο τελικός πίνακας αποθηκεύτηκε ως [`X_explore.csv`](file:///home/zoetsou/ntua/iatriki/colon-rectal-cancer_and_microbiome_study/data/crc_study_final_data/species_level/X_explore.csv) — 1,604 × 215 (214 features + στήλη Sample). Σκοπός του είναι η χρήση σε unsupervised μεθόδους (PCA, t-SNE, clustering) για κατανόηση της δομής των δεδομένων.

---

## Συνοπτικό Διάγραμμα

```
934 species (raw from curatedMetagenomicData)
    │
    ▼ Prevalence Filter ≥15%        ← κόψαμε 744 σπάνια species
190 species
    │
    ▼ Cross-Study Filter ≥0.1% σε ≥3 studies  ← κόψαμε 83 batch-specific species
107 species
    │
    ▼ Hurdle Encoding (PA + rCLR)   ← μετασχηματισμός, ΟΧΙ κόψιμο
214 features (107 PA + 107 rCLR)
```

> **Τα δείγματα (1,604) δεν κόπηκαν σε κανένα βήμα** — όλο το φιλτράρισμα αφορούσε μόνο species.

---

## Τελικά Αρχεία στον Φάκελο Data

| Αρχείο | Διαστάσεις | Ρόλος |
|---|---|---|
| [`metadata.csv`](file:///home/zoetsou/ntua/iatriki/colon-rectal-cancer_and_microbiome_study/data/crc_study_final_data/species_level/metadata.csv) | 1,604 × 6 | Κλινικά/δημογραφικά (Sample, Study, Condition, Gender, Age, BMI) |
| [`species_prevalence_stats.csv`](file:///home/zoetsou/ntua/iatriki/colon-rectal-cancer_and_microbiome_study/data/crc_study_final_data/species_level/species_prevalence_stats.csv) | 934 × 3 | Log 1ου φίλτρου |
| [`cross_study_filter_log.csv`](file:///home/zoetsou/ntua/iatriki/colon-rectal-cancer_and_microbiome_study/data/crc_study_final_data/species_level/cross_study_filter_log.csv) | 190 × 4 | Log 2ου φίλτρου |
| [`species_filtered_crossstudy.csv`](file:///home/zoetsou/ntua/iatriki/colon-rectal-cancer_and_microbiome_study/data/crc_study_final_data/species_level/species_filtered_crossstudy.csv) | 1,604 × 108 | Raw abundances μετά filtering (107 species + Sample) |
| [`X_explore.csv`](file:///home/zoetsou/ntua/iatriki/colon-rectal-cancer_and_microbiome_study/data/crc_study_final_data/species_level/X_explore.csv) | 1,604 × 215 | Τελικός πίνακας Hurdle (214 features + Sample) |
