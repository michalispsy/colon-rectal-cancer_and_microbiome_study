# ΒΗΜΑ 1.3 — Έλεγχος Κατανομών Μεταβλητών ανά Μελέτη

## 1. Sex Distribution ανά Study

| Study | Female (n) | Female (%) | Male (n) | Male (%) | Flag |
|---|---|---|---|---|---|
| FengQ_2015 | 67 | 43.5% | 87 | 56.5% | ✅ OK |
| GuptaA_2019 | 30 | 50.0% | 30 | 50.0% | ✅ OK |
| HanniganGD_2017 | 35 | 43.2% | 46 | 56.8% | ✅ OK |
| ThomasAM_2018a | 28 | 35.0% | 52 | 65.0% | ✅ OK |
| ThomasAM_2018b | 21 | 35.0% | 39 | 65.0% | ✅ OK |
| ThomasAM_2019_c | 35 | 43.8% | 45 | 56.2% | ✅ OK |
| VogtmannE_2016 | 30 | 28.8% | 74 | 71.2% | ⚠️ Borderline |
| WirbelJ_2018 | 52 | 41.6% | 73 | 58.4% | ✅ OK |
| YachidaS_2019 | 231 | 40.1% | 345 | 59.9% | ✅ OK |
| YuJ_2015 | 47 | 36.7% | 81 | 63.3% | ✅ OK |
| ZellerG_2014 | 69 | 44.2% | 87 | 55.8% | ✅ OK |

> [!NOTE]
> Καμία μελέτη δεν εμφανίζει extreme imbalance (< 25%). Η `VogtmannE_2016` (71.2% Male) είναι η πιο ανισόρροπη αλλά παραμένει εντός αποδεκτών ορίων.

---

## 2. Age Distribution ανά Study

| Study | Mean | Median | Min | Max | Missing | Flag |
|---|---|---|---|---|---|---|
| FengQ_2015 | 66.9 | 68.0 | 43 | 86 | 0 | ✅ OK |
| GuptaA_2019 | 50.6 | 57.0 | 22 | 75 | 0 | ✅ OK |
| HanniganGD_2017 | 58.6 | 59.0 | 35 | 88 | 0 | ✅ OK |
| ThomasAM_2018a | 67.5 | 67.0 | 49 | 84 | 0 | ✅ OK |
| ThomasAM_2018b | 58.2 | 59.0 | 38 | 70 | 1 | ✅ OK |
| ThomasAM_2019_c | 61.1 | 64.0 | 28 | 78 | 0 | ✅ OK |
| VogtmannE_2016 | 61.5 | 63.0 | 31 | 89 | 0 | ✅ OK |
| WirbelJ_2018 | 59.6 | 60.0 | 28 | 90 | 0 | ✅ OK |
| YachidaS_2019 | 61.9 | 64.0 | 21 | 79 | 0 | ✅ OK |
| YuJ_2015 | 64.2 | 64.0 | 34 | 89 | 0 | ✅ OK |
| ZellerG_2014 | 63.3 | 63.5 | 25 | 89 | 0 | ✅ OK |

> [!NOTE]
> Όλες οι μελέτες αφορούν ενήλικες (μέσος όρος 50–68 ετών). Η `GuptaA_2019` (Ινδία) είναι η πιο "νεαρή" κοόρτη (mean=50.6), ενώ η `ThomasAM_2018a` (Ιταλία) η πιο "ηλικιωμένη" (mean=67.5). Μόνο 1 τιμή ηλικίας λείπει σε όλο το dataset.

---

## 3. Country / Ethnicity ανά Study

| Study | Χώρα (ISO) | Ήπειρος |
|---|---|---|
| FengQ_2015 | 🇦🇹 Austria (AUT) | Ευρώπη |
| GuptaA_2019 | 🇮🇳 India (IND) | Ασία |
| HanniganGD_2017 | 🇺🇸 USA (87%) + 🇨🇦 Canada (13%) | Β. Αμερική |
| ThomasAM_2018a | 🇮🇹 Italy (ITA) | Ευρώπη |
| ThomasAM_2018b | 🇮🇹 Italy (ITA) | Ευρώπη |
| ThomasAM_2019_c | 🇯🇵 Japan (JPN) | Ασία |
| VogtmannE_2016 | 🇺🇸 USA | Β. Αμερική |
| WirbelJ_2018 | 🇩🇪 Germany (DEU) | Ευρώπη |
| YachidaS_2019 | 🇯🇵 Japan (JPN) | Ασία |
| YuJ_2015 | 🇨🇳 China (CHN) | Ασία |
| ZellerG_2014 | 🇫🇷 France (FRA) | Ευρώπη |

> [!NOTE]
> 9 χώρες, 3 ήπειροι. Κάθε μελέτη αντιστοιχεί σε 1 γεωγραφική περιοχή (εκτός HanniganGD_2017 που είναι μικτή USA/CAN). Η Ιαπωνία κυριαρχεί με 656 δείγματα (41% του dataset) από 2 μελέτες (YachidaS_2019 + ThomasAM_2019_c).

---

## 4. Συνοπτικό Πόρισμα

> [!IMPORTANT]
> **Δεν παρατηρείται κανένα extreme imbalance** σε φύλο, ηλικία ή γεωγραφική προέλευση ικανό να απαιτήσει αφαίρεση μελέτης. Το dataset θεωρείται κατάλληλο για downstream ανάλυση (Blocked Wilcoxon, ML models, LOSO cross-validation).
