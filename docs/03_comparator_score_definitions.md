# The comparator clinical scores: three SOFA definitions and two MELD-Na definitions

**Status: NOT unified, by decision of the authors.** The scripts in `analysis/` use
different definitions of "modified SOFA" and of MELD-Na. This document records exactly
what each definition is, which one produced the numbers published in the manuscript, and
what would change if a different one were adopted. Every claim below is reproducible with
`analysis/score_definition_impact.py`.

---

## 1. The three SOFA definitions in the code

All three omit the respiration component, which is why they are called "modified SOFA".
They differ in the cardiovascular subscore, in whether GCS is included, and in the renal
subscore.

| | `dual_cohort/reproduce_models.py` | `dual_cohort/compute_step3.py`, `compute_step4.py` | `single_center/run_single_center_v3.py`, `prep_fig_data.py` |
|---|---|---|---|
| **Coagulation** (platelets) | <20→4, <50→3, <100→2, <150→1 | same | same |
| **Liver** (bilirubin µmol/L ÷ 17.1) | ≥12→4, ≥6→3, ≥2→2, ≥1.2→1 | same | thresholds 205 / 102 / 34.2 / 20.5 µmol/L — **numerically identical** |
| **Cardiovascular** | vasopressor → **4**; else MAP<70 → 1 | `vaso==1 or MAP<70` → **1** | `vaso==1 or MAP<70` → **2** |
| **CNS** (GCS) | <6→4, <10→3, <13→2, <15→1 | same | **not included** |
| **Renal** (creatinine µmol/L ÷ 88.4) | ≥5→4, ≥3.5→3, ≥2→2, ≥1.2→1 | same | same, **plus CRRT → 4** |

Note that the standard SOFA cardiovascular subscore awards **4** points for vasopressors
— but at a *dose threshold* (dopamine > 15 µg/kg/min, or adrenaline / noradrenaline
> 0.1 µg/kg/min). In this dataset `vaso` is a **binary, non-dose-graded** indicator of any
vasopressor, so the standard rule cannot be applied as written: awarding 4 points assumes
every flagged patient was on high-dose vasopressors, which the data cannot confirm.

---

## 2. What each definition actually gives on the external cohort (n = 229)

The external cohort is unaffected by the constant-zero `vent`/`vaso`/`crrt` defect in the
archived MIMIC-IV export, so these numbers are sound.

| Definition | 7-day AUC | 14-day AUC | 28-day AUC | mean ± SD of score |
|---|---|---|---|---|
| **1 point** (`compute_step3/4`) | **0.692** | **0.654** | **0.668** | 7.75 ± 2.27 |
| 2 points, no GCS, CRRT renal (`single_center`) | 0.738 | 0.691 | 0.708 | 5.69 ± 2.33 |
| 3 points | 0.741 | 0.713 | 0.732 | — |
| **4 points** (standard SOFA rule) | 0.750 | 0.726 | 0.746 | 8.94 ± 3.14 |
| **Reported in the manuscript** | **0.700** | **0.654** | **0.673** | — |

(The single-centre variant's AUCs move by about 0.005 depending on whether missing
values are imputed with the training-set median, as `run_single_center_v3.py` does, or
with the median of the whole 229-patient frame, as this comparison script does. The
difference does not affect any conclusion.)

**The 1-point variant is the one that produced the published Table 3 SOFA AUCs.**
Refitting the 14-day and 28-day models and rerunning the DeLong test reproduces the
manuscript essentially digit-for-digit:

| | manuscript | recomputed with the 1-point variant |
|---|---|---|
| 14-day ML external AUC | 0.717 | 0.718 |
| 14-day SOFA AUC | 0.654 | 0.654 |
| 14-day DeLong P | 0.088 | **0.086** |
| 28-day ML external AUC | 0.763 | 0.760 |
| 28-day SOFA AUC | 0.673 | 0.668 |
| 28-day DeLong P | 0.011 | **0.012** |

---

## 3. Consequence of switching to the standard (4-point) SOFA

This is why the definitions were left as they are: changing them is **not** a cosmetic
edit. On the same external cohort and with the same refitted models:

| DeLong P, ML vs SOFA | 1 point (current) | 3 points | 4 points (standard) |
|---|---|---|---|
| 14 days | **0.086** | 0.876 | 0.808 |
| 28 days | **0.012** | 0.388 | 0.665 |

Under the 3- or 4-point rule the SOFA AUC rises to 0.71–0.75, i.e. to the level of the
machine-learning models themselves (0.717 at 14 days, 0.763 at 28 days), and the
28-day claim that the model outperforms SOFA (P = 0.011) would no longer be supported.
At 7 and 14 days the ordering would reverse.

**Therefore the cardiovascular subscore must be stated explicitly in the Methods.**
The score used for the published comparison awards **1 point for any vasopressor** and is
not the standard 4-point SOFA rule. This is defensible — the exposure variable is binary
and not dose-graded — but it has to be written down, because a reader who computes SOFA
the standard way will get different numbers and different P-values.

---

## 4. The two MELD-Na definitions in the code

| | `compute_step3.py`, `compute_step4.py` | `reproduce_models.py` | `run_single_center_v3.py`, `prep_fig_data.py` |
|---|---|---|---|
| Base MELD | 2016 OPTN: `10 × (0.957 ln Cr + 0.378 ln Tb + 1.120 ln INR + 0.643)` | same in closed form (`3.78 ln Tb + 11.2 ln INR + 9.57 ln Cr + 6.43`) — **algebraically identical** | same closed form |
| Sodium term | `1.32 × (137 − Na) − 0.033 × MELD × (137 − Na)` | `1.59 × (135 − Na)` | `1.32 × (135 − Na) − 0.033 × MELD × (135 − Na)` |
| Sodium clamp | 125–137 | 125–140 | 125–140 |

The sodium reference differs (137 vs 135) and so does the correction form. On the
external cohort:

| | 7-day | 14-day | 28-day |
|---|---|---|---|
| **MELD-Na, Na reference 137** (`compute_step3/4`) | **0.646** | **0.619** | **0.641** |
| MELD-Na, Na reference 135 (`reproduce_models`) | 0.601 | 0.590 | 0.616 |
| **Reported in the manuscript** | **0.656** | **0.622** | **0.636** |

Again, the definition in `compute_step3/4.py` — the 2016 OPTN form with sodium reference
137 — is the one behind the published numbers. The two definitions differ for **all 229**
patients.

---

## 5. What to do if you want to harmonise them

1. Pick the definition to keep. Preserving the published results means keeping the
   **1-point** SOFA and the **Na-137** MELD-Na, i.e. the versions in `compute_step3.py` /
   `compute_step4.py`.
2. Replace the `sofa` / `sofa_row` / `modified_sofa` and `meld_na` / `meldna` definitions
   in the other scripts with that one.
3. Recompute every SOFA and MELD-Na AUC and DeLong P in the manuscript and in the
   supplementary tables — adopting the standard 4-point SOFA changes Table 3 and removes
   the 28-day superiority claim, as shown above.
4. State the cardiovascular subscore explicitly in the Methods in any case.

`analysis/score_definition_impact.py` recomputes every number in this document.
