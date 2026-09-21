# Code audit — SALI reproducibility repository

Audit performed on 2026-09-20 against the analysis code archived in
`数据与代码_SALI_2026-09-09/02_分析代码/`, immediately before public release.
Seven Python scripts and one SQL extraction script were reviewed line by line.

Severity key: **[A]** must fix before release · **[B]** must disclose · **[C]** cosmetic

---

## 1. Findings

### A1 — Hard-coded absolute path to a personal directory **[A] — fixed**

`single_center/prep_fig_data.py` contained

```python
FIG = r'C:/Users/intelligence/Documents/kimi/workspace/frontiers_revision/figures'
os.makedirs(FIG, exist_ok=True)
```

This is a personal filesystem path that would have been published verbatim, and the
script would silently create that directory on anyone else's machine. **Removed.** The
script writes all output to `../../03_运行结果/` and now creates no directory outside the
repository. No other behaviour changed.

### A2 — Incorrect mid-rank function makes the DeLong P-values in `compute_step3.py` unreliable **[A] — disclosed, superseded**

`dual_cohort/compute_step3.py::midrank()` returns

```python
s - (cnt[inv] - 1) / 2.0      # wrong
```

where `s` is the sum of ranks inside a tie group and `cnt` the group size. The correct
mid-rank is `s / cnt`. Consequences:

* The **AUCs printed by the script are unaffected** — they come from
  `sklearn.metrics.roc_auc_score`.
* The **DeLong P-values** produced through this helper are **not trustworthy**.

`compute_step4.py` contains a correct, standard Sun–Xu fast DeLong implementation
(`compute_midrank` / `fastDeLong` / `delong_p`) and supersedes the step-3 version. A
warning block has been added to the top of `compute_step3.py`.

### B1 — Three mutually incompatible definitions of "modified SOFA" across the codebase **[B] — disclosed**

| Script | Cardiovascular rule | GCS | Maximum |
|---|---|---|---|
| `dual_cohort/reproduce_models.py` | vasopressor → **4** pts, else MAP < 70 → 1 pt | yes | standard SOFA |
| `dual_cohort/compute_step3.py`, `compute_step4.py` | `vaso == 1 or MAP < 70` → **1** pt | yes | non-standard |
| `single_center/run_single_center_v3.py` | `vaso == 1 or MAP < 70` → **2** pts | **no** | non-standard |

The step-3/step-4 variant is annotated in the source as a *"best-matching variant"*,
i.e. it was chosen from among candidates because it came closest to the values already
reported in the manuscript. **This must be stated in the paper / cover letter if the
comparator score is challenged**, because a reviewer running `reproduce_models.py` will
get different SOFA AUCs than the ones in the manuscript. The three definitions are not
interchangeable. Warning comments have been added at each definition.

**Decision (2026-09-20): the definitions were deliberately NOT harmonised.** They were
quantified instead — see [`docs/03_comparator_score_definitions.md`](docs/03_comparator_score_definitions.md),
every number in which is reproducible with `analysis/dual_cohort/score_definition_impact.py`.
In short: on the external cohort the 1-point variant gives SOFA AUCs of 0.692 / 0.654 /
0.668 against 0.750 / 0.726 / 0.746 for the standard 4-point rule, and refitting the
models reproduces the manuscript's DeLong P-values to three decimal places under the
1-point variant (0.086 vs 0.088 at 14 days; 0.012 vs 0.011 at 28 days) but gives
P = 0.81 and P = 0.67 under the standard rule. Harmonising onto standard SOFA would
therefore remove the 28-day superiority claim. The score used for the published
comparison must be spelled out in the Methods.

One reason the standard rule cannot simply be applied: `vaso` is a **binary,
non-dose-graded** indicator, whereas standard SOFA awards its 4 vasopressor points only
above a dose threshold. That is a legitimate reason not to award 4 points — but it is a
reason that has to be written down, not one that can be left implicit in the code.

### B2 — MELD-Na is computed with two different sodium references **[B] — disclosed**

* `compute_step3.py` / `compute_step4.py`: `1.32 * (137 - Na) - 0.033 * MELD * (137 - Na)`, Na clamped to 125–137 — the 2016 OPTN formula.
* `single_center/run_single_center_v3.py`: `(135 - Na)`, Na clamped to 125–140 — the older MELD-Na form.

Both are defensible historical variants, but they give different scores — in this cohort
they differ for **all 229** patients (mean 21.9 vs 18.1) — and the two manuscripts
therefore do not report the same MELD-Na. Should be declared.

Quantified alongside B1 in
[`docs/03_comparator_score_definitions.md`](docs/03_comparator_score_definitions.md): on
the external cohort the Na-137 form gives AUCs of 0.646 / 0.619 / 0.641 (manuscript
0.656 / 0.622 / 0.636) and the Na-135 form 0.601 / 0.590 / 0.616. **Also not harmonised**,
for the same reason as B1.

### B3 — `run_single_center_v3.py::delong_p()` is not a DeLong test **[B] — disclosed**

It is a paired bootstrap test of the AUC difference (5,000 resamples, seed 42). The
JSON keys it feeds (`p_vs_sofa`, `p_vs_meld`) are bootstrap P-values, and the output
printing them is labelled `P=`. A comment has been added in the source.

### B4 — Known defect in the archived development-data export **[B] — documented in README**

`vent`, `vaso` and `crrt` are constant at zero in the MIMIC-IV training/validation
extracts. This contradicts the 71.2 % / 60.1 % / 14.5 % reported in Table 1. Any model
retrained on the archived extract will **not** reproduce the manuscript's AUCs, and the
SHAP stability frequencies for those three features are not estimable. Run
`code/MIMIC4_SALI_cohort_extraction.sql` against MIMIC-IV v3.1 to regenerate inputs with
real treatment variables.

### C1 — Minor issues **[C] — left as-is, listed for completeness**

| Location | Issue |
|---|---|
| `reproduce_models.py` | `med_gcs = 1.0` is assigned and never used. |
| `reproduce_models.py` | GCS is imputed with a hard-coded `fillna(9)` while every other variable uses the training-set median. |
| `compute_step3.py` | Dead expression `v01_1 = (ty1 - a2 if False else (ty1 - a1))/m`. |
| `make_suppl_tables.py` | `TOP5` feature lists are hard-coded from the manuscript results rather than re-derived. |
| `requirements.txt` | Pins are extremely new (`numpy==2.4.6`, `matplotlib==3.11.0`); `>=` ranges are provided in the repository `requirements.txt` for portability. |

---

## 2. What is deliberately **not** in this repository

| Excluded | Reason |
|---|---|
| `01_原始数据/MIMIC队列/*.xlsx` | Derived from MIMIC-IV, which is released under the PhysioNet Credentialed Health Data License 1.5.0. Redistribution is prohibited; obtain MIMIC-IV directly from PhysioNet after credentialing. |
| `01_原始数据/单中心队列/*.xlsx` | Patient-level data from Xiaogan Central Hospital. Release requires institutional approval. |
| `03_运行结果/shap_values.json` | Contains per-patient raw feature values for the validation set. Patient-level; not publishable. |
| `03_运行结果/fig_data.json` | Contains per-patient predicted probabilities, outcomes and score values for the validation set. Patient-level; not publishable. |
| Any frozen model object (`.pkl`, `.joblib`) | **No such file exists anywhere in the source packages.** See below. |

## 3. Frozen model objects

The manuscript's Data Availability Statement offers a *"frozen 7-day LightGBM model
object"* on request. **No serialised model object existed in the archived packages** —
they contain only scripts and data, and no `.pkl` / `.joblib` / `.sav` / `.h5` file could
be found in any of them.

The objects in `models/` were therefore **generated from this repository**, by
`analysis/dual_cohort/freeze_models.py`, using the Supplementary Table S4 hyperparameters
and the preprocessing described in Methods Section 2.9 (`random_state=42` throughout):

| Object | Internal validation AUC | External AUC | Published (val / external) |
|---|---|---|---|
| `sali_7d_lightgbm.pkl` | 0.828 | 0.760 | 0.840 / 0.720 |
| `sali_14d_logreg.pkl` | 0.831 | 0.718 | 0.841 / 0.717 |
| `sali_28d_logreg.pkl` | 0.842 | 0.760 | 0.850 / 0.763 |

`preprocessing_parameters.json` holds the imputation medians and scaling means/SDs;
`MODEL_CARD.json` holds the hyperparameters, the feature order and the provenance.

**Read the provenance block in `MODEL_CARD.json` before using them.** They were refitted
from the archived development export, in which `vent`, `vaso` and `crrt` are constant at
zero (finding B4). They are therefore **not** the models that produced the published
performance estimates. The 14- and 28-day logistic regressions are unaffected by that
defect and reproduce the published external AUCs to within 0.003; the **7-day LightGBM
model is affected and differs by 0.040** (0.760 here against 0.720 published), because it
is the only one of the three that makes use of the three defective indicators.

The Data Availability Statement must not claim that the deposited objects produced the
reported predictions. Regenerating faithful objects requires running
`code/MIMIC4_SALI_cohort_extraction.sql` against MIMIC-IV v3.1 and re-running
`freeze_models.py`.

## 4. Files changed during the audit

| File | Change | Numerical effect |
|---|---|---|
| `analysis/single_center/prep_fig_data.py` | Removed hard-coded personal absolute path | none |
| `analysis/dual_cohort/compute_step3.py` | Added header warning (bad `midrank`) + SOFA-variant note | none |
| `analysis/dual_cohort/compute_step4.py` | Added SOFA-variant note | none |
| `analysis/single_center/run_single_center_v3.py` | Added note that `delong_p` is a bootstrap test | none |

No numerical output, hyperparameter, cohort definition or result was altered by this
audit.
