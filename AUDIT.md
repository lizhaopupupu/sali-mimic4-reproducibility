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

### B2 — MELD-Na is computed with two different sodium references **[B] — disclosed**

* `compute_step3.py` / `compute_step4.py`: `1.32 * (137 - Na) - 0.033 * MELD * (137 - Na)`, Na clamped to 125–137 — the 2016 OPTN formula.
* `single_center/run_single_center_v3.py`: `(135 - Na)`, Na clamped to 125–140 — the older MELD-Na form.

Both are defensible historical variants, but they give different scores, and the two
manuscripts therefore do not report the same MELD-Na. Should be declared.

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

## 3. Frozen model object

The manuscript's Data Availability Statement offers a *"frozen 7-day LightGBM model
object"* on request. **No serialised model object exists** — the archived packages
contain only scripts and data, and no `.pkl` / `.joblib` / `.sav` / `.h5` file could be
found in any of them.

* The model is fully reconstructible from this repository: the estimator and its
  hyperparameters are explicit in `dual_cohort/reproduce_models.py` (`LGBMClassifier`,
  `random_state=42`) and the preprocessing (training-set median imputation, training-set
  mean/SD standardisation) is in the same file.
* A model object regenerated from the **archived** extract would be trained with
  `vent`/`vaso`/`crrt` all zero, i.e. it would not be the model behind the published
  results. Generating one and presenting it as the manuscript's model would be misleading.
* Recommendation: keep the model object on request, and point the Data Availability
  Statement at this repository for the code and preprocessing parameters.

## 4. Files changed during the audit

| File | Change | Numerical effect |
|---|---|---|
| `analysis/single_center/prep_fig_data.py` | Removed hard-coded personal absolute path | none |
| `analysis/dual_cohort/compute_step3.py` | Added header warning (bad `midrank`) + SOFA-variant note | none |
| `analysis/dual_cohort/compute_step4.py` | Added SOFA-variant note | none |
| `analysis/single_center/run_single_center_v3.py` | Added note that `delong_p` is a bootstrap test | none |

No numerical output, hyperparameter, cohort definition or result was altered by this
audit.
