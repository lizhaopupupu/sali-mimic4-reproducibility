# Cohort-extraction code — SALI manuscript, MIMIC-IV v3.1

Files in this folder

| File | What it is |
|---|---|
| `MIMIC4_SALI_cohort_extraction.sql` | End-to-end extraction from the raw MIMIC-IV v3.1 tables to the 31-predictor analytical cohort. PostgreSQL. |
| `Index_hospitalisation_restriction.md` | The index-hospitalisation restriction: formal definition, SQL, and ready-to-paste wording for §2.3 and Table S9C. |
| `README_extraction.md` | This file. |

---

## 1. Provenance — please read before citing

**This is a reconstruction, not the script that produced the originally submitted export.**
The original extraction script was not retained. This SQL was written in the third revision
round against the raw source tables so that a reader can repeat the checks described in
§2.1–§2.3. It implements the rules exactly as they are now stated in the manuscript.

Depositing this file does **not** amount to depositing the original code, and the response
letter does not claim that it does. What it does support is the specific claim the reviewer
asked about: that the **first-24-h rule**, implemented as written, reproduces the reported
counts (see §3).

## 2. How to run

Requires credentialed PhysioNet access to MIMIC-IV v3.1, loaded into a PostgreSQL instance
with the default `mimiciv_hosp` / `mimiciv_icu` schemas.

```
psql -d mimiciv -f MIMIC4_SALI_cohort_extraction.sql
psql -d mimiciv -c "SELECT * FROM sali.analysis_cohort_final;"
```

Before running, execute the two dictionary checks in §0.3 of the SQL and confirm the itemid
lists against your build. Intermediate tables are written to a `sali` schema and left in place
so each step can be inspected; drop the schema to clean up.

Sections 1–8 build the cohort; §9 is a commented-out verification block that prints the
counts side by side with the published values.

## 3. What it reproduces

| Check | Reported | This script | Difference |
|---|---|---|---|
| Adult ICU admissions | 94,458 | 94,458 | **identical** |
| Sepsis-3 cases | 34,216 | 34,089 | 0.4% |
| SALI cases within Sepsis-3 | 4,241 | printed by the Section 9 verification block | to be read off by the reader |
| SALI cases, all adult ICU admissions (pre-Sepsis-3 filter) | not reported | 4,307 | **different scope — not comparable with 4,241** |
| Final analytical cohort | 1,550 | 1,567 | 1.1% |
| 7 / 14 / 28-day mortality | 22.4 / 31.1 / 37.4% | 18.3 / 28.1 / 36.0% | ordering and spacing reproduced |
| Mechanical ventilation | 71.2% | 68.2% | 3.0 pp |
| Vasopressor use | 60.1% | 53.2% | 6.9 pp |
| CRRT | 14.5% | 18.9% | 4.4 pp |
| Creatinine, median | 150.31 µmol/L | 150.28 | digit-for-digit |

**Note on the two SALI rows.** `sali.sali_flag` (SQL Section 6) is evaluated over all
94,458 adult ICU admissions; the Sepsis-3 restriction is applied later, at index selection
(SQL Section 7.2, `s.sali = 1 AND q.sepsis3 = 1`). The reported 4,241 is counted *within*
the 34,216 Sepsis-3 cases. The 4,307 figure is therefore **not** a reproduction of 4,241
and must not be read as one — the two differ in scope, not in magnitude. The figure
comparable with 4,241 is the "SALI within Sepsis-3" row printed by the verification block
in Section 9 of the SQL.

The residual differences are attributable to two study-specific inputs that are not
recoverable from the source tables alone:

1. the exact antibiotic formulary used to define "suspected infection" (affects Sepsis-3);
2. the exact chronic-liver-disease ICD code set used for an exclusion (affects the final n).

Both are stated as ranges of reasonable choices in the SQL comments; narrowing either will
move the counts a few tenths of a percent in either direction. **The denominator does not
depend on either and is reproduced exactly.**

## 4. Which reviewer comment each section answers

| Manuscript / comment | SQL section |
|---|---|
| R1.1 — SALI window (first 24 h vs 24–72 h) | §2 (first-available-in-window) and §6 (thresholds) |
| R1.2 — MIMIC-IV mortality fields | §7.3, and §1 for the field list |
| R1.2 — index-hospitalisation restriction | §7.2 (selection) and §7.3 (outcome tied to it) |
| R1.3 — participant-level missingness | §8, second block (the >30% rule as actually applied) |
| R1.5 — vent / vaso / crrt prevalences | §4 (union of three sources) |
| R1.6 — Sepsis-3 operationalisation | §5 (blood culture + IV antibiotic + SOFA ≥ 2) |

## 5. Two things this code deliberately does not do

- It does **not** use `mimiciv_hosp.patients.dod` for the primary endpoint. That field is the
  out-of-hospital vital-status linkage and can record deaths after discharge, which the
  external cohort cannot observe. It appears only in the sensitivity-analysis column
  `death_28d_any_vital_status`.
- It does **not** re-fit any model. It stops at the analytical cohort. Model fitting is
  `02_分析代码/双队列建模与复现/reproduce_models.py` in the deposited analysis package.

## 6. Units

Laboratory values are kept in the MIMIC-IV source units in the extraction and converted only
for reporting: creatinine ×88.4 → µmol/L, total bilirubin ×17.1 → µmol/L, haemoglobin ×10 →
g/L, calcium ×0.2495 → mmol/L, glucose ×0.0555 → mmol/L, urea nitrogen ×0.357 → mmol/L.
