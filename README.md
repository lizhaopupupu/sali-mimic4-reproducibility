# Cohort-extraction code — sepsis-associated liver injury (SALI), MIMIC-IV v3.1

Reproducibility material for the manuscript **"Development and External Validation of
Interpretable Machine-Learning Models for 7-, 14- and 28-Day Mortality in
Sepsis-Associated Liver Injury"**.

This repository contains the **cohort-extraction code only**: the PostgreSQL script that
rebuilds the MIMIC-IV development cohort (adult ICU admissions → Sepsis-3 → SALI →
index-hospitalisation selection → outcomes → participant-level exclusion → 31 predictors)
directly from the raw MIMIC-IV v3.1 source tables.

---

> ## ⚠️ Provenance — read before citing or reusing
>
> **This code is a RECONSTRUCTION, not the script that produced the originally submitted
> export.** The original extraction script was not retained. The SQL here was written
> during the third revision round, directly against the raw source tables, so that a
> reader can repeat the cohort-ascertainment checks described in Sections 2.1–2.3 of the
> manuscript.
>
> It faithfully implements the rules **as they are now stated in the manuscript**, and it
> reproduces the reported denominator exactly and the reported case counts to within
> 1–2%. It does **not** claim — and cannot claim — to be byte-identical to the originally
> executed code.

---

## Contents

| Path | What it is |
|---|---|
| `code/MIMIC4_SALI_cohort_extraction.sql` | The extraction script (9 sections, PostgreSQL). |
| `docs/01_index_hospitalisation_restriction.md` | The executable definition of the index-hospitalisation restriction (one record per patient; deaths must fall inside the index admission). |
| `docs/02_extraction_notes.md` | What the script reproduces, side by side with the published counts, and which reviewer comment each section answers. |

## Requirements

* **MIMIC-IV v3.1** (`mimiciv_hosp`, `mimiciv_icu` schemas), PostgreSQL.
* Access requires the **PhysioNet credentialed access agreement**. This repository
  contains no patient-level data and no derived data.
* The script writes intermediate tables into a scratch schema `sali`, which it creates.

## How to run

```bash
psql -d mimiciv -f code/MIMIC4_SALI_cohort_extraction.sql
```

Sections run in order and build on one another:

1. Adult ICU cohort + index-hospitalisation definition
2. Laboratory values — **first** available measurement inside the first 24 h
3. Vital signs and GCS (first 24 h)
4. Treatment / organ-support indicators — **whole ICU stay** (see caveat 3 below)
5. Sepsis-3 (blood culture + IV antibiotic + SOFA ≥ 2)
6. SALI (total bilirubin > 2.0 mg/dL **and** INR > 1.5, first 24 h)
7. Chronic-liver-disease exclusion, index-stay selection, outcomes
8. Analytical cohort: 31 predictors + endpoints + participant-level exclusion
9. Verification block — run it and compare with the manuscript

Section 0.3 of the script gives two pre-flight queries to re-verify the `itemid`
dictionaries in your build. **Run them first**: MIMIC-IV itemids are stable within a
version but have changed between versions.

## What it reproduces

| Check | Manuscript | This script |
|---|---|---|
| Adult ICU admissions | 94,458 | 94,458 (**identical**) |
| Sepsis-3 cases | 34,216 | ~34,089 (0.4%) |
| SALI cases, all adult ICU admissions | not reported | ~4,307 (see caveat 1) |
| SALI cases within Sepsis-3 | 4,241 | printed by the verification block |
| Final analytical cohort | 1,550 | ~1,567 (1.1%) |
| 7 / 14 / 28-day mortality | not tabulated in the manuscript | 18.3 / 28.1 / 36.0 % |
| Mechanical ventilation | 71.2 % | 68.2 % |
| Vasopressor use | 60.1 % | 53.2 % |
| CRRT | 14.5 % | 18.9 % |
| Creatinine, median | 150.3 µmol/L | 150.28 µmol/L (digit-for-digit) |

## Caveats — please read these before interpreting any count

1. **The `sali_flag` count is not the published 4,241.** Section 6 evaluates SALI over
   **all** 94,458 adult ICU admissions; the Sepsis-3 restriction is applied later, at
   index selection (Section 7.2, `s.sali = 1 AND q.sepsis3 = 1`). The published 4,241 is
   counted **within** the 34,216 Sepsis-3 cases. The two are different quantities with
   different scopes. The figure comparable with 4,241 is printed by the verification
   block in Section 9.

2. **Residual differences have known causes.** The exact antibiotic formulary used to
   define "suspected infection" and the exact chronic-liver-disease ICD code set used for
   an exclusion are not recoverable from the source tables alone. Both are stated as
   ranges of reasonable choices in the SQL comments; narrowing either moves the counts a
   few tenths of a percent in either direction. **The denominator does not depend on
   either, and is reproduced exactly.**

3. **The three treatment indicators are whole-stay, not first-24-h.** Table 1 reports
   whole-stay prevalences (71.2 % / 60.1 % / 14.5 %). Section 4 therefore uses the whole
   ICU stay. A first-24-h variant gives materially lower prevalences and is *not* what
   Table 1 reports; the script documents how to produce it if you want to compare.

4. **No model fitting is included.** See "What is not in this repository" below.

## What is **not** in this repository

* the model-fitting / analysis code (pipeline search, cross-validation, SHAP);
* the frozen model objects (7-day LightGBM) and the 14-/28-day logistic-regression
  coefficient vectors — these are printed in Supplementary Table S8B of the manuscript;
* the preprocessing parameters (imputation medians, scaling means and SDs) — printed in
  Supplementary Table S8A;
* the external (Xiaogan) cohort, which cannot be shared under its data-governance terms.

The manuscript's Data Availability Statement records that these are available from the
corresponding author on reasonable request.

## Citation

If you use this code, please cite the manuscript. The code itself may be reused under the
MIT licence (see `LICENSE`), but note the provenance statement above: it is a
reconstruction of a published rule set, not the originally executed analysis.

## 中文说明

本仓库只含**队列提取代码**（PostgreSQL），用于从 MIMIC-IV v3.1 原始表重建开发队列：成人 ICU
住院 → Sepsis-3 → SALI → index 住院选择 → 结局 → 参与者层面排除 → 31 个预测变量。

**重要**：这是**重建实现**，不是当初执行并产生投稿导出文件的那份脚本——原脚本未予保留。它按稿件
现在所写的规则实现，分母 94,458 完全一致，各病例计数复现到 1–2% 以内，但不声称与原始代码逐字节
相同。

三点容易误读之处：

1. 第 6 节的 SALI 计数覆盖**全部**成人 ICU 住院（约 4,307），而稿件报告的 4,241 是在 34,216 例
   Sepsis-3 病例**内部**计数的；二者口径不同，不是同一个量。与 4,241 可比的数值由第 9 节核验块
   输出。
2. 第 4 节三项治疗指标是**整个 ICU 住院期间**口径，与表 1（71.2% / 60.1% / 14.5%）一致；
   首 24 小时口径会明显更低，不是表 1 所用口径。
3. 本仓库**不含**建模/分析代码、冻结模型对象与预处理参数；稿件数据可用性声明中说明这些可向通讯
   作者合理索取。

运行前请先跑脚本第 0.3 节的两条预检查询，确认你所用 MIMIC-IV 版本的 itemid 字典一致。
