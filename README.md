# Reproducibility code — sepsis-associated liver injury (SALI), MIMIC-IV v3.1

Reproducibility material for the manuscript **"Development and External Validation of
Interpretable Machine-Learning Models for 7-, 14- and 28-Day Mortality in
Sepsis-Associated Liver Injury"**.

Two things are here:

1. **`code/`** — the PostgreSQL script that rebuilds the MIMIC-IV development cohort
   (adult ICU admissions → Sepsis-3 → SALI → index-hospitalisation selection → outcomes →
   participant-level exclusion → 31 predictors) directly from the raw MIMIC-IV v3.1
   source tables.
2. **`analysis/`** — the seven Python scripts that carry out the modelling, validation,
   comparator-score and supplementary-table computations, as archived from the original
   analysis.

**No patient-level data is included, and none may be redistributed from here.** See
"What is deliberately not in this repository".

> **Read [`AUDIT.md`](AUDIT.md) before you trust a number produced by this code.** The
> scripts were reviewed line by line before release; four issues are documented there,
> including a broken mid-rank function that invalidates the DeLong P-values in
> `compute_step3.py`, three different definitions of "modified SOFA", and a personal
> filesystem path that was removed from `prep_fig_data.py`.

---

> ## ⚠️ Provenance — read before citing or reusing
>
> **The SQL is a RECONSTRUCTION, not the script that produced the originally submitted
> export.** The original extraction script was not retained. The SQL here was written
> during the third revision round, directly against the raw source tables, so that a
> reader can repeat the cohort-ascertainment checks described in Sections 2.1–2.3 of the
> manuscript.
>
> It faithfully implements the rules **as they are now stated in the manuscript**, and it
> reproduces the reported denominator exactly and the reported case counts to within
> 1–2%. It does **not** claim — and cannot claim — to be byte-identical to the originally
> executed code.
>
> The Python scripts in `analysis/` **are** the original analysis scripts, released with
> only the four annotations described in `AUDIT.md`; no numerical result was altered.

---

## Contents

| Path | What it is |
|---|---|
| `code/MIMIC4_SALI_cohort_extraction.sql` | Cohort-extraction script (9 sections, PostgreSQL). |
| `analysis/dual_cohort/reproduce_models.py` | Loads and harmonises the two cohorts, checks Table 1, rebuilds the final models with the Supplementary Table S4 hyperparameters, and computes modified SOFA / MELD-Na. |
| `analysis/dual_cohort/compute_step3.py` | External AUCs with bootstrap CIs, comparator scores, decision-curve and classification metrics, SHAP ranking stability. |
| `analysis/dual_cohort/compute_step4.py` | Corrected Sun–Xu DeLong test, fair-baseline logistic regression, 14-day recalibration, DCA threshold range. |
| `analysis/dual_cohort/make_suppl_tables.py` | Writes Supplementary Tables S5 / S6 / S7 as `.docx`. |
| `analysis/single_center/run_single_center_v3.py` | Single-centre pipeline: the three selected models with AUC, calibration, DCA and comparator tests. |
| `analysis/single_center/prep_fig_data.py` | Figure data: ROC / calibration / DCA coordinates plus the 100-pipeline AUC matrix. |
| `analysis/single_center/extract_shap_values.py` | Per-sample SHAP values (KernelExplainer for the naive-Bayes models). |
| `docs/01_index_hospitalisation_restriction.md` | Executable definition of the index-hospitalisation restriction (one record per patient; deaths must fall inside the index admission). |
| `docs/02_extraction_notes.md` | What the SQL reproduces, side by side with the published counts, and which reviewer comment each section answers. |
| `AUDIT.md` | Line-by-line code audit performed before release. |

## Requirements

For the SQL:

* **MIMIC-IV v3.1** (`mimiciv_hosp`, `mimiciv_icu` schemas), PostgreSQL.
* Access requires the **PhysioNet credentialed access agreement**.
* The script writes intermediate tables into a scratch schema `sali`, which it creates.

For the Python scripts:

```bash
pip install -r requirements.txt
```

The scripts expect to be run from their own directory with the input data at
`../../01_原始数据/` — that is, place `analysis/` one level below the directory holding
`01_原始数据/` (or symlink it). Nothing else needs editing. See "What is deliberately not
in this repository" for how to obtain the inputs.

## How to run

```bash
# 1. rebuild the development cohort from the raw MIMIC-IV tables
psql -d mimiciv -f code/MIMIC4_SALI_cohort_extraction.sql

# 2. run the analysis, in this order
cd analysis/dual_cohort   && python reproduce_models.py && python compute_step3.py \
                          && python compute_step4.py   && python make_suppl_tables.py
cd ../single_center       && python run_single_center_v3.py && python prep_fig_data.py \
                          && python extract_shap_values.py
```

SQL sections run in order and build on one another:

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

## What the SQL reproduces

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

4. **The archived development-data extract has a known defect.** In the Excel extracts
   used for the original analysis, `vent`, `vaso` and `crrt` are constant at zero, which
   contradicts the 71.2 % / 60.1 % / 14.5 % in Table 1. A model retrained on those files
   will not reproduce the published AUCs. Regenerate the inputs with
   `code/MIMIC4_SALI_cohort_extraction.sql` first.

## What is deliberately **not** in this repository

| Excluded | Why |
|---|---|
| `01_原始数据/MIMIC队列/*.xlsx` | Derived from MIMIC-IV, released under the PhysioNet Credentialed Health Data License 1.5.0. Redistribution is prohibited — obtain MIMIC-IV from PhysioNet after credentialing. |
| `01_原始数据/单中心队列/*.xlsx` | Patient-level data from Xiaogan Central Hospital; release requires institutional approval. |
| `03_运行结果/shap_values.json`, `fig_data.json` | Contain per-patient feature values, predicted probabilities and outcomes for the validation set. Patient-level. |
| Frozen model objects (`.pkl` / `.joblib`) | **No such file exists** — none was ever produced. The models are fully reconstructible from `analysis/dual_cohort/reproduce_models.py`, which states the estimators, hyperparameters and preprocessing explicitly. |

The manuscript's Data Availability Statement records that items not publishable here are
available from the corresponding author on reasonable request.

## Citation

If you use this code, please cite the manuscript. The code itself may be reused under the
MIT licence (see `LICENSE`), but note the provenance statement above: the SQL is a
reconstruction of a published rule set, not the originally executed script.

## 中文说明

本仓库含两部分：

1. **`code/`** —— PostgreSQL 队列提取脚本，从 MIMIC-IV v3.1 原始表重建开发队列：成人 ICU 住院 →
   Sepsis-3 → SALI → index 住院选择 → 结局 → 参与者层面排除 → 31 个预测变量。
2. **`analysis/`** —— 7 个 Python 分析脚本，覆盖双队列建模与复现、单中心分析流程，即原分析代码本体。

**不含任何患者层面数据，也不得从此处再分发。**

**重要**：请先读 [`AUDIT.md`](AUDIT.md)。发布前已逐行审核，记录 4 个问题，其中两个必须知道：
`compute_step3.py` 的 midrank 函数有错，使其 DeLong P 值不可信（正确实现在 `compute_step4.py`）；
"改良 SOFA" 在不同脚本里有三套互不兼容的定义（`AUDIT.md` B1）。`prep_fig_data.py` 中原先硬编码的
本机绝对路径已删除。审核未改动任何数值结果。

**溯源**：SQL 是**重建实现**，不是当初执行并产生投稿导出文件的那份脚本——原脚本未予保留。它按稿件现在
所写的规则实现，分母 94,458 完全一致，各病例计数复现到 1–2% 以内，但不声称与原始代码逐字节相同。
`analysis/` 下的 Python 脚本则是原分析脚本，仅加了 `AUDIT.md` 所述的 4 处说明性注释。

三点容易误读之处：

1. 第 6 节 SALI 计数覆盖**全部**成人 ICU 住院（约 4,307），而稿件报告的 4,241 是在 34,216 例
   Sepsis-3 病例**内部**计数的；二者口径不同，不是同一个量。与 4,241 可比的数值由第 9 节核验块输出。
2. 第 4 节三项治疗指标是**整个 ICU 住院期间**口径，与表 1（71.2% / 60.1% / 14.5%）一致；首 24 小时
   口径会明显更低，不是表 1 所用口径。
3. 归档的 Excel 提取文件中 `vent`/`vaso`/`crrt` 三列**全为 0**，与表 1 矛盾；直接用它复跑得不到稿件
   的 AUC。请先用 `code/` 下的 SQL 重新生成输入。

**冻结模型对象不存在**：从未生成过任何 `.pkl`/`.joblib` 文件。模型可由
`analysis/dual_cohort/reproduce_models.py` 完整重建（估计器、超参数、预处理均在脚本中显式给出）。

运行前请先跑脚本第 0.3 节的两条预检查询，确认你所用 MIMIC-IV 版本的 itemid 字典一致。
