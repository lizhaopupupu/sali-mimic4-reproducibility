# The index-hospitalisation restriction

Draft wording + implementation for Reviewer 1 comment R1.2. Three parts:
(a) the definition, (b) the SQL, (c) ready-to-paste prose for the manuscript and the letter.

---

## (a) Definition

> The index hospitalisation is the single hospital admission during which the patient first
> satisfied the SALI criteria. Each patient contributes exactly one record: among all ICU stays
> satisfying the eligibility criteria, the stay with the earliest ICU admission time is retained
> and all later stays and later admissions are discarded. Deaths are counted only when the
> in-hospital death timestamp falls within the same hospital admission as the index ICU stay
> and within 7, 14 or 28 days of index ICU admission. Patients discharged alive before the
> horizon elapsed are classified as survivors at that horizon.

Two consequences worth stating explicitly, because they are what the reviewer actually asked
about:

1. **No post-discharge death enters the endpoint.** Because the death must fall inside the
   index admission, a death recorded only in the out-of-hospital vital-status linkage
   (`patients.dod`) cannot create an event.
2. **The two cohorts target the same estimand.** The external (Xiaogan) cohort has no
   out-of-hospital linkage at all, so restricting MIMIC-IV to the index hospitalisation is
   precisely what makes the development and external endpoints comparable. This is the
   substantive reason the restriction exists, not a data-processing convenience.

---

## (b) SQL

Selection of the index stay (one row per patient):

```sql
SELECT DISTINCT ON (w.subject_id)
       w.subject_id, w.hadm_id, w.stay_id
FROM sali.icu_adult w
JOIN sali.sali_flag s ON s.stay_id = w.stay_id      -- SALI criteria met
JOIN sali.sepsis3   q ON q.stay_id = w.stay_id      -- Sepsis-3 met
LEFT JOIN sali.cld  c ON c.stay_id = w.stay_id      -- chronic liver disease
WHERE s.sali = 1
  AND q.sepsis3 = 1
  AND c.stay_id IS NULL
ORDER BY w.subject_id, w.intime;                    -- earliest eligible stay wins
```

Endpoint, tied to the index admission:

```sql
CASE WHEN w.hospital_expire_flag = 1
      AND w.deathtime IS NOT NULL
      AND w.deathtime >= w.intime
      AND w.deathtime <  w.intime + INTERVAL '28 days'
      AND w.deathtime <= w.dischtime      -- same admission, made explicit
     THEN 1 ELSE 0 END AS death_28d
```

`w.deathtime <= w.dischtime` is satisfied by construction once the index `hadm_id` has been
selected — `deathtime` and `dischtime` both live on that admission row. It is written out
anyway so that the restriction is visible to a reader rather than implicit in a join.

---

## (c) Ready-to-paste prose

### For Manuscript §2.3 (English)

> Deaths were restricted to the index hospitalisation, defined as the hospital admission
> during which the patient first met the SALI criteria. Where a patient had more than one
> eligible ICU stay, the stay with the earliest ICU admission time was retained and all later
> stays were discarded, so each patient contributed one record. An event at the 7-, 14- or
> 28-day horizon required that the in-hospital death timestamp fall within the index
> hospitalisation and within the corresponding horizon of index ICU admission; patients
> discharged alive before the horizon elapsed were classified as survivors. Deaths recorded
> only in the out-of-hospital vital-status linkage were therefore not counted as events.

### For Manuscript §2.3（中文）

> 死亡限定在 index 住院期间；index 住院指患者首次满足 SALI 标准的那一次住院。若患者有不止一次
> 符合条件的 ICU 住院，则保留 ICU 入院时间最早的那一次，其后各次均剔除，因此每位患者仅贡献一条
> 记录。7、14 或 28 天时间窗内的事件，要求院内死亡时间戳落在 index 住院期间之内，且在 index ICU
> 入院后相应的时间窗之内；在时间窗届满前存活出院者归为存活者。因此，仅在院外生存状态联动中记录
> 到的死亡不计为事件。

### For the response letter, R1.2（英文，追加在"index 住院限定"一句之后）

> Deaths were restricted to the index hospitalisation in the sense that (i) each patient
> contributed one record — the earliest ICU stay at which the eligibility criteria were met —
> and (ii) an event required the in-hospital death timestamp to fall inside that same hospital
> admission and within the horizon of index ICU admission. Because the external cohort has no
> out-of-hospital vital-status linkage, this restriction is what makes the two cohorts target
> the same estimand; it is implemented in Section 7.2–7.3 of the deposited extraction code.

### For the response letter, R1.2（中文）

> “限定在 index 住院期间”的具体含义是：（i）每位患者只贡献一条记录——即首次满足入组标准的那次
> ICU 住院；（ii）事件要求院内死亡时间戳落在同一次住院之内，且在 index ICU 入院后的相应时间窗内。
> 由于外部队列完全没有院外生存状态联动，这一限定正是两个队列能针对同一估计目标的原因；其实现见
> 随附提取代码的第 7.2–7.3 节。
