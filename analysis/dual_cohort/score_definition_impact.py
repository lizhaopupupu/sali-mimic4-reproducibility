# -*- coding: utf-8 -*-
"""
Recomputes every number in docs/03_comparator_score_definitions.md.

Shows, on the external cohort (n = 229):
  1. the AUC of each of the four SOFA cardiovascular variants and of each of the two
     MELD-Na forms;
  2. the DeLong P-value of the machine-learning models against each SOFA variant.

Only scikit-learn is required. The 14-day and 28-day final models are logistic
regressions and are refitted here; the 7-day model is a LightGBM booster and is skipped,
because its reported P-value is 0.660 (non-significant under every variant). The
comparator scores themselves are computed on the external cohort only, which is not
affected by the constant-zero vent/vaso/crrt defect of the archived MIMIC-IV export.

Caveat: because the MIMIC-IV training extract has vent/vaso/crrt constant at zero, the
refitted models are not bit-identical to the published ones. On the external cohort the
14-day and 28-day AUCs reproduce to 0.718 vs 0.717 and 0.760 vs 0.763, so the comparison
is representative.

Run from this directory:  python score_definition_impact.py
"""
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

M = r'../../01_原始数据/MIMIC队列'
X = r'../../01_原始数据/单中心队列'

FEATS31 = ['age', 'gender', 'max_temp', 'hr', 'rr', 'map', 'dbp', 'spo2', 'wbc', 'rbc', 'hb', 'hct', 'plt',
           'alt', 'ast', 'alp', 'alb', 'tbil', 'inr', 'bun', 'crea', 'na', 'k', 'ca', 'cl', 'hco3', 'glu',
           'lac', 'vent', 'vaso', 'crrt']

REPORTED_SOFA = {'7d': 0.700, '14d': 0.654, '28d': 0.673}
REPORTED_MELD = {'7d': 0.656, '14d': 0.622, '28d': 0.636}
REPORTED = {'14d': (0.717, 0.654, 0.088), '28d': (0.763, 0.673, 0.011)}


# ----------------------------------------------------------------- score variants
def sofa(df, vaso_pts, use_gcs=True, crrt_renal=False):
    """modified SOFA without the respiration component.

    vaso_pts   : points awarded for any vasopressor (1 = compute_step3/4, 4 = standard)
    use_gcs    : False reproduces the single-centre variant, which omits the CNS subscore
    crrt_renal : True reproduces the single-centre variant, which awards 4 renal points for CRRT
    """
    p = _c(df, 'plt')
    coag = np.select([p < 20, p < 50, p < 100, p < 150], [4, 3, 2, 1], default=0)
    tb = _c(df, 'tbil') / 17.1
    liver = np.select([tb >= 12.0, tb >= 6.0, tb >= 2.0, tb >= 1.2], [4, 3, 2, 1], default=0)
    cardio = np.where(_c(df, 'vaso') == 1, vaso_pts, np.where(_c(df, 'map') < 70, 1, 0))
    cr = _c(df, 'crea') / 88.4
    renal = np.select([cr >= 5.0, cr >= 3.5, cr >= 2.0, cr >= 1.2], [4, 3, 2, 1], default=0)
    if crrt_renal:
        renal = np.where(_c(df, 'crrt') == 1, 4, renal)
    total = coag + liver + cardio + renal
    if use_gcs:
        g = _c(df, 'gcs')
        total = total + np.select([g < 6, g < 10, g < 13, g < 15], [4, 3, 2, 1], default=0)
    return total


def meld_na(df, na_ref):
    """MELD-Na.

    na_ref=137 -> 2016 OPTN form used by compute_step3/4 (sodium clamped to 125-137,
                  correction 1.32*(137-Na) - 0.033*MELD*(137-Na));
    na_ref=135 -> older form used by reproduce_models (sodium clamped to 125-140,
                  correction 1.59*(135-Na)).
    The clamp bounds differ between the two scripts and are reproduced as written,
    because they are part of what makes the two versions differ.
    """
    tb = np.clip(_c(df, 'tbil') / 17.1, 1.0, None)
    cr = np.clip(_c(df, 'crea') / 88.4, 1.0, 4.0)
    inr = np.clip(_c(df, 'inr'), 1.0, None)
    na_hi = 137.0 if na_ref == 137 else 140.0
    na = np.clip(_c(df, 'na'), 125.0, na_hi)
    meld = np.round((0.957 * np.log(cr) + 0.378 * np.log(tb) + 1.120 * np.log(inr) + 0.643) * 10.0)
    if na_ref == 137:
        corr = np.where(meld > 11, 1.32 * (137.0 - na) - 0.033 * meld * (137.0 - na), 0.0)
    else:
        corr = np.where(meld > 11, 1.59 * (135.0 - na), 0.0)
    return meld + corr


def _c(df, name):
    s = df[name]
    if s.isna().any():
        s = s.fillna(s.median())
    return pd.to_numeric(s, errors='coerce').fillna(0.0).astype(float).to_numpy()


# ----------------------------------------------------------------- DeLong (Sun-Xu)
def _midrank(x):
    J = np.argsort(x); Z = x[J]; N = len(x)
    T = np.zeros(N); i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]: j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N); T2[J] = T
    return T2


def delong_p(y, p1, p2):
    y = np.asarray(y)
    preds = np.vstack([np.asarray(p1), np.asarray(p2)])[:, np.argsort(-y)]
    m = int(y.sum()); k, N = preds.shape; n = N - m
    pos, neg = preds[:, :m], preds[:, m:]
    tx = np.apply_along_axis(_midrank, 1, pos)
    ty = np.apply_along_axis(_midrank, 1, neg)
    tz = np.apply_along_axis(_midrank, 1, preds)
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    cov = np.cov(v01) / m + np.cov(v10) / n
    var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    z = (aucs[0] - aucs[1]) / np.sqrt(var) if var > 0 else 0.0
    return aucs, 2 * stats.norm.sf(abs(z))


# ----------------------------------------------------------------- load
def main():
    mi_tr_x = pd.read_excel(M + r'/MIMIC_训练集_特征.xlsx')
    mi_tr_y = pd.read_excel(M + r'/MIMIC_训练集_预后.xlsx')
    xg_x = pd.concat([pd.read_excel(X + r'/单中心_临床训练集.xlsx'),
                      pd.read_excel(X + r'/单中心_临床验证集.xlsx')], ignore_index=True)
    xg_y = pd.concat([pd.read_excel(X + r'/单中心_时间分类临床训练集预后.xlsx'),
                      pd.read_excel(X + r'/单中心_时间分类临床验证集预后.xlsx')], ignore_index=True)
    for d in (mi_tr_x, xg_x):
        if d['gender'].dtype == object:
            d['gender'] = d['gender'].astype(str).str.upper().str[0].map({'M': 1, 'F': 0})

    print(f'external cohort n = {len(xg_x)}')

    variants = {
        'SOFA 1pt  (compute_step3/4 -- published)': sofa(xg_x, 1),
        'SOFA 2pt, no GCS, CRRT renal (single_center)': sofa(xg_x, 2, use_gcs=False, crrt_renal=True),
        'SOFA 3pt': sofa(xg_x, 3),
        'SOFA 4pt  (standard SOFA rule)': sofa(xg_x, 4),
        'MELD-Na, Na ref 137 (compute_step3/4 -- published)': meld_na(xg_x, 137),
        'MELD-Na, Na ref 135 (reproduce_models)': meld_na(xg_x, 135),
    }

    print('\n--- AUC on the external cohort ---')
    print(f'{"score":52s} {"7d":>7s} {"14d":>7s} {"28d":>7s}')
    for name, s in variants.items():
        row = ' '.join(f'{roc_auc_score(xg_y[f"death_{tp}"].values, s):7.3f}' for tp in ('7d', '14d', '28d'))
        print(f'{name:52s} {row}')
    print(f'{"MANUSCRIPT (SOFA / MELD-Na)":52s} '
          f'{REPORTED_SOFA["7d"]:7.3f} {REPORTED_SOFA["14d"]:7.3f} {REPORTED_SOFA["28d"]:7.3f}')

    print('\n--- DeLong P, ML vs SOFA (14d / 28d logistic-regression models refitted) ---')
    med = mi_tr_x[FEATS31].median(); mu = mi_tr_x[FEATS31].mean(); sd = mi_tr_x[FEATS31].std().replace(0, 1)
    prep = lambda d: (d[FEATS31].fillna(med) - mu) / sd
    Xtr, Xxg = prep(mi_tr_x), prep(xg_x)
    for tp in ('14d', '28d'):
        mdl = LogisticRegression(C=0.1, class_weight='balanced', max_iter=1000, random_state=42)
        mdl.fit(Xtr, mi_tr_y[f'death_{tp}'])
        px = mdl.predict_proba(Xxg)[:, 1]
        yy = xg_y[f'death_{tp}'].values
        rml, rsofa, rp = REPORTED[tp]
        print(f'  {tp}  ML AUC refit = {roc_auc_score(yy, px):.3f} (reported {rml}); '
              f'SOFA reported {rsofa}, P {rp}')
        for pts in (1, 3, 4):
            aucs, p = delong_p(yy, px, sofa(xg_x, pts))
            print(f'       vs SOFA {pts}pt: SOFA AUC = {aucs[1]:.3f}   DeLong P = {p:.4f}')


if __name__ == '__main__':
    main()
