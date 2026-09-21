# -*- coding: utf-8 -*-
"""
Reproduction & completion pipeline for the Frontiers SALI manuscript.
Step 1: load data, harmonize, verify Table 1 stats and cohort sizes/mortality.
Step 2: rebuild final models with S4 hyperparameters; verify internal/external AUCs.
Step 3: compute modified SOFA / MELD-Na on external cohort; verify AUCs & DeLong P.
Outputs: results_reproduction.txt
"""
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings('ignore')

FEATS31 = ['age','gender','max_temp','hr','rr','map','dbp','spo2','wbc','rbc','hb','hct','plt',
           'alt','ast','alp','alb','tbil','inr','bun','crea','na','k','ca','cl','hco3','glu',
           'lac','vent','vaso','crrt']
OUT = open('results_reproduction.txt','w',encoding='utf-8')
def pr(*a):
    print(*a); OUT.write(' '.join(str(x) for x in a)+'\n')

# ---------- load ----------
mi_tr_x = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_训练集_特征.xlsx')
mi_tr_y = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_训练集_预后.xlsx')
mi_va_x = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_验证集_特征.xlsx')
mi_va_y = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_验证集_预后.xlsx')
xg_tr_x = pd.read_excel('../../01_原始数据/单中心队列/单中心_临床训练集.xlsx')
xg_tr_y = pd.read_excel('../../01_原始数据/单中心队列/单中心_时间分类临床训练集预后.xlsx')
xg_va_x = pd.read_excel('../../01_原始数据/单中心队列/单中心_临床验证集.xlsx')
xg_va_y = pd.read_excel('../../01_原始数据/单中心队列/单中心_时间分类临床验证集预后.xlsx')

def harmonize(df):
    df = df.copy()
    if df['gender'].dtype == object:
        df['gender'] = df['gender'].astype(str).str.upper().str[0].map({'M':1,'F':0})
    return df
mi_tr_x, mi_va_x = harmonize(mi_tr_x), harmonize(mi_va_x)
xg_x = harmonize(pd.concat([xg_tr_x, xg_va_x], ignore_index=True))
xg_y = pd.concat([xg_tr_y, xg_va_y], ignore_index=True)

pr('MIMIC train', mi_tr_x.shape, 'MIMIC val', mi_va_x.shape, 'XG external', xg_x.shape)
pr('Mortality rates % (7/14/28d):')
for nm, y in [('MIMIC train', mi_tr_y), ('MIMIC val', mi_va_y), ('XG external', xg_y)]:
    r = (y.mean()*100).round(1).tolist()
    pr(f'  {nm}: {r}  deaths={y.sum().tolist()}')

# ---------- Table 1 verification (mean±sd) ----------
pr('\n--- Table 1 check: MIMIC train vs XG external (mean ± sd) ---')
rows = [('age','Age'),('max_temp','Peak temp'),('hr','HR'),('rr','RR'),('map','MAP'),('dbp','DBP'),
        ('spo2','SpO2'),('wbc','WBC'),('rbc','RBC'),('hb','Hb'),('hct','Hct'),('plt','Plt'),
        ('alt','ALT'),('ast','AST'),('alp','ALP'),('alb','Alb'),('tbil','Tbili'),('inr','INR'),
        ('bun','BUN'),('crea','Crea'),('na','Na'),('k','K'),('ca','Ca'),('cl','Cl'),
        ('hco3','HCO3'),('glu','Glu'),('lac','Lac')]
for c, lab in rows:
    a, b = mi_tr_x[c].dropna(), xg_x[c].dropna()
    pr(f'  {lab:10s} {a.mean():8.1f} ± {a.std():8.1f}  |  {b.mean():7.1f} ± {b.std():7.1f}')
for c, lab in [('gender','Male%'),('vent','Vent%'),('vaso','Vaso%'),('crrt','CRRT%')]:
    pr(f'  {lab:10s} {mi_tr_x[c].mean()*100:8.1f}          |  {xg_x[c].mean()*100:7.1f}')

# ---------- preprocessing ----------
med = mi_tr_x[FEATS31].median()
mu, sd = mi_tr_x[FEATS31].mean(), mi_tr_x[FEATS31].std().replace(0,1)
def prep(df):
    X = df[FEATS31].copy()
    X = X.fillna(med)
    return (X - mu) / sd
Xtr, Xva, Xxg = prep(mi_tr_x), prep(mi_va_x), prep(xg_x)

# ---------- final models (S4 hyperparameters) ----------
from lightgbm import LGBMClassifier
from sklearn.linear_model import LogisticRegression
models = {
    '7d':  LGBMClassifier(boosting_type='gbdt', class_weight='balanced', colsample_bytree=1.0,
                          learning_rate=0.05, max_depth=3, min_child_samples=20, n_estimators=200,
                          num_leaves=31, objective='binary', random_state=42, verbose=-1),
    '14d': LogisticRegression(C=0.1, class_weight='balanced', max_iter=1000, random_state=42),
    '28d': LogisticRegression(C=0.1, class_weight='balanced', max_iter=1000, random_state=42),
}
pr('\n--- Final model AUCs (target: val 0.840/0.841/0.850; ext 0.720/0.717/0.763) ---')
preds = {}
for tp, m in models.items():
    ytr = mi_tr_y[f'death_{tp}']
    m.fit(Xtr, ytr)
    pv, px = m.predict_proba(Xva)[:,1], m.predict_proba(Xxg)[:,1]
    preds[tp] = (pv, px)
    av = roc_auc_score(mi_va_y[f'death_{tp}'], pv)
    ax = roc_auc_score(xg_y[f'death_{tp}'], px)
    pr(f'  {tp}: internal-val AUC={av:.3f}   external AUC={ax:.3f}')

# ---------- modified SOFA (no respiration) & MELD-Na on raw values ----------
def sofa_row(r):
    s = 0
    p = r['plt']
    s += 0 if p>=150 else 1 if p>=100 else 2 if p>=50 else 3 if p>=20 else 4
    tb = r['tbil']/17.1
    s += 0 if tb<1.2 else 1 if tb<2.0 else 2 if tb<6.0 else 3 if tb<12.0 else 4
    if r['vaso']==1: s += 4
    elif r['map']<70: s += 1
    g = r['gcs']
    s += 0 if g>=15 else 1 if g>=13 else 2 if g>=10 else 3 if g>=6 else 4
    cr = r['crea']/88.4
    s += 0 if cr<1.2 else 1 if cr<2.0 else 2 if cr<3.5 else 3 if cr<5.0 else 4
    return s

def meld_na(r):
    tb = max(r['tbil']/17.1, 1.0); cr = max(min(r['crea']/88.4, 4.0), 1.0); inr = max(r['inr'], 1.0)
    meld = 3.78*np.log(tb) + 11.2*np.log(inr) + 9.57*np.log(cr) + 6.43
    meld = round(meld)
    if meld > 11:
        na = min(max(r['na'],125),140)
        meld = meld + 1.59*(135-na)
    return meld

OUT.close()

# ================= STEP 2: scores on XG + calibration fidelity =================
import numpy as np
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.linear_model import LogisticRegression as LR2

OUT = open('results_reproduction.txt','a',encoding='utf-8')
med_gcs = 1.0
xg2 = xg_x.copy()
for c in ['plt','tbil','vaso','map','crea','na','inr']:
    xg2[c] = xg2[c].fillna(med[c])
xg2['gcs'] = xg2['gcs'].fillna(9)
xg2['SOFA_mod'] = xg2.apply(sofa_row, axis=1)
xg2['MELD_Na'] = xg2.apply(meld_na, axis=1)

def boot_ci(y, s, n=1000, seed=42):
    rng = np.random.default_rng(seed); y=np.asarray(y); s=np.asarray(s); aucs=[]
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(y[i]))<2: continue
        aucs.append(roc_auc_score(y[i], s[i]))
    return np.percentile(aucs,[2.5,97.5])

pr('\n--- Score AUCs + 95%CI on external (n=229) ---')
pr('Reported Table3: 7d SOFA 0.700(0.625-0.770) MELD 0.656(0.582-0.727); 14d 0.654(0.579-0.722)/0.622(0.549-0.693); 28d 0.673(0.605-0.736)/0.636(0.564-0.703)')
pr('Figures S3/S4: 7d MELD 0.638; 14d MELD 0.612')
for tp in ['7d','14d','28d']:
    y = xg_y[f'death_{tp}'].values
    for sc in ['SOFA_mod','MELD_Na']:
        a = roc_auc_score(y, xg2[sc]); ci = boot_ci(y, xg2[sc].values)
        pr(f'  {tp} {sc}: AUC={a:.3f} ({ci[0]:.3f}-{ci[1]:.3f})')

pr('\n--- Calibration fidelity check (reported: 7d slope 0.94 int 0.33; 14d slope 0.59; 28d slope 0.88 int 0.08) ---')
for tp,(pv,px) in preds.items():
    y = xg_y[f'death_{tp}'].values
    eps = 1e-15
    logit = np.log(np.clip(px,eps,1-eps)/(1-np.clip(px,eps,1-eps)))
    cal = LR2().fit(logit.reshape(-1,1), y)
    br = brier_score_loss(y, px)
    pr(f'  {tp}: ext slope={cal.coef_[0][0]:.2f} intercept={cal.intercept_[0]:.2f} Brier={br:.3f}')
OUT.close()
