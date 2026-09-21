# -*- coding: utf-8 -*-
"""Step 4: corrected DeLong + refined placeholder computations."""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.linear_model import LogisticRegression
from lightgbm import LGBMClassifier
from scipy import stats

FEATS31 = ['age','gender','max_temp','hr','rr','map','dbp','spo2','wbc','rbc','hb','hct','plt',
           'alt','ast','alp','alb','tbil','inr','bun','crea','na','k','ca','cl','hco3','glu',
           'lac','vent','vaso','crrt']
OUT = open('results_step4.txt','w',encoding='utf-8')
def pr(*a):
    print(*a); OUT.write(' '.join(str(x) for x in a)+'\n')

mi_tr_x = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_训练集_特征.xlsx')
mi_tr_y = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_训练集_预后.xlsx')
xg_x = pd.concat([pd.read_excel('../../01_原始数据/单中心队列/单中心_临床训练集.xlsx'),
                  pd.read_excel('../../01_原始数据/单中心队列/单中心_临床验证集.xlsx')], ignore_index=True)
xg_y = pd.concat([pd.read_excel('../../01_原始数据/单中心队列/单中心_时间分类临床训练集预后.xlsx'),
                  pd.read_excel('../../01_原始数据/单中心队列/单中心_时间分类临床验证集预后.xlsx')], ignore_index=True)
med = mi_tr_x[FEATS31].median(); mu = mi_tr_x[FEATS31].mean(); sd = mi_tr_x[FEATS31].std().replace(0,1)
def prep(df, feats=FEATS31):
    return (df[feats].fillna(med[feats]) - mu[feats]) / sd[feats]
Xtr, Xxg = prep(mi_tr_x), prep(xg_x)

def mk(tp):
    if tp=='7d':
        return LGBMClassifier(boosting_type='gbdt', class_weight='balanced', colsample_bytree=1.0,
                              learning_rate=0.05, max_depth=3, min_child_samples=20, n_estimators=200,
                              num_leaves=31, objective='binary', random_state=42, verbose=-1)
    return LogisticRegression(C=0.1, class_weight='balanced', max_iter=1000, random_state=42)

# --- standard fast DeLong (Sun-Xu) ---
def compute_midrank(x):
    J = np.argsort(x); Z = x[J]; N = len(x)
    T = np.zeros(N); i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]: j += 1
        T[i:j] = 0.5*(i+j-1)+1
        i = j
    T2 = np.empty(N); T2[J] = T
    return T2

def fastDeLong(preds_sorted, m):
    k, N = preds_sorted.shape
    n = N - m
    pos = preds_sorted[:, :m]; neg = preds_sorted[:, m:]
    tx = np.apply_along_axis(compute_midrank, 1, pos)
    ty = np.apply_along_axis(compute_midrank, 1, neg)
    tz = np.apply_along_axis(compute_midrank, 1, preds_sorted)
    aucs = tz[:, :m].sum(axis=1)/m/n - (m+1.0)/2.0/n
    v01 = (tz[:, :m] - tx)/n
    v10 = 1.0 - (tz[:, m:] - ty)/m
    sx = np.cov(v01); sy = np.cov(v10)
    return aucs, sx/m + sy/n

def delong_p(y, p1, p2):
    y = np.asarray(y)
    order = np.argsort(-y)  # positives first
    preds = np.vstack([np.asarray(p1), np.asarray(p2)])[:, order]
    m = int(y.sum())
    aucs, cov = fastDeLong(preds, m)
    var = cov[0,0] + cov[1,1] - 2*cov[0,1]
    z = (aucs[0]-aucs[1])/np.sqrt(var) if var>0 else 0
    return aucs, 2*stats.norm.sf(abs(z))

# sanity check of DeLong: ML vs identical predictor should be P=1
y28 = xg_y['death_28d'].values
m28 = mk('28d'); m28.fit(Xtr, mi_tr_y['death_28d']); p28 = m28.predict_proba(Xxg)[:,1]
m14 = mk('14d'); m14.fit(Xtr, mi_tr_y['death_14d']); p14 = m14.predict_proba(Xxg)[:,1]
m7  = mk('7d');  m7.fit(Xtr, mi_tr_y['death_7d']);  p7  = m7.predict_proba(Xxg)[:,1]
aucs, p = delong_p(y28, p28, p28); pr(f'DeLong self-check (expect P=1): P={p:.4f}')

# scores on XG
# NOTE (repository release): this is the same non-standard SOFA variant as in
# `compute_step3.py` - vasopressor use contributes 1 point, not the standard 4.  It was
# selected because it reproduced the manuscript's reported values most closely.
# `reproduce_models.py` (standard SOFA, 4 points for vasopressors) and
# `../single_center/run_single_center_v3.py` (2 points, no GCS) differ from it.
# See AUDIT.md.
def sofa(r):
    s=0
    p=r['plt']; s += 0 if p>=150 else 1 if p>=100 else 2 if p>=50 else 3 if p>=20 else 4
    tb=r['tbil']/17.1; s += 0 if tb<1.2 else 1 if tb<2.0 else 2 if tb<6.0 else 3 if tb<12.0 else 4
    if r['vaso']==1 or r['map']<70: s += 1
    g=r['gcs']; s += 0 if g>=15 else 1 if g>=13 else 2 if g>=10 else 3 if g>=6 else 4
    cr=r['crea']/88.4; s += 0 if cr<1.2 else 1 if cr<2.0 else 2 if cr<3.5 else 3 if cr<5.0 else 4
    return s
def meldna(r):
    tb=max(r['tbil']/17.1,1.0); cr=max(min(r['crea']/88.4,4.0),1.0); inr=max(r['inr'],1.0)
    m=(0.957*np.log(cr)+0.378*np.log(tb)+1.120*np.log(inr)+0.643)*10; m=round(m)
    if m>11:
        na=min(max(r['na'],125),137)
        m = m + 1.32*(137-na) - 0.033*m*(137-na)
    return m
xgs = xg_x.copy()
for c in ['plt','tbil','vaso','map','crea','na','inr','gcs']:
    xgs[c] = xgs[c].fillna(xgs[c].median())
xgs['SOFA_mod']=xgs.apply(sofa,axis=1); xgs['MELD_Na']=xgs.apply(meldna,axis=1)

pr('\n=== DeLong P: ML vs scores (reported: vs SOFA 0.660/0.088/0.011; vs MELD 0.132/0.017/0.002) ===')
for tp, px in [('7d',p7),('14d',p14),('28d',p28)]:
    y = xg_y[f'death_{tp}'].values
    for sc in ['SOFA_mod','MELD_Na']:
        aucs, p = delong_p(y, px, xgs[sc].values)
        pr(f'  {tp} vs {sc}: ML={aucs[0]:.3f} score={aucs[1]:.3f} P={p:.4f}')

pr('\n=== Fair-baseline LR (28d) ===')
SOFA_VARS = ['plt','tbil','map','vaso','crea']
medf=mi_tr_x[SOFA_VARS].median(); muf=mi_tr_x[SOFA_VARS].mean(); sdf=mi_tr_x[SOFA_VARS].std().replace(0,1)
fb = LogisticRegression(C=0.1, class_weight='balanced', max_iter=1000, random_state=42)
fb.fit((mi_tr_x[SOFA_VARS].fillna(medf)-muf)/sdf, mi_tr_y['death_28d'])
pfb = fb.predict_proba((xg_x[SOFA_VARS].fillna(medf)-muf)/sdf)[:,1]
aucs, p = delong_p(y28, p28, pfb)
rng=np.random.default_rng(42); boots=[]
for _ in range(1000):
    i=rng.integers(0,len(y28),len(y28))
    if len(np.unique(y28[i]))<2: continue
    boots.append(roc_auc_score(y28[i],pfb[i]))
ci=np.percentile(boots,[2.5,97.5])
pr(f'  full ML={aucs[0]:.3f} fair-baseline={aucs[1]:.3f} ({ci[0]:.3f}-{ci[1]:.3f})  DeLong P={p:.4f}')

pr('\n=== 14d recalibration exact values ===')
y14 = xg_y['death_14d'].values
eps=1e-15
logit = np.log(np.clip(p14,eps,1-eps)/(1-np.clip(p14,eps,1-eps)))
cal = LogisticRegression().fit(logit.reshape(-1,1), y14)
p_cal = cal.predict_proba(logit.reshape(-1,1))[:,1]
b0 = brier_score_loss(y14,p14); b1 = brier_score_loss(y14,p_cal)
pr(f'  before: Brier={b0:.4f}; after: Brier={b1:.4f}; recal params a={cal.intercept_[0]:.3f} b={cal.coef_[0][0]:.3f}')

pr('\n=== DCA fine grid 28d ext ===')
n=len(y28); prev=y28.mean()
region=[]
for t in np.arange(0.05,0.55,0.01):
    pred = p28>=t
    tp=((pred)&(y28==1)).sum(); fp=((pred)&(y28==0)).sum()
    nb = tp/n - fp/n*t/(1-t)
    nb_all = prev-(1-prev)*t/(1-t)
    if nb>nb_all and nb>0: region.append(t)
pr(f'  model beats treat-all & treat-none over thresholds {min(region):.2f}-{max(region):.2f}')
t=0.30; pred=p28>=t; tp=((pred)&(y28==1)).sum(); fp=((pred)&(y28==0)).sum()
nb=tp/n-fp/n*t/(1-t); nb_all=prev-(1-prev)*t/(1-t)
pr(f'  at 30%: model NB={nb:.4f}, treat-all NB={nb_all:.4f}, gain vs treat-all={nb-nb_all:.4f}')

pr('\n=== Bootstrap per-feature top-3 occurrence (200 resamples) ===')
import shap
rng = np.random.default_rng(42)
Xtr_raw = mi_tr_x[FEATS31].fillna(med)
for tp, ref in [('7d',['spo2','inr','lac']), ('14d',['hct','vaso','map']), ('28d',['hct','vaso','map'])]:
    ytr = mi_tr_y[f'death_{tp}'].values
    counts = {f:0 for f in ref}
    B=200
    for b in range(B):
        i = rng.integers(0,len(Xtr),len(Xtr))
        Xb = (Xtr_raw.iloc[i].reset_index(drop=True)-mu)/sd
        yb = ytr[i]
        if len(np.unique(yb))<2: continue
        m = mk(tp); m.fit(Xb, yb)
        if tp=='7d':
            sv = shap.TreeExplainer(m).shap_values(Xb)
            if isinstance(sv, list): sv = sv[1]
            imp = np.abs(sv).mean(0)
        else:
            imp = np.abs(m.coef_[0]) * np.abs(Xb.values).mean(0)
        top3 = set(pd.Series(imp, index=FEATS31).sort_values(ascending=False).index[:3])
        for f in ref:
            if f in top3: counts[f]+=1
    pr(f'  {tp}: ' + ', '.join(f'{f} in top-3 {counts[f]/B*100:.0f}%' for f in ref))
OUT.close()
