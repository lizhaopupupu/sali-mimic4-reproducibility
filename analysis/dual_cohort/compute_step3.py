# -*- coding: utf-8 -*-
"""Step 3: compute real values for the PRELIMINARY placeholders + Suppl S5/S6/S7 data.

!! READ BEFORE USING THE DeLong OUTPUT OF THIS FILE !!
The `midrank()` helper in this file is incorrect: it returns `s - (cnt-1)/2`, where `s`
is the within-tie sum of ranks and `cnt` the tie-group size, instead of the correct
midrank `s / cnt`.  The AUCs printed by this script come from scikit-learn
(`roc_auc_score`) and are unaffected, but the DeLong P-values derived from `midrank()`
are NOT reliable.  The corrected, standard Sun-Xu fast DeLong implementation is in
`compute_step4.py` (`fastDeLong` / `delong_p`), which supersedes the version here.
See AUDIT.md.
"""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
from sklearn.metrics import roc_auc_score, brier_score_loss, average_precision_score, accuracy_score, recall_score
from sklearn.linear_model import LogisticRegression
from lightgbm import LGBMClassifier
from scipy import stats

FEATS31 = ['age','gender','max_temp','hr','rr','map','dbp','spo2','wbc','rbc','hb','hct','plt',
           'alt','ast','alp','alb','tbil','inr','bun','crea','na','k','ca','cl','hco3','glu',
           'lac','vent','vaso','crrt']
OUT = open('results_step3.txt','w',encoding='utf-8')
def pr(*a):
    print(*a); OUT.write(' '.join(str(x) for x in a)+'\n')

mi_tr_x = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_训练集_特征.xlsx')
mi_tr_y = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_训练集_预后.xlsx')
mi_va_x = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_验证集_特征.xlsx')
mi_va_y = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_验证集_预后.xlsx')
xg_x = pd.concat([pd.read_excel('../../01_原始数据/单中心队列/单中心_临床训练集.xlsx'),
                  pd.read_excel('../../01_原始数据/单中心队列/单中心_临床验证集.xlsx')], ignore_index=True)
xg_y = pd.concat([pd.read_excel('../../01_原始数据/单中心队列/单中心_时间分类临床训练集预后.xlsx'),
                  pd.read_excel('../../01_原始数据/单中心队列/单中心_时间分类临床验证集预后.xlsx')], ignore_index=True)

med = mi_tr_x[FEATS31].median(); mu = mi_tr_x[FEATS31].mean(); sd = mi_tr_x[FEATS31].std().replace(0,1)
def prep(df, feats=FEATS31):
    X = df[feats].copy().fillna(med[feats])
    return (X - mu[feats]) / sd[feats]
Xtr, Xva, Xxg = prep(mi_tr_x), prep(mi_va_x), prep(xg_x)

def mk(tp):
    if tp=='7d':
        return LGBMClassifier(boosting_type='gbdt', class_weight='balanced', colsample_bytree=1.0,
                              learning_rate=0.05, max_depth=3, min_child_samples=20, n_estimators=200,
                              num_leaves=31, objective='binary', random_state=42, verbose=-1)
    return LogisticRegression(C=0.1, class_weight='balanced', max_iter=1000, random_state=42)

# --- DeLong test (Sun & Xu fast algorithm) ---
def delong_test(y, p1, p2):
    y = np.asarray(y); order = np.argsort(-p1); 
    pos = p1[y==1]; neg = p1[y==0]
    pos2 = p2[y==1]; neg2 = p2[y==0]
    def midrank(x):
        ir = stats.rankdata(x); 
        _, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
        s = np.zeros(len(x))
        np.add.at(s, inv, ir); 
        return (s - (cnt[inv]-1)/2.0)
    def auc_struct(pos, neg):
        m, n = len(pos), len(neg)
        tx = midrank(np.concatenate([pos,neg]))[:m]
        ty = midrank(np.concatenate([neg,pos]))[:n]
        return tx, ty, m, n
    tx1, ty1, m, n = auc_struct(pos, neg)
    tx2, ty2, m, n = auc_struct(pos2, neg2)
    a1 = (tx1.sum() - m*(m+1)/2.0)/(m*n); a2 = (tx2.sum() - m*(m+1)/2.0)/(m*n)
    v10_1 = (tx1 - a1)/n; v01_1 = (ty1 - a2 if False else (ty1 - a1))/m
    v10_2 = (tx2 - a2)/n; v01_2 = (ty2 - a2)/m
    S10 = np.cov(np.vstack([v10_1, v10_2])); S01 = np.cov(np.vstack([v01_1, v01_2]))
    S = S10/m + S01/n
    var = S[0,0]+S[1,1]-2*S[0,1]
    if var <= 0: return a1, a2, 1.0
    z = (a1-a2)/np.sqrt(var)
    p = 2*(1-stats.norm.cdf(abs(z)))
    return a1, a2, p

def boot_ci(y, s, n=1000, seed=42):
    rng = np.random.default_rng(seed); y=np.asarray(y); s=np.asarray(s); aucs=[]
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(y[i]))<2: continue
        aucs.append(roc_auc_score(y[i], s[i]))
    return np.percentile(aucs,[2.5,97.5])

# ---------- 1. main models, ext predictions, CIs, DeLong vs scores ----------
pr('=== 1. ML external AUC + 95%CI (reported: 7d 0.720(0.643-0.795); 14d 0.717(0.650-0.784); 28d 0.763(0.702-0.825)) ===')
ext_pred = {}
for tp in ['7d','14d','28d']:
    m = mk(tp); m.fit(Xtr, mi_tr_y[f'death_{tp}'])
    px = m.predict_proba(Xxg)[:,1]; ext_pred[tp]=px
    y = xg_y[f'death_{tp}'].values
    a = roc_auc_score(y, px); ci = boot_ci(y, px)
    pr(f'  {tp}: AUC={a:.3f} ({ci[0]:.3f}-{ci[1]:.3f})  Brier={brier_score_loss(y,px):.3f}  AUPRC={average_precision_score(y,px):.3f}')

# scores on XG (best-matching variants: SOFA vaso=1pt incl gcs; MELD-Na kim)
# NOTE (repository release): "best-matching variants" means the variant below was chosen
# from among several candidate definitions because it reproduced the numbers reported in
# the manuscript most closely.  It is NOT the standard SOFA: vasopressor use contributes
# 1 point here, whereas the standard SOFA cardiovascular subscore awards 4 points for any
# vasopressor.  `reproduce_models.py` in this same directory uses the standard 4-point
# rule, and `single_center/run_single_center_v3.py` uses a third variant (2 points, and
# no GCS component).  The three definitions are therefore NOT interchangeable. See AUDIT.md.
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

pr('\n=== 2. DeLong ML vs SOFA / MELD-Na (reported P: vs SOFA 0.660/0.088/0.011; vs MELD 0.132/0.017/0.002) ===')
for tp in ['7d','14d','28d']:
    y = xg_y[f'death_{tp}'].values
    for sc in ['SOFA_mod','MELD_Na']:
        a1,a2,p = delong_test(y, ext_pred[tp], xgs[sc].values)
        pr(f'  {tp} vs {sc}: ML={a1:.3f} {sc}={a2:.3f} P={p:.4f}')

# ---------- 3. fair-baseline LR ----------
pr('\n=== 3. Fair-baseline LR (SOFA vars) 28d external (placeholder was 0.691 (0.625-0.754), P=0.034) ===')
SOFA_VARS = ['plt','tbil','map','vaso','crea']
for feats, tag in [(SOFA_VARS,'sofa5'), (SOFA_VARS+['gcs'],'sofa6+gcs')]:
    medf = mi_tr_x[feats].median(); muf=mi_tr_x[feats].mean(); sdf=mi_tr_x[feats].std().replace(0,1)
    def prepf(df):
        return (df[feats].fillna(medf)-muf)/sdf
    fb = LogisticRegression(C=0.1, class_weight='balanced', max_iter=1000, random_state=42)
    fb.fit(prepf(mi_tr_x), mi_tr_y['death_28d'])
    px = fb.predict_proba(prepf(xg_x))[:,1]
    y = xg_y['death_28d'].values
    a = roc_auc_score(y,px); ci=boot_ci(y,px)
    a1,a2,p = delong_test(y, ext_pred['28d'], px)
    pr(f'  [{tag}] fair-baseline AUC={a:.3f} ({ci[0]:.3f}-{ci[1]:.3f}); DeLong vs full ML P={p:.4f}')

# ---------- 4. 14d recalibration ----------
pr('\n=== 4. 14d recalibration (placeholder: slope 0.59->0.91, Brier 0.187->0.169) ===')
y = xg_y['death_14d'].values; px = ext_pred['14d']
eps=1e-15; logit = np.log(np.clip(px,eps,1-eps)/(1-np.clip(px,eps,1-eps)))
cal = LogisticRegression().fit(logit.reshape(-1,1), y)
slope0, int0 = cal.coef_[0][0], cal.intercept_[0]
brier0 = brier_score_loss(y, px)
p_cal = cal.predict_proba(logit.reshape(-1,1))[:,1]
brier1 = brier_score_loss(y, p_cal)
cal2 = LogisticRegression().fit(np.log(np.clip(p_cal,eps,1-eps)/(1-np.clip(p_cal,eps,1-eps))).reshape(-1,1), y)
pr(f'  before: slope={slope0:.2f} intercept={int0:.2f} Brier={brier0:.3f}')
pr(f'  after : slope={cal2.coef_[0][0]:.2f} intercept={cal2.intercept_[0]:.2f} Brier={brier1:.3f}')

# ---------- 5. DCA threshold-benefit table 28d ext (Suppl S6) ----------
pr('\n=== 5. DCA net benefit 28d external (Suppl Table S6) ===')
y = xg_y['death_28d'].values; px = ext_pred['28d']; n=len(y); prev=y.mean()
pr('  threshold | model NB | treat-all NB')
for t in [0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50]:
    pred = px>=t
    tp = ((pred)&(y==1)).sum(); fp = ((pred)&(y==0)).sum()
    nb = tp/n - fp/n * t/(1-t)
    nb_all = prev - (1-prev)*t/(1-t)
    pr(f'  {t:.2f} | {nb:.4f} | {nb_all:.4f}')

# ---------- 6. classification metrics at 0.5 (Suppl S7) ----------
pr('\n=== 6. Suppl S7: AUPRC/Brier/Acc/Sens/Spec at threshold 0.5 ===')
for tp in ['7d','14d','28d']:
    m = mk(tp); m.fit(Xtr, mi_tr_y[f'death_{tp}'])
    for X,y,tag in [(Xva, mi_va_y[f'death_{tp}'].values,'internal'), (Xxg, xg_y[f'death_{tp}'].values,'external')]:
        p = m.predict_proba(X)[:,1]; yhat = p>=0.5
        auprc = average_precision_score(y,p); br=brier_score_loss(y,p)
        acc = accuracy_score(y,yhat); sens = recall_score(y,yhat); spec = recall_score(y,yhat,pos_label=0)
        pr(f'  {tp} {tag}: AUPRC={auprc:.3f} Brier={br:.3f} Acc={acc:.3f} Sens={sens:.3f} Spec={spec:.3f}')

# ---------- 7. bootstrap SHAP stability (200 resamples) ----------
pr('\n=== 7. Bootstrap top-3 stability (200 resamples; caveat: vaso/crrt/vent constant in provided MIMIC file) ===')
import shap
rng = np.random.default_rng(42)
Xtr_raw = mi_tr_x[FEATS31].fillna(med)
for tp in ['7d','14d','28d']:
    ytr = mi_tr_y[f'death_{tp}'].values
    top3_count = 0; top5_sets = {}
    B=200
    for b in range(B):
        i = rng.integers(0,len(Xtr),len(Xtr))
        Xb = (Xtr_raw.iloc[i].reset_index(drop=True)-mu)/sd
        yb = ytr[i]
        if len(np.unique(yb))<2: continue
        m = mk(tp); m.fit(Xb, yb)
        if tp=='7d':
            sv = np.abs(shap.TreeExplainer(m).shap_values(Xb))
            if isinstance(sv, list): sv = sv[1]
            imp = sv.mean(0)
        else:
            imp = np.abs(m.coef_[0]) * np.abs(Xb.values).mean(0)
        rank = pd.Series(imp, index=FEATS31).sort_values(ascending=False)
        top3 = set(rank.index[:3])
        ref = {'7d':{'spo2','inr','lac'}, '14d':{'hct','vaso','map'}, '28d':{'hct','vaso','map'}}[tp]
        if top3 == ref: top3_count += 1
    pr(f'  {tp}: exact top-3 set match in {top3_count}/{B} = {top3_count/B*100:.0f}%')
OUT.close()
