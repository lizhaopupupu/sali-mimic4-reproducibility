# -*- coding: utf-8 -*-
"""Generate Supplementary Tables S5, S6, S7 (docx) with recomputed values."""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
from sklearn.metrics import roc_auc_score, brier_score_loss, average_precision_score, accuracy_score, recall_score
from sklearn.linear_model import LogisticRegression
from lightgbm import LGBMClassifier
import shap, os
from docx import Document
from docx.shared import Pt

FEATS31 = ['age','gender','max_temp','hr','rr','map','dbp','spo2','wbc','rbc','hb','hct','plt',
           'alt','ast','alp','alb','tbil','inr','bun','crea','na','k','ca','cl','hco3','glu',
           'lac','vent','vaso','crrt']
os.makedirs('submission/supplementary', exist_ok=True)

mi_tr_x = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_训练集_特征.xlsx')
mi_tr_y = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_训练集_预后.xlsx')
mi_va_x = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_验证集_特征.xlsx')
mi_va_y = pd.read_excel('../../01_原始数据/MIMIC队列/MIMIC_验证集_预后.xlsx')
xg_x = pd.concat([pd.read_excel('../../01_原始数据/单中心队列/单中心_临床训练集.xlsx'),
                  pd.read_excel('../../01_原始数据/单中心队列/单中心_临床验证集.xlsx')], ignore_index=True)
xg_y = pd.concat([pd.read_excel('../../01_原始数据/单中心队列/单中心_时间分类临床训练集预后.xlsx'),
                  pd.read_excel('../../01_原始数据/单中心队列/单中心_时间分类临床验证集预后.xlsx')], ignore_index=True)
med = mi_tr_x[FEATS31].median(); mu = mi_tr_x[FEATS31].mean(); sd = mi_tr_x[FEATS31].std().replace(0,1)
def prep(df): return (df[FEATS31].fillna(med)-mu)/sd
Xtr, Xva, Xxg = prep(mi_tr_x), prep(mi_va_x), prep(xg_x)

def mk(tp, seed=42):
    if tp=='7d':
        return LGBMClassifier(boosting_type='gbdt', class_weight='balanced', colsample_bytree=1.0,
                              learning_rate=0.05, max_depth=3, min_child_samples=20, n_estimators=200,
                              num_leaves=31, objective='binary', random_state=seed, verbose=-1)
    return LogisticRegression(C=0.1, class_weight='balanced', max_iter=1000, random_state=seed)

TOP5 = {'7d':['spo2','inr','lac','alt','crrt'], '14d':['hct','vaso','map','na','hr'], '28d':['hct','vaso','map','rbc','na']}

# ---------------- S5: bootstrap stability ----------------
rng = np.random.default_rng(42)
Xtr_raw = mi_tr_x[FEATS31].fillna(med)
B = 200
stab = {}
for tp in ['7d','14d','28d']:
    ytr = mi_tr_y[f'death_{tp}'].values
    c3 = {f:0 for f in TOP5[tp]}; c5 = {f:0 for f in TOP5[tp]}
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
        ranks = pd.Series(imp, index=FEATS31).sort_values(ascending=False)
        t3, t5 = set(ranks.index[:3]), set(ranks.index[:5])
        for f in TOP5[tp]:
            if f in t3: c3[f]+=1
            if f in t5: c5[f]+=1
    stab[tp] = (c3, c5)

doc = Document()
doc.add_paragraph('Table S5: Bootstrap-based ranking stability of the leading SHAP features (200 resamples of the MIMIC-IV training set)')
t = doc.add_table(rows=1, cols=4)
t.style = 'Light Grid Accent 1'
hdr = t.rows[0].cells
for j,h in enumerate(['Timepoint','Feature','Frequency in top-3 ranks (%)','Frequency in top-5 ranks (%)']):
    hdr[j].text = h
for tp in ['7d','14d','28d']:
    c3, c5 = stab[tp]
    for f in TOP5[tp]:
        row = t.add_row().cells
        row[0].text = tp; row[1].text = f
        row[2].text = f'{c3[f]/B*100:.0f}'; row[3].text = f'{c5[f]/B*100:.0f}'
p = doc.add_paragraph('Note: frequencies are the proportion of 200 bootstrap resamples in which the feature appeared among the top three (or top five) features ranked by mean |SHAP|. '
    'In the development-data export used for this analysis, vasopressor use, mechanical ventilation and CRRT were constant in the MIMIC-IV training file; '
    'frequencies for these features (vaso) are therefore not estimable (shown as 0) and will be updated on the complete export.')
doc.save('submission/supplementary/Table_S5.docx')
print('S5 saved')
for tp in ['7d','14d','28d']:
    print(tp, 'top3:', {k:f'{v/2:.0f}%' for k,v in stab[tp][0].items()}, 'top5:', {k:f'{v/2:.0f}%' for k,v in stab[tp][1].items()})

# ---------------- S6: DCA threshold-benefit ----------------
m28 = mk('28d'); m28.fit(Xtr, mi_tr_y['death_28d'])
p28 = m28.predict_proba(Xxg)[:,1]
y = xg_y['death_28d'].values; n=len(y); prev=y.mean()
doc = Document()
doc.add_paragraph('Table S6: Threshold–net-benefit pairs for decision-curve analysis of the 28-day model in external validation (n = 229)')
t = doc.add_table(rows=1, cols=4); t.style='Light Grid Accent 1'
for j,h in enumerate(['Threshold probability','Net benefit: model','Net benefit: treat all','Net benefit: treat none']):
    t.rows[0].cells[j].text = h
for th in np.arange(0.05,0.55,0.05):
    pred = p28>=th
    tp_=((pred)&(y==1)).sum(); fp=((pred)&(y==0)).sum()
    nb = tp_/n - fp/n*th/(1-th); nb_all = prev-(1-prev)*th/(1-th)
    row = t.add_row().cells
    row[0].text=f'{th:.2f}'; row[1].text=f'{nb:.4f}'; row[2].text=f'{nb_all:.4f}'; row[3].text='0.0000'
doc.add_paragraph('Net benefit = TP/n − FP/n × threshold/(1 − threshold). Treat-none net benefit is zero by definition.')
doc.save('submission/supplementary/Table_S6.docx')
print('S6 saved')

# ---------------- S7: classification metrics ----------------
doc = Document()
doc.add_paragraph('Table S7: AUPRC, Brier score, accuracy, sensitivity and specificity at threshold 0.5 in internal and external validation')
t = doc.add_table(rows=1, cols=7); t.style='Light Grid Accent 1'
for j,h in enumerate(['Timepoint','Dataset','AUPRC','Brier score','Accuracy','Sensitivity','Specificity']):
    t.rows[0].cells[j].text = h
for tp in ['7d','14d','28d']:
    m = mk(tp); m.fit(Xtr, mi_tr_y[f'death_{tp}'])
    for X, yy, tag in [(Xva, mi_va_y[f'death_{tp}'].values,'Internal validation (n = 431)'), (Xxg, xg_y[f'death_{tp}'].values,'External validation (n = 229)')]:
        p = m.predict_proba(X)[:,1]; yhat = p>=0.5
        row = t.add_row().cells
        row[0].text=tp; row[1].text=tag
        row[2].text=f'{average_precision_score(yy,p):.3f}'
        row[3].text=f'{brier_score_loss(yy,p):.3f}'
        row[4].text=f'{accuracy_score(yy,yhat):.3f}'
        row[5].text=f'{recall_score(yy,yhat):.3f}'
        row[6].text=f'{recall_score(yy,yhat,pos_label=0):.3f}'
doc.add_paragraph('AUPRC, area under the precision-recall curve.')
doc.save('submission/supplementary/Table_S7.docx')
print('S7 saved')
