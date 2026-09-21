# -*- coding: utf-8 -*-
"""
v3: 仅评估已选定的最优模型（不重跑 100 组合）
7d=chi2+NB, 14d=mutual_info+NB, 28d=none+RF
修复: DeLong 用 numpy 数组; SHAP 用 KernelExplainer 兼容 NB
"""
import pandas as pd, numpy as np, json, warnings, os
warnings.filterwarnings('ignore')
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.calibration import CalibratedClassifierCV
from sklearn.feature_selection import SelectKBest, chi2, f_classif, mutual_info_classif

def _mi_score(X, y):
    return mutual_info_classif(X, y, random_state=42)
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import roc_auc_score, brier_score_loss, average_precision_score, accuracy_score, recall_score, f1_score
from scipy.stats import norm, linregress
from scipy.special import logit

BASE = r'../../01_原始数据/单中心队列'
OUT = r'../../03_运行结果'
FEATS31 = ['age','gender','max_temp','hr','rr','map','dbp','spo2','wbc','rbc','hb','hct','plt',
           'alt','ast','alp','alb','tbil','inr','bun','crea','na','k','ca','cl','hco3','glu',
           'lac','vent','vaso','crrt']

def load():
    tr_x = pd.read_excel(os.path.join(BASE,'单中心_临床训练集.xlsx'))
    va_x = pd.read_excel(os.path.join(BASE,'单中心_临床验证集.xlsx'))
    tr_y = pd.read_excel(os.path.join(BASE,'单中心_时间分类临床训练集预后.xlsx'))
    va_y = pd.read_excel(os.path.join(BASE,'单中心_时间分类临床验证集预后.xlsx'))
    tr_x, va_x = tr_x[FEATS31].copy(), va_x[FEATS31].copy()
    med = tr_x.median()
    return tr_x.fillna(med), tr_y, va_x.fillna(med), va_y

def modified_sofa(df):
    circ = np.where((df['map'].fillna(70)<70)|(df['vaso'].fillna(0)==1), 2, 0).astype(float)
    plt_ = df['plt'].fillna(150)
    coag = np.select([plt_<20, plt_<50, plt_<100, plt_<150], [4,3,2,1], default=0).astype(float)
    tb = df['tbil'].fillna(10)
    liver = np.select([tb>205, tb>102, tb>34.2, tb>20.5], [4,3,2,1], default=0).astype(float)
    cr = df['crea'].fillna(80)/88.4
    kidney = np.select([df['crrt'].fillna(0)==1, cr>5, cr>3.5, cr>2, cr>1.2], [4,4,3,2,1], default=0).astype(float)
    return circ + coag + liver + kidney

def meld_na(df):
    tb_mg = np.clip(df['tbil'].fillna(20),1,400)/17.1
    inr = np.clip(df['inr'].fillna(1.5),1,5)
    cr = np.clip(df['crea'].fillna(100)/88.4,1,4)
    meld = 3.78*np.log(tb_mg)+11.2*np.log(inr)+9.57*np.log(cr)+6.43
    na = np.clip(df['na'].fillna(138),125,140)
    return np.clip(meld+1.32*(135-na)-0.033*meld*(135-na),6,40)

def delong_p(a1, a2, y, p1, p2, n_boot=5000, seed=42):
    # NOTE (repository release): despite its name this is NOT the DeLong test.  It is a
    # paired bootstrap test of the AUC difference (5000 resamples, seed 42).  The DeLong
    # test proper is implemented in ../dual_cohort/compute_step4.py.  The JSON keys this
    # function feeds (`p_vs_sofa`, `p_vs_meld`) are therefore bootstrap P-values.
    # See AUDIT.md.
    y, p1, p2 = np.asarray(y), np.asarray(p1), np.asarray(p2)
    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2: continue
        try:
            diffs.append(roc_auc_score(y[idx], p1[idx]) - roc_auc_score(y[idx], p2[idx]))
        except Exception:
            continue
    diffs = np.array(diffs)
    if len(diffs) < 200: return float('nan')
    sd = diffs.std(); z = (a1-a2)/sd if sd > 0 else 0
    return float(2*(1-norm.cdf(abs(z))))

def boot_ci(y, p, n_boot=5000, seed=42):
    y, p = np.asarray(y), np.asarray(p)
    rng = np.random.default_rng(seed)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2: continue
        try: aucs.append(roc_auc_score(y[idx], p[idx]))
        except Exception: continue
    q = np.percentile(aucs, [2.5, 97.5])
    return round(float(q[0]),3), round(float(q[1]),3)

def run():
    tr_x, tr_y, va_x, va_y = load()
    all_x = pd.concat([tr_x, va_x], ignore_index=True)
    sofa_all, meld_all = modified_sofa(all_x), meld_na(all_x)
    n_tr = len(tr_x)
    configs = {'death_7d': ('chi2','NB'), 'death_14d': ('mutual_info','NB'), 'death_28d': ('none','RF')}
    out = {}
    for h, (fs, clf) in configs.items():
        y_tr, y_va = tr_y[h].values, va_y[h].values
        sofa_va = np.asarray(sofa_all[n_tr:]); meld_va = np.asarray(meld_all[n_tr:])
        # 特征选择 + 标准化
        scaler = StandardScaler().fit(tr_x)
        Xt, Xv = scaler.transform(tr_x), scaler.transform(va_x)
        sel = None
        if fs == 'chi2':
            mm = MinMaxScaler().fit(Xt); Xt, Xv = mm.transform(Xt), mm.transform(Xv)
            sel = SelectKBest(chi2, k=6)
        elif fs == 'mutual_info':
            sel = SelectKBest(_mi_score, k=6)
        feat_names = list(tr_x.columns)
        if sel is not None:
            sel.fit(Xt, y_tr)
            mask = sel.get_support()
            feat_names = [n for n,m in zip(feat_names, mask) if m]
            Xt, Xv = sel.transform(Xt), sel.transform(Xv)
        base = RandomForestClassifier(n_estimators=500, random_state=42, class_weight='balanced') if clf=='RF' else GaussianNB()
        base.fit(Xt, y_tr)
        # Platt 校准
        cal = CalibratedClassifierCV(base, method='sigmoid', cv=5)
        cal.fit(Xt, y_tr)
        p_cal = cal.predict_proba(Xv)[:,1]
        auc_cal = roc_auc_score(y_va, p_cal)
        lo, hi = boot_ci(y_va, p_cal)
        brier = brier_score_loss(y_va, p_cal)
        auprc = average_precision_score(y_va, p_cal)
        pv = np.clip(p_cal, 1e-4, 1-1e-4)
        lr_res = linregress(logit(pv), y_va)
        pred = (p_cal >= 0.5).astype(int)
        acc = accuracy_score(y_va, pred); sens = recall_score(y_va, pred)
        spec = recall_score(y_va, pred, pos_label=0); f1 = f1_score(y_va, pred)
        auc_sofa = roc_auc_score(y_va, sofa_va); auc_meld = roc_auc_score(y_va, meld_va)
        p_sofa = delong_p(auc_cal, auc_sofa, y_va, p_cal, sofa_va)
        p_meld = delong_p(auc_cal, auc_meld, y_va, p_cal, meld_va)
        nb = {}
        for thr in [0.1,0.2,0.3,0.4,0.5]:
            tp = ((p_cal>=thr)&(y_va==1)).sum(); fp = ((p_cal>=thr)&(y_va==0)).sum()
            n = len(y_va)
            nb[str(thr)] = round(float(tp/n-(fp/n)*(thr/(1-thr))),4)
        # SHAP: RF->TreeExplainer, NB->KernelExplainer
        shap_top = {}
        try:
            import shap
            if clf == 'RF':
                expl = shap.TreeExplainer(base)
                sv = expl.shap_values(Xv)
                if isinstance(sv, list):
                    sv = sv[1]
                elif sv.ndim == 3:   # shap>=0.45 RF 二分类返回 (n, p, 2)
                    sv = sv[..., 1]
            else:
                # NB 等：包装为正类概率一维输出
                pred_fn = lambda X: base.predict_proba(X)[:, 1]
                expl = shap.KernelExplainer(pred_fn, Xt[:50])
                sv = expl.shap_values(Xv)
                if isinstance(sv, list):
                    sv = sv[1]
            imp = np.abs(sv).mean(0)
            shap_top = {str(n): round(float(v),4) for n,v in sorted(zip(feat_names, imp), key=lambda x:-x[1])[:10]}
        except Exception as e:
            shap_top = {'error': str(e)[:150]}
        out[h] = {
            'best_fs': fs, 'best_clf': clf, 'selected_features': feat_names,
            'val_auc_cal': round(float(auc_cal),3), 'auc_95ci': [lo,hi],
            'sofa_auc': round(float(auc_sofa),3), 'meld_auc': round(float(auc_meld),3),
            'p_vs_sofa': round(float(p_sofa),4), 'p_vs_meld': round(float(p_meld),4),
            'brier': round(float(brier),3), 'auprc': round(float(auprc),3),
            'calib_slope': round(float(lr_res.slope),2), 'calib_intercept': round(float(lr_res.intercept),2),
            'acc': round(float(acc),3), 'sens': round(float(sens),3), 'spec': round(float(spec),3), 'f1': round(float(f1),3),
            'dca_nb': nb,
            'deaths_train': int(y_tr.sum()), 'deaths_val': int(y_va.sum()),
            'n_train': len(y_tr), 'n_val': len(y_va),
            'mortality_val_pct': round(float(y_va.mean()*100),1),
            'shap_top10': shap_top,
        }
        print(f'[{h}] {fs}+{clf} ValAUC={auc_cal:.3f}[{lo}-{hi}] SOFA={auc_sofa:.3f}(P={p_sofa:.4f}) MELD={auc_meld:.3f}(P={p_meld:.4f}) Brier={brier:.3f} slope={lr_res.slope:.2f}')
        print(f'   SHAP top5: {list(shap_top.items())[:5]}')
    with open(os.path.join(OUT,'results_single_center_v3.json'),'w',encoding='utf-8') as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print('DONE')

if __name__ == '__main__':
    run()
