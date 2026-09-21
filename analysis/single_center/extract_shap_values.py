# -*- coding: utf-8 -*-
"""
提取逐样本 SHAP 值（蜂群图数据）
基于 run_single_center_v3.py 的最优模型配置
"""
import pandas as pd, numpy as np, json, os, warnings
warnings.filterwarnings('ignore')
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.feature_selection import SelectKBest, chi2, f_classif, mutual_info_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.calibration import CalibratedClassifierCV
import shap

BASE = r'../../01_原始数据/单中心队列'
OUT  = r'../../03_运行结果'
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

def _mi_score(X, y):
    return mutual_info_classif(X, y, random_state=42)

def run():
    tr_x, tr_y, va_x, va_y = load()
    n_tr = len(tr_x)
    configs = {'death_7d': ('chi2','NB'), 'death_14d': ('mutual_info','NB'), 'death_28d': ('none','RF')}
    out = {}
    for h, (fs, clf) in configs.items():
        y_tr, y_va = tr_y[h].values, va_y[h].values
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
        # SHAP values
        if clf == 'RF':
            expl = shap.TreeExplainer(base)
            sv = expl.shap_values(Xv)
            if isinstance(sv, list):
                sv = np.array(sv[1])
            elif sv.ndim == 3:
                sv = sv[..., 1]
        else:
            pred_fn = lambda X: base.predict_proba(X)[:, 1]
            expl = shap.KernelExplainer(pred_fn, Xt[:50])
            sv = expl.shap_values(Xv)
            if isinstance(sv, list):
                sv = np.array(sv[1])
        sv = np.asarray(sv)
        # 保存：SHAP值矩阵、特征名、原始特征值（验证集标准化前，用于着色）
        # 对于着色，用验证集的原始特征值（但只保留选中的特征）
        va_x_sel = va_x[feat_names].values if sel is not None else va_x.values
        # 计算平均|SHAP|排序
        imp = np.abs(sv).mean(0)
        top_idx = np.argsort(-imp)
        top_names = [feat_names[i] for i in top_idx]
        top_imp = imp[top_idx].tolist()
        out[h] = {
            'feature_names': feat_names,
            'shap_values': sv.tolist(),           # [n_samples, n_features]
            'feature_values': va_x_sel.tolist(),  # [n_samples, n_features] 原始值
            'top_indices': top_idx.tolist(),
            'top_names': top_names,
            'top_mean_abs': [round(float(v),4) for v in top_imp],
        }
        print(f'[{h}] sv shape={sv.shape}, features={feat_names}, top5={top_names[:5]}')
    with open(os.path.join(OUT,'shap_values.json'),'w',encoding='utf-8') as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print('SHAP values saved to shap_values.json')

if __name__ == '__main__':
    run()
