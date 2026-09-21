# -*- coding: utf-8 -*-
"""
图 3/4 数据准备：
1) 三时点最终模型在验证集的预测概率、校准箱数据、ROC 坐标、DCA 坐标、SOFA/MELD 评分
2) 100 组合在验证集的 AUROC 矩阵（7/14/28 d 各一张热力图数据）
输出到 single_center_figs/fig_data.json + 100x100 AUC 矩阵 npy
"""
import pandas as pd, numpy as np, json, os, warnings
warnings.filterwarnings('ignore')
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.calibration import CalibratedClassifierCV
from sklearn.feature_selection import SelectKBest, chi2, mutual_info_classif
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.feature_selection import RFE, SelectFromModel, VarianceThreshold
from sklearn.linear_model import LassoCV
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss
import lightgbm as lgb
import xgboost as xgb

BASE = r'../../01_原始数据/单中心队列'
OUT = r'../../03_运行结果'
# NOTE (repository release): an absolute path to the original author's local working
# directory was hard-coded here and was used to create a figures directory.  It has been
# removed: this script writes all of its output to OUT (../../03_运行结果) and creates no
# directories outside the repository.  Nothing else in the script was changed.

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

def selectors():
    return {
        'all': lambda: None,
        'F-test': lambda: SelectKBest(f_classif_, k=6),
        'mutual_info': lambda: SelectKBest(_mi_score, k=6),
        'chi2': lambda: SelectKBest(chi2, k=6),
        'RFE-RF': lambda: RFE(RandomForestClassifier(n_estimators=50, random_state=42), n_features_to_select=6),
        'RFE-LR': lambda: RFE(LogisticRegression(max_iter=2000), n_features_to_select=6),
        'LASSO': lambda: SelectFromModel(LassoCV(cv=3, random_state=42, max_iter=5000)),
        'TreeImp': lambda: SelectFromModel(RandomForestClassifier(n_estimators=50, random_state=42)),
        'GBImp': lambda: SelectFromModel(GradientBoostingClassifier(n_estimators=50, random_state=42)),
        'Variance': lambda: VarianceThreshold(threshold=0.01),
    }
from sklearn.feature_selection import f_classif as f_classif_
from functools import partial

def _mi_score(X, y):
    from sklearn.feature_selection import mutual_info_classif
    return mutual_info_classif(X, y, random_state=42)

def classifiers():
    return {
        'LR': lambda: LogisticRegression(max_iter=3000, class_weight='balanced'),
        'RF': lambda: RandomForestClassifier(n_estimators=200, random_state=42, class_weight='balanced'),
        'XGB': lambda: xgb.XGBClassifier(n_estimators=200, eval_metric='logloss', random_state=42),
        'LGBM': lambda: lgb.LGBMClassifier(n_estimators=200, random_state=42, verbose=-1),
        'SVM': lambda: SVC(probability=True, class_weight='balanced', random_state=42),
        'KNN': lambda: KNeighborsClassifier(),
        'DT': lambda: DecisionTreeClassifier(random_state=42, class_weight='balanced'),
        'ET': lambda: ExtraTreesClassifier(n_estimators=200, random_state=42, class_weight='balanced'),
        'NB': lambda: GaussianNB(),
        'MLP': lambda: MLPClassifier(max_iter=2000, random_state=42),
    }

def transform(tr_x, y_tr, va_x, fs_name):
    scaler = StandardScaler().fit(tr_x)
    Xt, Xv = scaler.transform(tr_x), scaler.transform(va_x)
    sel = selectors()[fs_name]()
    if sel is not None:
        if fs_name == 'chi2':
            mm = MinMaxScaler().fit(Xt); Xt, Xv = mm.transform(Xt), mm.transform(Xv)
        sel.fit(Xt, y_tr)
        Xt, Xv = sel.transform(Xt), sel.transform(Xv)
    return Xt, Xv

def main():
    tr_x, tr_y, va_x, va_y = load()
    all_x = pd.concat([tr_x, va_x], ignore_index=True)
    sofa_all = modified_sofa(all_x); meld_all = meld_na(all_x)
    n_tr = len(tr_x)
    configs = {'death_7d': ('chi2','NB'), 'death_14d': ('mutual_info','NB'), 'death_28d': ('all','RF')}
    fig_data = {}
    auc_matrix = {}
    for h, (fs, clf) in configs.items():
        y_tr, y_va = tr_y[h].values, va_y[h].values
        sofa_va = np.asarray(sofa_all[n_tr:]); meld_va = np.asarray(meld_all[n_tr:])
        # 最终模型
        Xt, Xv = transform(tr_x, y_tr, va_x, fs)
        base = classifiers()[clf]()
        base.fit(Xt, y_tr)
        cal = CalibratedClassifierCV(base, method='sigmoid', cv=5)
        cal.fit(Xt, y_tr)
        p_cal = cal.predict_proba(Xv)[:,1]
        auc = roc_auc_score(y_va, p_cal)
        fpr, tpr, _ = roc_curve(y_va, p_cal)
        # 校准箱（10 等频箱）
        order = np.argsort(p_cal)
        bins = np.array_split(order, 10)
        cal_bins = []
        for b in bins:
            cal_bins.append({'pred': float(p_cal[b].mean()), 'obs': float(y_va[b].mean()), 'n': int(len(b))})
        # DCA
        dca = []
        for pt in np.arange(0.01, 0.61, 0.01):
            tp = ((p_cal>=pt)&(y_va==1)).sum(); fp = ((p_cal>=pt)&(y_va==0)).sum()
            n = len(y_va)
            nb_model = tp/n - (fp/n)*(pt/(1-pt))
            nb_all = y_va.mean() - (1-y_va.mean())*(pt/(1-pt))
            nb_none = 0.0
            dca.append({'thr': round(float(pt),2), 'model': round(float(nb_model),4),
                       'treat_all': round(float(nb_all),4), 'treat_none': nb_none})
        # SOFA/MELD ROC + DCA 参考
        sofa_auc = roc_auc_score(y_va, sofa_va); meld_auc = roc_auc_score(y_va, meld_va)
        fig_data[h] = {
            'fs': fs, 'clf': clf, 'val_auc': round(float(auc),3),
            'sofa_auc': round(float(sofa_auc),3), 'meld_auc': round(float(meld_auc),3),
            'fpr': [round(float(x),4) for x in fpr], 'tpr': [round(float(x),4) for x in tpr],
            'cal_bins': cal_bins, 'dca': dca,
            'sofa': [round(float(x),3) for x in sofa_va], 'meld': [round(float(x),3) for x in meld_va],
            'p': [round(float(x),4) for x in p_cal], 'y': [int(x) for x in y_va],
            'n_train': len(y_tr), 'n_val': len(y_va), 'mort_pct': round(float(y_va.mean()*100),1),
        }
        # 100 组合验证集 AUC 矩阵
        sel_names = list(selectors().keys()); clf_names = list(classifiers().keys())
        M = np.zeros((len(sel_names), len(clf_names)))
        for i, sn in enumerate(sel_names):
            for j, cn in enumerate(clf_names):
                try:
                    Xt2, Xv2 = transform(tr_x, y_tr, va_x, sn)
                    m2 = classifiers()[cn]()
                    m2.fit(Xt2, y_tr)
                    pp = m2.predict_proba(Xv2)[:,1]
                    M[i,j] = roc_auc_score(y_va, pp)
                except Exception:
                    M[i,j] = np.nan
        auc_matrix[h] = {'sel': sel_names, 'clf': clf_names, 'matrix': M.tolist()}
        print(f'[{h}] done, AUC={auc:.3f}, 100-combo matrix shape={M.shape}', flush=True)
    with open(os.path.join(OUT,'fig_data.json'),'w',encoding='utf-8') as f:
        json.dump({'fig3': fig_data, 'fig4': auc_matrix}, f, ensure_ascii=False, indent=1)
    print('ALL DONE -> fig_data.json', flush=True)

if __name__ == '__main__':
    main()
