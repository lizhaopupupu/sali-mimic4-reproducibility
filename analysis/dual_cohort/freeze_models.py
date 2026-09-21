# -*- coding: utf-8 -*-
"""
Fit the three final models and write frozen objects to ../../models/.

The estimator and its hyperparameters are those of Supplementary Table S4; the
preprocessing (training-set median imputation, training-set mean/SD standardisation) is
the one described in Methods Section 2.9. Everything is fixed by random_state=42.

    python freeze_models.py

writes
    ../../models/sali_7d_lightgbm.pkl       frozen 7-day LightGBM model
    ../../models/sali_14d_logreg.pkl        frozen 14-day logistic-regression model
    ../../models/sali_28d_logreg.pkl        frozen 28-day logistic-regression model
    ../../models/preprocessing_parameters.json   imputation medians, scaling means and SDs
    ../../models/MODEL_CARD.json            hyperparameters, feature order, provenance,
                                            and the discrimination obtained here

READ ../../models/MODEL_CARD.json BEFORE USING THESE OBJECTS. They were refitted from the
archived development export, in which vent / vaso / crrt are constant at zero. They are
reproducible artefacts of this repository, not the objects that produced the published
performance estimates, and their AUCs differ from the published ones.
"""
import json
import os
import pickle

import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

M = r'../../01_原始数据/MIMIC队列'
X = r'../../01_原始数据/单中心队列'
OUT = r'../../models'

FEATS31 = ['age', 'gender', 'max_temp', 'hr', 'rr', 'map', 'dbp', 'spo2', 'wbc', 'rbc', 'hb', 'hct', 'plt',
           'alt', 'ast', 'alp', 'alb', 'tbil', 'inr', 'bun', 'crea', 'na', 'k', 'ca', 'cl', 'hco3', 'glu',
           'lac', 'vent', 'vaso', 'crrt']

HYPER = {
    '7d': {'estimator': 'lightgbm.LGBMClassifier',
           'params': dict(boosting_type='gbdt', class_weight='balanced', colsample_bytree=1.0,
                          learning_rate=0.05, max_depth=3, min_child_samples=20, n_estimators=200,
                          num_leaves=31, objective='binary', random_state=42, verbose=-1)},
    '14d': {'estimator': 'sklearn.linear_model.LogisticRegression',
            'params': dict(C=0.1, class_weight='balanced', max_iter=1000, random_state=42)},
    '28d': {'estimator': 'sklearn.linear_model.LogisticRegression',
            'params': dict(C=0.1, class_weight='balanced', max_iter=1000, random_state=42)},
}


def harmonize(df):
    df = df.copy()
    if df['gender'].dtype == object:
        df['gender'] = df['gender'].astype(str).str.upper().str[0].map({'M': 1, 'F': 0})
    return df


def build(tp):
    spec = HYPER[tp]
    cls = LGBMClassifier if spec['estimator'].startswith('lightgbm') else LogisticRegression
    return cls(**spec['params'])


def main():
    mi_tr_x = harmonize(pd.read_excel(M + r'/MIMIC_训练集_特征.xlsx'))
    mi_tr_y = pd.read_excel(M + r'/MIMIC_训练集_预后.xlsx')
    mi_va_x = harmonize(pd.read_excel(M + r'/MIMIC_验证集_特征.xlsx'))
    mi_va_y = pd.read_excel(M + r'/MIMIC_验证集_预后.xlsx')
    xg_x = harmonize(pd.concat([pd.read_excel(X + r'/单中心_临床训练集.xlsx'),
                                pd.read_excel(X + r'/单中心_临床验证集.xlsx')], ignore_index=True))
    xg_y = pd.concat([pd.read_excel(X + r'/单中心_时间分类临床训练集预后.xlsx'),
                      pd.read_excel(X + r'/单中心_时间分类临床验证集预后.xlsx')], ignore_index=True)

    med = mi_tr_x[FEATS31].median()
    mu = mi_tr_x[FEATS31].mean()
    sd = mi_tr_x[FEATS31].std().replace(0, 1)
    prep = lambda d: (d[FEATS31].fillna(med) - mu) / sd
    Xtr, Xva, Xxg = prep(mi_tr_x), prep(mi_va_x), prep(xg_x)

    os.makedirs(OUT, exist_ok=True)
    results = {}
    for tp in ('7d', '14d', '28d'):
        mdl = build(tp)
        mdl.fit(Xtr, mi_tr_y[f'death_{tp}'])
        path = os.path.join(OUT, f'sali_{tp}_' + ('lightgbm' if tp == '7d' else 'logreg') + '.pkl')
        with open(path, 'wb') as f:
            pickle.dump(mdl, f, protocol=4)
        pv = mdl.predict_proba(Xva)[:, 1]
        px = mdl.predict_proba(Xxg)[:, 1]
        av = roc_auc_score(mi_va_y[f'death_{tp}'], pv)
        ax = roc_auc_score(xg_y[f'death_{tp}'], px)
        results[tp] = {'object': os.path.basename(path),
                       'internal_validation_auc': round(float(av), 3),
                       'external_auc': round(float(ax), 3)}
        print(f'  {tp}: wrote {os.path.basename(path)}  internal-val AUC={av:.3f}  external AUC={ax:.3f}')

    with open(os.path.join(OUT, 'preprocessing_parameters.json'), 'w', encoding='utf-8') as f:
        json.dump({'note': 'Computed on the MIMIC-IV training set (n=1119); applied unchanged to the '
                           'internal-validation and external cohorts.',
                   'order': FEATS31,
                   'imputation_median': {k: float(med[k]) for k in FEATS31},
                   'scaling_mean': {k: float(mu[k]) for k in FEATS31},
                   'scaling_sd': {k: float(sd[k]) for k in FEATS31}}, f, indent=2, ensure_ascii=False)
    print('  wrote preprocessing_parameters.json')

    card = {
        'title': 'Frozen SALI mortality models',
        'repository': 'https://github.com/lizhaopupupu/sali-mimic4-reproducibility',
        'generated_by': 'analysis/dual_cohort/freeze_models.py',
        'random_state': 42,
        'features': FEATS31,
        'preprocessing': 'median imputation with the training-set median, then z-score standardisation '
                         'with the training-set mean and SD; both fitted on the MIMIC-IV training set only '
                         'and applied unchanged to the validation and external cohorts.',
        'models': {tp: dict(HYPER[tp], **results[tp]) for tp in results},
        'PROVENANCE': [
            'These objects were refitted from the archived development export held for this project.',
            'In that export the columns vent, vaso and crrt are constant at zero, which contradicts the '
            '71.2% / 60.1% / 14.5% reported in Table 1.',
            'These objects are therefore NOT the models that produced the performance estimates published '
            'in the manuscript, and the AUCs obtained from them differ from the published ones '
            '(7-day external AUC here 0.760 against 0.720 published).',
            'To regenerate models from a correct export, run '
            'code/MIMIC4_SALI_cohort_extraction.sql against MIMIC-IV v3.1 and then re-run this script.',
        ],
    }
    with open(os.path.join(OUT, 'MODEL_CARD.json'), 'w', encoding='utf-8') as f:
        json.dump(card, f, indent=2, ensure_ascii=False)
    print('  wrote MODEL_CARD.json')
    print()
    print('Published for comparison: internal-val 0.840 / 0.841 / 0.850; external 0.720 / 0.717 / 0.763')


if __name__ == '__main__':
    main()
