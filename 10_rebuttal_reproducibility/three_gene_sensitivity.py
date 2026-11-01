#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Three-gene (CD3G/CD6/THBS1) sensitivity analysis relative to the four-gene
(CD3G/PRTN3/CD6/THBS1) diagnostic signature for moderate-to-severe BPD.

The four-gene model is the pre-specified primary model with fixed coefficients
(CD3G -0.820, PRTN3 +0.610, CD6 -1.402, THBS1 +1.560; intercept -3.375).
The three-gene model is refitted by maximum likelihood in the training set and
then applied WITHOUT re-estimation to the Day-28 within-cohort subset and to
GSE220135, matching the fixed-coefficient strategy used for the primary model.

Reported per dataset: AUC (4-gene and 3-gene), delta-AUC, DeLong test,
continuous NRI and IDI (4-gene as "new" vs 3-gene as "old"; >0 means the
4-gene model is better, i.e. the loss incurred by removing PRTN3), with
95% bootstrap confidence intervals.

Usage
-----
    python three_gene_sensitivity.py --data-dir ./example_data --out ./output

Inputs (CSV, one row per sample; required columns: gsm, outcome, CD3G, PRTN3, CD6, THBS1)
---------------------------------------------------------------------------------------
    gse32472_day14_training_4genes.csv   training set, n = 63 (24 cases / 39 controls)
    gse32472_day28_4genes.csv            Day-28 within-cohort subset, n = 62 (24 / 38);
                                         optional column 'p4_stored' (four-gene probability
                                         as reported in the manuscript) is used if present
    gse220135_day14_4genes.csv           independent cohort, day 14, n = 68 (14 / 54)
    gse220135_day28_4genes.csv           independent cohort, day 28, n = 68 (13 / 55)

Data source: NCBI GEO GSE32472 and GSE220135. Expression values are log2
(microarray, GSE32472) or variance-stabilised counts (RNA-seq, GSE220135) for
the four genes only; see Supplementary Table S10 of the manuscript.

Requirements: python >= 3.9, numpy, pandas, scikit-learn, scipy
Outputs: three_gene_sensitivity.csv and .xlsx in the --out directory.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

# --------------------------------------------------------------------------
# model coefficients (pre-specified four-gene model, fitted in the training set)
# --------------------------------------------------------------------------
COEF4 = {'CD3G': -0.820, 'PRTN3': 0.610, 'CD6': -1.402, 'THBS1': 1.560}
INT4 = -3.375
GENES4 = ['CD3G', 'PRTN3', 'CD6', 'THBS1']
GENES3 = ['CD3G', 'CD6', 'THBS1']

B = 2000                       # bootstrap iterations
RNG = np.random.default_rng(20261101)


# --------------------------------------------------------------------------
# statistics: DeLong test, continuous NRI, IDI
# --------------------------------------------------------------------------
def _midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1)
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T + 1
    return T2


def _fast_delong(pred_sorted, m):
    n = pred_sorted.shape[1] - m
    k = pred_sorted.shape[0]
    pos, neg = pred_sorted[:, :m], pred_sorted[:, m:]
    tx = np.empty([k, m]); ty = np.empty([k, n]); tz = np.empty([k, m + n])
    for r in range(k):
        tx[r, :] = _midrank(pos[r, :])
        ty[r, :] = _midrank(neg[r, :])
        tz[r, :] = _midrank(pred_sorted[r, :])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    cov = np.cov(v01) / m + np.cov(v10) / n
    return aucs, cov


def delong(y, p1, p2):
    """DeLong test comparing AUC(p1) with AUC(p2). Returns (auc1, auc2, z, p)."""
    y = np.asarray(y); p1 = np.asarray(p1); p2 = np.asarray(p2)
    order = np.argsort(-y)
    m = int(y.sum())
    preds = np.vstack((p1, p2))[:, order]
    aucs, cov = _fast_delong(preds, m)
    L = np.array([[1.0, -1.0]])
    var = float((L @ cov @ L.T).item())
    if var <= 0:
        return aucs[0], aucs[1], np.nan, np.nan
    z = float((L @ aucs).item()) / np.sqrt(var)
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return aucs[0], aucs[1], z, p


def nri_idi(y, p_old, p_new):
    """Continuous NRI and IDI; >0 means the new model is better."""
    y = np.asarray(y); ev = y == 1; ne = ~ev
    up_ev = float((p_new[ev] > p_old[ev]).mean()); dn_ev = float((p_new[ev] < p_old[ev]).mean())
    up_ne = float((p_new[ne] > p_old[ne]).mean()); dn_ne = float((p_new[ne] < p_old[ne]).mean())
    nri = (up_ev - dn_ev) + (dn_ne - up_ne)
    idi = float((p_new[ev].mean() - p_old[ev].mean()) - (p_new[ne].mean() - p_old[ne].mean()))
    return nri, idi


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def linear_predictor(coefs, intercept, df):
    lp = np.full(len(df), float(intercept))
    for g, c in coefs.items():
        lp = lp + c * df[g].values.astype(float)
    return lp


def fit_three_gene(train):
    X = train[GENES3].values.astype(float)
    y = train['outcome'].values.astype(int)
    mdl = LogisticRegression(C=np.inf, solver='lbfgs', max_iter=20000)
    mdl.fit(X, y)
    return dict(zip(GENES3, mdl.coef_[0])), float(mdl.intercept_[0])


def load(path):
    df = pd.read_csv(path)
    need = ['outcome'] + GENES4
    missing = [c for c in need if c not in df.columns]
    if missing:
        sys.exit('ERROR: %s is missing columns: %s' % (path, missing))
    return df


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description='Three-gene sensitivity analysis')
    ap.add_argument('--data-dir', default='example_data', help='folder with the four input CSV files')
    ap.add_argument('--out', default='output', help='output folder')
    ap.add_argument('--boot', type=int, default=B, help='bootstrap iterations (default 2000)')
    args = ap.parse_args()

    nboot = args.boot
    ddir = args.data_dir
    os.makedirs(args.out, exist_ok=True)

    # ---- training set ----
    train = load(os.path.join(ddir, 'gse32472_day14_training_4genes.csv'))
    y_tr = train['outcome'].values.astype(int)

    c4_fit, i4_fit = None, None
    m4 = LogisticRegression(C=np.inf, solver='lbfgs', max_iter=20000)
    m4.fit(train[GENES4].values.astype(float), y_tr)
    print('Sanity check - four-gene coefficients refitted in the training set:')
    for g, v in zip(GENES4, m4.coef_[0]):
        print('   %-6s %.3f   (published %.3f)' % (g, v, COEF4[g]))
    print('   intercept %.3f   (published %.3f)' % (m4.intercept_[0], INT4))

    c3, i3 = fit_three_gene(train)
    print('Three-gene coefficients (training set): %s  intercept = %.3f\n'
          % ({g: round(v, 3) for g, v in c3.items()}, i3))

    p4_tr = sigmoid(linear_predictor(COEF4, INT4, train))
    p3_tr = sigmoid(i3 + sum(c * train[g].values.astype(float) for g, c in c3.items()))
    datasets = [('Training set (GSE32472, day 14)', y_tr, p4_tr, p3_tr)]

    # ---- Day-28 within-cohort ----
    d28 = load(os.path.join(ddir, 'gse32472_day28_4genes.csv'))
    y28 = d28['outcome'].values.astype(int)
    p3_28 = sigmoid(i3 + sum(c * d28[g].values.astype(float) for g, c in c3.items()))
    if 'p4_stored' in d28.columns and d28['p4_stored'].notna().all():
        p4_28 = d28['p4_stored'].values.astype(float)
        print('Day-28: using the four-gene probabilities stored in the input file')
    else:
        p4_28 = sigmoid(linear_predictor(COEF4, INT4, d28))
        print('Day-28: four-gene probabilities recomputed from the fixed coefficients')
    datasets.append(('Day-28 within-cohort (GSE32472, partially overlapping)', y28, p4_28, p3_28))

    # ---- GSE220135 ----
    for fn, label in (('gse220135_day14_4genes.csv', 'GSE220135, day 14'),
                      ('gse220135_day28_4genes.csv', 'GSE220135, day 28')):
        df = load(os.path.join(ddir, fn))
        y = df['outcome'].values.astype(int)
        p4 = sigmoid(linear_predictor(COEF4, INT4, df))
        p3 = sigmoid(i3 + sum(c * df[g].values.astype(float) for g, c in c3.items()))
        datasets.append((label, y, p4, p3))

    # ---- metrics ----
    rows = []
    print()
    for label, y, p4, p3 in datasets:
        y = np.asarray(y); p4 = np.asarray(p4, float); p3 = np.asarray(p3, float)
        ok = ~(np.isnan(p4) | np.isnan(p3) | np.isnan(y))
        y, p4, p3 = y[ok], p4[ok], p3[ok]
        n = len(y)
        a4 = roc_auc_score(y, p4); a3 = roc_auc_score(y, p3)
        _, _, z, pv = delong(y, p4, p3)
        nri, idi = nri_idi(y, p3, p4)

        nris, idis, daucs, a3s = [], [], [], []
        for _ in range(nboot):
            idx = RNG.integers(0, n, n)
            if len(np.unique(y[idx])) < 2:
                continue
            a_, b_ = nri_idi(y[idx], p3[idx], p4[idx])
            nris.append(a_); idis.append(b_)
            daucs.append(roc_auc_score(y[idx], p4[idx]) - roc_auc_score(y[idx], p3[idx]))
            a3s.append(roc_auc_score(y[idx], p3[idx]))
        q = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
        nlo, nhi = q(nris); ilo, ihi = q(idis); dlo, dhi = q(daucs); a3lo, a3hi = q(a3s)

        rows.append(dict(dataset=label, n=n, cases=int(y.sum()),
                         auc_4gene=round(a4, 4), auc_3gene=round(a3, 4),
                         auc3_lo=round(a3lo, 4), auc3_hi=round(a3hi, 4),
                         delta_auc=round(a4 - a3, 4), delta_lo=round(dlo, 4), delta_hi=round(dhi, 4),
                         delong_p=round(float(pv), 4) if not np.isnan(pv) else None,
                         nri=round(nri, 3), nri_lo=round(nlo, 3), nri_hi=round(nhi, 3),
                         idi=round(idi, 3), idi_lo=round(ilo, 3), idi_hi=round(ihi, 3),
                         n_boot=len(nris)))

        print('%-56s n=%3d' % (label, n))
        print('   AUC 4-gene %.4f | 3-gene %.4f (95%% CI %.4f-%.4f) | dAUC %+.4f (%.4f to %.4f)'
              % (a4, a3, a3lo, a3hi, a4 - a3, dlo, dhi))
        print('   DeLong P = %.4f | NRI %+.3f (%.3f, %.3f) | IDI %+.3f (%.3f, %.3f)'
              % (pv, nri, nlo, nhi, idi, ilo, ihi))

    df = pd.DataFrame(rows)
    out_csv = os.path.join(args.out, 'three_gene_sensitivity.csv')
    df.to_csv(out_csv, index=False, encoding='utf-8-sig')
    print('\n[saved]', out_csv)
    try:
        df.to_excel(os.path.join(args.out, 'three_gene_sensitivity.xlsx'), index=False)
        print('[saved]', os.path.join(args.out, 'three_gene_sensitivity.xlsx'))
    except Exception:
        pass
    print('\n' + df.to_string(index=False))


if __name__ == '__main__':
    main()
