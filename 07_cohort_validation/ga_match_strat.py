# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/ga_match_strat.py
#  Status : ORIGINAL (as used in the analysis)
#  Paths  : the authors' working copy used absolute paths on a local Windows
#           machine. They have been replaced here by P("<relative path>")
#           calls. Set BPD_REPO to the repository root (or run the script from
#           the repository root) and place inputs as described in docs/DATA.md.
#  Note   : the statistical code itself is unchanged.
# ==========================================================================
import os as _os
def P(rel):
    """Resolve a repository-relative path (placeholder for local absolute paths)."""
    root = _os.environ.get("BPD_REPO") or _os.path.abspath(
        _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))
    return _os.path.join(root, *rel.replace("\\", "/").split("/"))

"""GA-matched / GA-stratified / restricted / IPW sensitivity analyses"""
import pandas as pd, numpy as np, json
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from scipy import stats

m = pd.read_csv('training_set_clinical_genotype.csv')
G = ['CD3G', 'PRTN3', 'CD6', 'THBS1']
y = m['outcome'].values.astype(int)


def auc_boot(sub, n_boot=1000, seed=1):
    """AUC of the four-gene model fitted and evaluated within a subset (apparent), with bootstrap CI."""
    yy = sub['outcome'].values.astype(int)
    X = sub[G].values
    if len(np.unique(yy)) < 2:
        return np.nan, (float('nan'), float('nan')), len(sub), int((yy == 1).sum()), int((yy == 0).sum())
    mod = LogisticRegression(C=np.inf, max_iter=8000).fit(X, yy)
    a = roc_auc_score(yy, mod.predict_proba(X)[:, 1])
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(yy), len(yy))
        if len(np.unique(yy[idx])) < 2:
            continue
        try:
            mm = LogisticRegression(C=np.inf, max_iter=8000).fit(X[idx], yy[idx])
            bs.append(roc_auc_score(yy[idx], mm.predict_proba(X[idx])[:, 1]))
        except Exception:
            continue
    ci = np.percentile(bs, [2.5, 97.5]) if len(bs) > 50 else (np.nan, np.nan)
    return a, (float(ci[0]), float(ci[1])), len(sub), int((yy == 1).sum()), int((yy == 0).sum())


def weighted_auc(y_, scores, w):
    """Weighted c-statistic (IPW-adjusted AUC)."""
    y_ = np.asarray(y_); scores = np.asarray(scores); w = np.asarray(w)
    pos, neg = scores[y_ == 1], scores[y_ == 0]
    wp, wn = w[y_ == 1], w[y_ == 0]
    num = 0.0; den = 0.0
    for s, ws in zip(pos, wp):
        cmp_ = (s > neg).astype(float) + 0.5 * (s == neg).astype(float)
        num += ws * np.sum(wn * cmp_)
        den += ws * np.sum(wn)
    return num / den


res = []

# ---------- 0) overlap diagnostics ----------
c, n = m[m['outcome'] == 1], m[m['outcome'] == 0]
res.append({'analysis': 'GA distribution (cases)', 'detail': 'weeks: ' + ', '.join(map(str, sorted(c['ga'].astype(int))))})
res.append({'analysis': 'GA distribution (controls)', 'detail': 'weeks: ' + ', '.join(map(str, sorted(n['ga'].astype(int))))})
res.append({'analysis': 'Cases with GA <= 25 weeks', 'detail': '%d of %d cases; 0 controls in this range' % ((c['ga'] <= 25).sum(), len(c))})
res.append({'analysis': 'Controls with GA >= 28 weeks', 'detail': '%d of %d controls; %d cases in this range' % ((n['ga'] >= 28).sum(), len(n), (c['ga'] >= 28).sum())})

# ---------- 1) restricted (overlap) subgroups ----------
for thr in [26, 27, 28]:
    sub = m[m['ga'] >= thr]
    a, ci, nsub, ncase, ncc = auc_boot(sub, seed=100 + thr)
    res.append({'analysis': 'Restricted to GA >= %d weeks' % thr,
                'n': nsub, 'n_case': ncase, 'n_ctrl': ncc,
                'auc': None if np.isnan(a) else round(a, 3),
                'ci': '%.3f-%.3f' % ci if not np.isnan(ci[0]) else 'NA'})

# also restrict to GA <= 27 (case-dominant side)
sub = m[m['ga'] <= 27]
a, ci, nsub, ncase, ncc = auc_boot(sub, seed=200)
res.append({'analysis': 'Restricted to GA <= 27 weeks', 'n': nsub, 'n_case': ncase, 'n_ctrl': ncc,
            'auc': None if np.isnan(a) else round(a, 3),
            'ci': '%.3f-%.3f' % ci if not np.isnan(ci[0]) else 'NA (only %d controls)' % ncc})

# ---------- 2) GA-matched (nearest neighbour, caliper) ----------
def match(cases, ctrls, caliper):
    used = set(); pairs = []
    for _, cc in cases.sort_values('ga').iterrows():
        cand = ctrls[~ctrls.index.isin(used)].copy()
        cand['d'] = (cand['ga'] - cc['ga']).abs()
        cand = cand[cand['d'] <= caliper].sort_values('d')
        if len(cand):
            j = cand.index[0]; used.add(j)
            pairs.append((cc.name, j, cc['ga'], ctrls.loc[j, 'ga']))
    return pairs

for cal in [1, 2]:
    pairs = match(c, n, cal)
    if pairs:
        idx = [p[0] for p in pairs] + [p[1] for p in pairs]
        sub = m.loc[idx]
        a, ci, nsub, ncase, ncc = auc_boot(sub, seed=300 + cal)
        res.append({'analysis': 'GA-matched pairs (caliper %d week)' % cal,
                    'n': nsub, 'n_case': ncase, 'n_ctrl': ncc,
                    'auc': None if np.isnan(a) else round(a, 3),
                    'ci': '%.3f-%.3f' % ci if not np.isnan(ci[0]) else 'NA',
                    'detail': 'pairs=%d; matched GA pairs: %s' % (len(pairs), '; '.join('%d-%d' % (p[2], p[3]) for p in pairs))})
    else:
        res.append({'analysis': 'GA-matched pairs (caliper %d week)' % cal, 'n': 0,
                    'detail': 'no case-control pairs within caliper'})

# ---------- 3) BW-matched ----------
for cal in [150]:
    pairs = match(c.assign(ga=c['bw']), n.assign(ga=n['bw']), cal)
    if pairs:
        idx = [p[0] for p in pairs] + [p[1] for p in pairs]
        sub = m.loc[idx]
        a, ci, nsub, ncase, ncc = auc_boot(sub, seed=400)
        res.append({'analysis': 'BW-matched pairs (caliper %d g)' % cal,
                    'n': nsub, 'n_case': ncase, 'n_ctrl': ncc,
                    'auc': None if np.isnan(a) else round(a, 3),
                    'ci': '%.3f-%.3f' % ci if not np.isnan(ci[0]) else 'NA',
                    'detail': 'pairs=%d' % len(pairs)})

# ---------- 4) IPW-adjusted AUC ----------
Z = np.column_stack([m['ga'].values, m['bw'].values, (m['sex'] == 'male').astype(int).values])
ps_mod = LogisticRegression(C=np.inf, max_iter=8000).fit(Z, y)
ps = np.clip(ps_mod.predict_proba(Z)[:, 1], 1e-6, 1 - 1e-6)
w = np.where(y == 1, 1 / ps, 1 / (1 - ps))
gene_mod = LogisticRegression(C=np.inf, max_iter=8000).fit(m[G].values, y)
score = gene_mod.predict_proba(m[G].values)[:, 1]
auc_ipw = weighted_auc(y, score, w)
res.append({'analysis': 'IPW-adjusted AUC (weights from GA + BW + sex)',
            'n': len(m), 'auc': round(auc_ipw, 3),
            'detail': 'PS range %.3f-%.3f; cases with PS<0.1: %d; controls with PS>0.9: %d (limited overlap)'
                      % (ps.min(), ps.max(), int(((y == 1) & (ps < 0.1)).sum()), int(((y == 0) & (ps > 0.9)).sum()))})

out = pd.DataFrame(res)
out.to_csv(P("data/derived/GA_matching_stratified_results.csv"), index=False, encoding='utf-8-sig')
print(out.to_string(index=False))
print()
print('PS distribution: cases median %.3f (IQR %.3f-%.3f) | controls median %.3f (IQR %.3f-%.3f)'
      % (np.median(ps[y == 1]), *np.percentile(ps[y == 1], [25, 75]),
         np.median(ps[y == 0]), *np.percentile(ps[y == 0], [25, 75])))
