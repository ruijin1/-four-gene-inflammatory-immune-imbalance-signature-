# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/baseline_adjust.py
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

"""Table 1 (baseline) + GA/BW/sex-adjusted sensitivity analysis"""
import pandas as pd, numpy as np, re, json
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

OUT = {}

# ---------------- 1) GSE32472 ----------------
d72 = pd.read_csv('gse32472_samples.csv')
d72['ga'] = pd.to_numeric(d72['gestational age (weeks)'], errors='coerce')
d72['bw'] = pd.to_numeric(d72['birth weight (grams)'], errors='coerce')
d72['sex'] = d72['gender'].str.extract(r'(female|male)')
d72['days'] = pd.to_numeric(d72['days of life'], errors='coerce')
d72['bpd36'] = d72['oxygen 36 weeks postmenstrual age (wpma)'].str.startswith('1')
d72['ox28'] = d72['oxygen 28 days'].str.startswith('1')


def subset(df, tp):
    s = df[df['timepoint'] == tp].copy()
    s = s[s['bpd36'] | ((~s['bpd36']) & (~s['ox28']))].copy()
    s['case'] = s['bpd36'].astype(int)
    return s


t14 = subset(d72, '12-17')   # training
t28 = subset(d72, '24-33')   # temporal
assert len(t14) == 63 and len(t28) == 62, (len(t14), len(t28))

# ---------------- 2) GSE220135 ----------------
d135 = pd.read_csv('gse220135_samples.csv')
d135['case'] = (d135['disease'] == 'BPD').astype(int)
g14 = d135[d135['time'] == 'Day14'].copy()
g28 = d135[d135['time'] == 'Day28'].copy()


def fmt_iqr(x, dec=0):
    x = pd.Series(x).dropna()
    if len(x) == 0:
        return 'NA'
    f = ('%.' + str(dec) + 'f')
    return '%s (%s-%s)' % (f % x.median(), f % x.quantile(.25), f % x.quantile(.75))


def mwu(a, b):
    a, b = pd.Series(a).dropna(), pd.Series(b).dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return stats.mannwhitneyu(a, b, alternative='two-sided').pvalue


def sextest(case_sex, ctrl_sex):
    a = (pd.Series(case_sex).dropna() == 'male').sum(), (pd.Series(case_sex).dropna() == 'female').sum()
    b = (pd.Series(ctrl_sex).dropna() == 'male').sum(), (pd.Series(ctrl_sex).dropna() == 'female').sum()
    tab = np.array([a, b])
    if (tab.sum(axis=0) > 0).all() and (tab.sum(axis=1) > 0).all():
        exp = stats.chi2_contingency(tab)[3]
        if (exp < 5).any():
            return stats.fisher_exact(tab)[1], 'Fisher'
        return stats.chi2_contingency(tab)[0:1][0] and stats.chi2_contingency(tab)[1], 'chi2'
    return np.nan, 'NA'


def pfmt(p, dec=3):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return 'NA'
    return ('<0.001' if p < 0.001 else ('%.' + str(dec) + 'f') % p)


def group_stats(df, label, ga_dec=1, with_days=False):
    c = df[df['case'] == 1]
    n = df[df['case'] == 0]
    ncm, ncf = (c['sex'] == 'male').sum(), (c['sex'] == 'female').sum()
    nnm, nnf = (n['sex'] == 'male').sum(), (n['sex'] == 'female').sum()
    row = {
        'group': label,
        'n_case': int(len(c)), 'n_ctrl': int(len(n)),
        'case_n_cells': len(c), 'ctrl_n_cells': len(n),
        'case_ga': fmt_iqr(c['ga'], ga_dec), 'ctrl_ga': fmt_iqr(n['ga'], ga_dec),
        'ga_p': pfmt(mwu(c['ga'], n['ga'])),
        'case_bw': fmt_iqr(c['bw'], 0), 'ctrl_bw': fmt_iqr(n['bw'], 0),
        'bw_p': pfmt(mwu(c['bw'], n['bw'])),
        'case_male': '%d (%d%%)' % (ncm, round(100 * ncm / max(1, ncm + ncf))),
        'ctrl_male': '%d (%d%%)' % (nnm, round(100 * nnm / max(1, nnm + nnf))),
        'sex_p': pfmt(sextest(c['sex'], n['sex'])[0]),
        'sex_test': sextest(c['sex'], n['sex'])[1],
    }
    if with_days:
        row['case_days'] = fmt_iqr(c['days'], 0)
        row['ctrl_days'] = fmt_iqr(n['days'], 0)
        row['days_p'] = pfmt(mwu(c['days'], n['days']))
    return row


rows = [
    group_stats(t14, 'GSE32472 day 14 (training)', with_days=True),
    group_stats(t28, 'GSE32472 day 28 (temporal)', with_days=True),
    group_stats(g14, 'GSE220135 day 14 (external)'),
    group_stats(g28, 'GSE220135 day 28 (external)'),
]
print('=== Baseline characteristics ===')
for r in rows:
    print(r['group'], '| n %d vs %d' % (r['n_case'], r['n_ctrl']))
    print('   GA: %s vs %s  p=%s | BW: %s vs %s  p=%s | male: %s vs %s  p=%s (%s)'
          % (r['case_ga'], r['ctrl_ga'], r['ga_p'], r['case_bw'], r['ctrl_bw'], r['bw_p'],
             r['case_male'], r['ctrl_male'], r['sex_p'], r['sex_test']))
    if 'case_days' in r:
        print('   days of life: %s vs %s  p=%s' % (r['case_days'], r['ctrl_days'], r['days_p']))
OUT['baseline'] = rows

# ---------------- 3) adjusted analysis (training set) ----------------
EX = P("data/raw/GSE32472/GSE32472_series_matrix.txt/正常组与中重度BPD的比较/12-17d/炎症并集/Nomogram/预测值-系数.xlsx")
e = pd.read_excel(EX, sheet_name='glm模型')
assert list(e['outcome'].astype(int)) == list(t14['case']), 'order mismatch!'
G = ['CD3G', 'PRTN3', 'CD6', 'THBS1']
Xg = e[G].values
cli = t14[['ga', 'bw']].copy()
cli['sex_male'] = (t14['sex'] == 'male').astype(int)
Xc = cli.values.astype(float)
y = e['outcome'].astype(int).values


def fit_auc(X, y, seed=0):
    m = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs').fit(X, y)
    p = m.predict_proba(X)[:, 1]
    return m, roc_auc_score(y, p), p


m_base, auc_base, p_base = fit_auc(Xg, y)
m_adj, auc_adj, p_adj = fit_auc(np.hstack([Xg, Xc]), y)
m_cli, auc_cli, p_cli = fit_auc(Xc, y)

rng = np.random.default_rng(20260910)
bs_base, bs_adj, bs_delta = [], [], []
for _ in range(1000):
    idx = rng.integers(0, len(y), len(y))
    if len(np.unique(y[idx])) < 2:
        continue
    try:
        _, a1, _ = fit_auc(Xg[idx], y[idx])
        _, a2, _ = fit_auc(np.hstack([Xg, Xc])[idx], y[idx])
    except Exception:
        continue
    bs_base.append(a1); bs_adj.append(a2); bs_delta.append(a1 - a2)
bs_base, bs_adj, bs_delta = map(np.array, (bs_base, bs_adj, bs_delta))
ci_base = np.percentile(bs_base, [2.5, 97.5])
ci_adj = np.percentile(bs_adj, [2.5, 97.5])
ci_delta = np.percentile(bs_delta, [2.5, 97.5])

print()
print('=== Adjusted sensitivity analysis (day-14 training set, n=63) ===')
print('A  four-gene signature only      : AUC = %.4f (bootstrap 95%% CI %.3f-%.3f)' % (auc_base, *ci_base))
print('B  four-gene + GA + BW + sex      : AUC = %.4f (bootstrap 95%% CI %.3f-%.3f)' % (auc_adj, *ci_adj))
print('C  GA + BW + sex only             : AUC = %.4f' % auc_cli)
print('delta AUC (A-B) = %+.4f (95%% CI %+.4f to %+.4f)' % (auc_base - auc_adj, *ci_delta))
print()
print('Model B coefficients:')
names = G + ['gestational age (weeks)', 'birth weight (g)', 'male sex']
for n_, c_, in zip(names, m_adj.coef_[0]):
    print('   %-28s beta = %+.4f' % (n_, c_))
print('   intercept = %+.4f' % m_adj.intercept_[0])

# univariable AUC (each clinical variable separately)
for k, arr in [('gestational age', cli['ga'].values), ('birth weight', cli['bw'].values), ('male sex', cli['sex_male'].values)]:
    a = roc_auc_score(y, arr) if len(np.unique(arr)) > 1 else np.nan
    a = max(a, 1 - a)
    print('single-variable AUC %-18s = %.3f' % (k, a))

OUT['adjusted'] = {
    'auc_gene_only': auc_base, 'ci_gene_only': list(ci_base),
    'auc_gene_plus_clinical': auc_adj, 'ci_gene_plus_clinical': list(ci_adj),
    'auc_clinical_only': auc_cli,
    'delta_auc': auc_base - auc_adj, 'ci_delta': list(ci_delta),
    'coef_model_B': dict(zip(names, [float(v) for v in m_adj.coef_[0]])),
    'intercept_model_B': float(m_adj.intercept_[0]),
}
json.dump(OUT, open('baseline_and_adjusted.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=2)

# write the merged table used for modelling
merged = pd.DataFrame({'gsm': t14['gsm'].values, 'outcome': y,
                       'ga': cli['ga'].values, 'bw': cli['bw'].values, 'sex': t14['sex'].values})
for g in G:
    merged[g] = Xg[:, G.index(g)]
merged.to_csv('training_set_clinical_genotype.csv', index=False, encoding='utf-8-sig')
print()
print('saved: baseline_and_adjusted.json, training_set_clinical_genotype.csv')
