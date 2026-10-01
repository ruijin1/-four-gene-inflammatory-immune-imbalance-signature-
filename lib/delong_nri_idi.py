# -*- coding: utf-8 -*-
"""
Reusable statistics utilities: DeLong test / continuous NRI / IDI / VIF / all-subset AUC
Depends only on numpy, pandas, scikit-learn and scipy.

As a module:
    from delong_nri_idi import delong, nri_idi, vif_table, all_subsets_auc, bootstrap_ci

As a script: fill in CONFIG below and run directly.
"""
import numpy as np
import pandas as pd
from itertools import combinations
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_auc_score
from scipy import stats

RNG = np.random.default_rng(20260911)


# ----------------------------------------------------------------------
# DeLong test (Sun & Xu 2014 fast implementation)
# ----------------------------------------------------------------------
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
    """Compare the AUCs of two prediction scores. Returns (auc1, auc2, z, p). p2 is the reference (old) model."""
    y = np.asarray(y); p1 = np.asarray(p1); p2 = np.asarray(p2)
    order = np.argsort(-y)
    m = int(y.sum())
    preds = np.vstack((p1, p2))[:, order]
    aucs, cov = _fast_delong(preds, m)
    L = np.array([[1.0, -1.0]])
    # numpy 2.x: .item() is required, otherwise float() raises "only 0-dimensional arrays..."
    var = float((L @ cov @ L.T).item())
    if var <= 0:
        return aucs[0], aucs[1], np.nan, np.nan
    z = float((L @ aucs).item()) / np.sqrt(var)
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return aucs[0], aucs[1], z, p


# ----------------------------------------------------------------------
# Continuous NRI / IDI
# ----------------------------------------------------------------------
def nri_idi(y, p_old, p_new):
    """Continuous NRI and IDI. Values > 0 mean the new model is better."""
    y = np.asarray(y); ev = y == 1; ne = ~ev
    up_ev = float((p_new[ev] > p_old[ev]).mean()); dn_ev = float((p_new[ev] < p_old[ev]).mean())
    up_ne = float((p_new[ne] > p_old[ne]).mean()); dn_ne = float((p_new[ne] < p_old[ne]).mean())
    nri = (up_ev - dn_ev) + (dn_ne - up_ne)
    idi = float((p_new[ev].mean() - p_old[ev].mean()) - (p_new[ne].mean() - p_old[ne].mean()))
    return nri, idi


def bootstrap_ci(fn, y, n=1000, rng=None):
    """Bootstrap fn(idx) and return the (lo, hi) percentile CI. fn must return a scalar."""
    rng = rng or RNG
    N = len(y)
    vals = []
    for _ in range(n):
        idx = rng.integers(0, N, N)
        if len(np.unique(np.asarray(y)[idx])) < 2:
            continue
        try:
            vals.append(fn(idx))
        except Exception:
            continue
    if not vals:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def bootstrap_auc_ci(y, score, n=2000, rng=None):
    rng = rng or RNG
    y = np.asarray(y); score = np.asarray(score)
    bs = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        bs.append(roc_auc_score(y[idx], score[idx]))
    return roc_auc_score(y, score), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


# ----------------------------------------------------------------------
# Collinearity
# ----------------------------------------------------------------------
def vif_table(df, cols):
    """Returns {column: VIF}. VIF > 5-10 is taken as collinear."""
    out = {}
    for c in cols:
        others = [x for x in cols if x != c]
        if not others:
            out[c] = 1.0
            continue
        X = df[others].values.astype(float)
        yv = df[c].values.astype(float)
        r2 = LinearRegression().fit(X, yv).score(X, yv)
        out[c] = float(1 / (1 - r2)) if r2 < 1 else np.inf
    return out


def collinearity_report(df, cols):
    """Print pairwise Pearson r and VIF."""
    for a, b in combinations(cols, 2):
        sub = df[[a, b]].dropna()
        r, p = stats.pearsonr(sub[a], sub[b])
        print("  r(%s,%s) = %.3f (P=%.3g)" % (a, b, r, p))
    v = vif_table(df[cols].dropna(), cols)
    print("  VIF: " + "  ".join("%s=%.2f" % (k, val) for k, val in v.items()))
    return v


# ----------------------------------------------------------------------
# All-subset AUC
# ----------------------------------------------------------------------
def all_subsets_auc(df, clin_cols, y_col, fit_df=None, y_fit=None, predict_df=None,
                    seed=20260911):
    """
    Enumerate all non-empty subsets of the clinical variables and compute their AUC.
    By default fit_df = df: coefficients are fitted on fit_df and evaluated on predict_df (default df).
    For an external cohort, pass fit_df = training set and predict_df = external cohort.
    Returns DataFrame[subset, n_vars, AUC, n].
    """
    fit_df = df if fit_df is None else fit_df
    y_fit = fit_df[y_col].values.astype(int) if y_fit is None else np.asarray(y_fit)
    predict_df = fit_df if predict_df is None else predict_df
    y_pr = predict_df[y_col].values.astype(int)
    rows = []
    for k in range(1, len(clin_cols) + 1):
        for s in combinations(clin_cols, k):
            mdl = LogisticRegression(C=np.inf, solver="lbfgs", max_iter=8000)
            mdl.fit(fit_df[list(s)].values.astype(float), y_fit)
            auc = roc_auc_score(y_pr, mdl.predict_proba(
                predict_df[list(s)].values.astype(float))[:, 1])
            rows.append(("+".join(s), k, auc, len(y_pr)))
    return pd.DataFrame(rows, columns=["subset", "n_vars", "AUC", "n"]).sort_values(
        "AUC", ascending=False).reset_index(drop=True)


# ----------------------------------------------------------------------
# Combined incremental-value report
# ----------------------------------------------------------------------
def incremental_value_report(y, p_clin, p_comb, label="", n_boot=1000):
    """Print AUC / dAUC / DeLong / continuous NRI / IDI with bootstrap CIs."""
    y = np.asarray(y); p_clin = np.asarray(p_clin); p_comb = np.asarray(p_comb)
    a0 = roc_auc_score(y, p_clin); a1 = roc_auc_score(y, p_comb)
    _, _, z, pv = delong(y, p_comb, p_clin)
    N_ = len(y)
    def _f(i, which):
        return lambda idx: nri_idi(y[idx], p_clin[idx], p_comb[idx])[which]
    lo_n, hi_n = bootstrap_ci(_f(None, 0), y, n_boot)
    lo_i, hi_i = bootstrap_ci(_f(None, 1), y, n_boot)
    nri, idi = nri_idi(y, p_clin, p_comb)
    print("  %s" % label)
    print("    clinical model AUC = %.4f -> + marker AUC = %.4f  (dAUC %+.4f)" % (a0, a1, a1 - a0))
    print("    DeLong: z = %.3f, P = %.4f" % (z, pv))
    print("    continuous NRI = %+.3f (95%%CI %+.3f to %+.3f)" % (nri, lo_n, hi_n))
    print("    IDI      = %+.3f (95%%CI %+.3f to %+.3f)" % (idi, lo_i, hi_i))
    return dict(auc_clin=a0, auc_comb=a1, d_auc=a1 - a0, delong_p=pv,
                nri=nri, nri_ci=(lo_n, hi_n), idi=idi, idi_ci=(lo_i, hi_i))


# ----------------------------------------------------------------------
# CONFIG -- fill in to run as a script
# ----------------------------------------------------------------------
CONFIG = {
    "train_csv": None,          # path to the training-set csv
    "y_col": "outcome",         # outcome column (0/1)
    "clin_cols": ["ga", "bw", "sexM"],   # clinical variable columns
    "ext": None,                # {"day14": {"csv": path, "y_col": ...}, ...}
}

if __name__ == "__main__":
    if not CONFIG["train_csv"]:
        print(__doc__)
        print("CONFIG is not filled in; import as a module instead:")
        print("  from delong_nri_idi import delong, nri_idi, vif_table, all_subsets_auc")
    else:
        tr = pd.read_csv(CONFIG["train_csv"])
        print("=== Collinearity ===")
        collinearity_report(tr, CONFIG["clin_cols"])
        print("\n=== Training set, all-subset AUC ===")
        print(all_subsets_auc(tr, CONFIG["clin_cols"], CONFIG["y_col"]).to_string(index=False))
