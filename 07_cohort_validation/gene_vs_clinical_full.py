# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/gene_vs_clinical_full.py
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

"""
Four-gene signature vs routine perinatal variables: systematic assessment of
Output:
    * GA-BW collinearity (Pearson r, VIF)
    * AUC of every subset of clinical variables in the training and external cohorts
    * continuous NRI + IDI: perinatal model -> + four genes
    * DeLong test for delta-AUC
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from itertools import combinations
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_auc_score
from scipy import stats

RNG = np.random.default_rng(20260911)
B = 2000

TRAIN = P("data/derived/training_set_clinical_genotype.csv")
CLIN = P("data/derived/gse220135_samples.csv")
S9 = P("data/raw/supporting_information/S9 Data.xlsx")

COEF = {"CD3G": -0.820, "PRTN3": 0.610, "CD6": -1.402, "THBS1": 1.560}
INTERCEPT = -3.375
CLLAB = {"ga": "GA", "bw": "BW", "sexM": "Sex"}


# ---------------- DeLong ----------------
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
    pos = pred_sorted[:, :m]
    neg = pred_sorted[:, m:]
    tx = np.empty([k, m]); ty = np.empty([k, n]); tz = np.empty([k, m + n])
    for r in range(k):
        tx[r, :] = _midrank(pos[r, :])
        ty[r, :] = _midrank(neg[r, :])
        tz[r, :] = _midrank(pred_sorted[r, :])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    cov = sx / m + sy / n
    return aucs, cov


def delong(y, p1, p2):
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


# ---------------- NRI / IDI ----------------
def nri_idi(y, p_old, p_new):
    y = np.asarray(y); ev = y == 1; ne = ~ev
    up_ev = float((p_new[ev] > p_old[ev]).mean()); dn_ev = float((p_new[ev] < p_old[ev]).mean())
    up_ne = float((p_new[ne] > p_old[ne]).mean()); dn_ne = float((p_new[ne] < p_old[ne]).mean())
    nri = (up_ev - dn_ev) + (dn_ne - up_ne)
    idi = float((p_new[ev].mean() - p_old[ev].mean()) - (p_new[ne].mean() - p_old[ne].mean()))
    return nri, idi


def boot_ci(fn, n=1000):
    vals = []
    for _ in range(n):
        idx = RNG.integers(0, N_, N_)
        if len(np.unique(Y_[idx])) < 2:
            continue
        try:
            vals.append(fn(idx))
        except Exception:
            continue
    if not vals:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def vif_table(df, cols):
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


# ---------------- Data ----------------
tr = pd.read_csv(TRAIN)
print("TRAIN columns:", list(tr.columns))
tr["sexM"] = tr["sex"].astype(str).str.lower().str.startswith("m").astype(int)
ytr = tr["outcome"].values.astype(int)
lp_tr = INTERCEPT + sum(COEF[g] * tr[g].values for g in COEF)

clin = pd.read_csv(CLIN)
clin["gsm"] = clin["gsm"].astype(str)
clin["case"] = (clin["disease"] == "BPD").astype(int)
clin["sexM"] = (clin["sex"] == "M").astype(int)

ext = {}
for sheet, tag in (("Day14d", "Day14"), ("Day28d", "Day28")):
    e = pd.read_excel(S9, sheet_name=sheet)
    e["gsm"] = e["Sample"].astype(str).str.split("_").str[0]
    m = e.merge(clin[["gsm", "case", "ga", "bw", "sexM", "time"]], on="gsm", how="inner")
    m = m[m["time"] == tag].copy()
    m["lp_gene"] = INTERCEPT + sum(COEF[g] * m[g + "_vst"].values for g in COEF)
    ext[tag] = m

lines = []


def P(s=""):
    print(s)
    lines.append(s)


print("=" * 92)
print("PART 1  Collinearity of perinatal variables (do GA and BW overlap heavily? -> is a parsimonious model viable)")
print("=" * 92)
for name, df in [("Training set GSE32472 day14", tr),
                 ("GSE220135 day14", ext["Day14"].dropna(subset=["ga", "bw", "sexM"])),
                 ("GSE220135 day28", ext["Day28"].dropna(subset=["ga", "bw", "sexM"]))]:
    sub = df[["ga", "bw", "sexM"]].dropna()
    r, p = stats.pearsonr(sub["ga"], sub["bw"])
    v = vif_table(sub, ["ga", "bw", "sexM"])
    print("%-24s n=%3d | Pearson r(GA,BW)=%.3f (P=%.3g) | VIF GA=%.2f  BW=%.2f  Sex=%.2f"
      % (name, len(sub), r, p, v["ga"], v["bw"], v["sexM"]))
print("Note: VIF > 5 is taken as collinear. If GA and BW are collinear, dropping one yields a pre-specified parsimonious model, but both models must be reported.")

# ---------------- Training set, all subsets ----------------
print("")
print("=" * 92)
print("PART 2  Training set GSE32472 day14 (n=63): every combination of clinical variables vs the four-gene signature")
print("=" * 92)
subsets = []
for k in (1, 2, 3):
    subsets += list(combinations(["ga", "bw", "sexM"], k))
rows2 = []
for s in subsets:
    X = tr[list(s)].values.astype(float)
    mdl = LogisticRegression(penalty=None, solver="lbfgs", max_iter=8000).fit(X, ytr)
    a = roc_auc_score(ytr, mdl.predict_proba(X)[:, 1])
    lab = "+".join(CLLAB[c] for c in s)
    rows2.append((lab, a, len(s)))
rows2.sort(key=lambda t: -t[1])
for lab, a, k in rows2:
    flag = ""
    if a < 0.9284188:
        flag = "   <-- below the four-gene signature (0.928)"
    print("  clinical %-14s (%d vars)  AUC %.4f%s" % (lab, k, a, flag))
print("  %-20s        AUC %.4f" % ("four-gene signature (4 vars)", 0.9284))

# ---------------- External cohort ----------------
print("")
print("=" * 92)
print("PART 3  External cohort GSE220135: coefficients of the clinical subset models are fitted on the TRAINING set and extrapolated as-is (not re-estimated)")
print("=" * 92)
ext_rows = []
for tag in ("Day14", "Day28"):
    m = ext[tag].dropna(subset=["ga", "bw", "sexM"]).copy()
    y = m["case"].values.astype(int)
    print("")
    print("--- GSE220135 %s | n=%d (cases %d / controls %d) ---" % (tag, len(m), y.sum(), (1 - y).sum()))
    a, lo, hi = roc_auc_score(y, m["lp_gene"].values), 0, 0
    bs = [roc_auc_score(y[i], m["lp_gene"].values[i]) for i in
          (RNG.integers(0, len(y), len(y)) for _ in range(B)) if len(np.unique(y[i])) > 1]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    print("  four-gene signature (training-set coefficients)            AUC %.3f (95%%CI %.3f-%.3f)" % (a, lo, hi))
    ext_rows.append((tag, "Four-gene signature", a, lo, hi, len(m)))
    for s in subsets:
        Xtr = tr[list(s)].values.astype(float)
        mdl = LogisticRegression(penalty=None, solver="lbfgs", max_iter=8000).fit(Xtr, ytr)
        Xe = m[list(s)].values.astype(float)
        pe = mdl.predict_proba(Xe)[:, 1]
        aa = roc_auc_score(y, pe)
        bs = []
        for _ in range(B):
            i = RNG.integers(0, len(y), len(y))
            if len(np.unique(y[i])) < 2:
                continue
            bs.append(roc_auc_score(y[i], pe[i]))
        lo2, hi2 = np.percentile(bs, [2.5, 97.5])
        lab = "+".join(CLLAB[c] for c in s)
        mark = "  <-- signature higher" if a > aa else ""
        print("  clinical %-14s (%d vars)        AUC %.3f (95%%CI %.3f-%.3f)%s" % (lab, len(s), aa, lo2, hi2, mark))
        ext_rows.append((tag, "Clinical " + lab, aa, lo2, hi2, len(m)))

# ---------------- NRI / IDI ----------------
print("")
print("=" * 92)
print("PART 4  Incremental value: adding the four genes on top of the perinatal model (continuous NRI / IDI)")
print("=" * 92)
print("  NRI > 0 with IDI > 0 indicates that the four genes add discriminatory information (the standard evidence of 'value' expected by reviewers)")
print("")
# Training set
Xc = tr[["ga", "bw", "sexM"]].values.astype(float)
mC = LogisticRegression(penalty=None, solver="lbfgs", max_iter=8000).fit(Xc, ytr)
p_old = mC.predict_proba(Xc)[:, 1]
mD = LogisticRegression(penalty=None, solver="lbfgs", max_iter=8000).fit(
    np.column_stack([Xc, lp_tr]), ytr)
p_new = mD.predict_proba(np.column_stack([Xc, lp_tr]))[:, 1]
N_, Y_ = len(ytr), ytr
nri, idi = nri_idi(ytr, p_old, p_new)
lo_n, hi_n = boot_ci(lambda i: nri_idi(ytr[i], p_old[i], p_new[i])[0])
lo_i, hi_i = boot_ci(lambda i: nri_idi(ytr[i], p_old[i], p_new[i])[1])
auc_old = roc_auc_score(ytr, p_old); auc_new = roc_auc_score(ytr, p_new)
_, _, z, pv = delong(ytr, p_new, p_old)
print("  Training set day14 (n=63)")
print("    perinatal model AUC = %.4f  ->  + four genes AUC = %.4f   dAUC = %+.4f" % (auc_old, auc_new, auc_new - auc_old))
print("    DeLong delta AUC: z = %.3f, P = %.4f" % (z, pv))
print("    continuous NRI = %+.3f (95%%CI %+.3f to %+.3f)" % (nri, lo_n, hi_n))
print("    IDI      = %+.3f (95%%CI %+.3f to %+.3f)" % (idi, lo_i, hi_i))
print("")
for tag in ("Day14", "Day28"):
    m = ext[tag].dropna(subset=["ga", "bw", "sexM"]).reset_index(drop=True)
    y = m["case"].values.astype(int)
    Xe = m[["ga", "bw", "sexM"]].values.astype(float)
    lp_e = m["lp_gene"].values
    po = mC.predict_proba(Xe)[:, 1]
    pn = mD.predict_proba(np.column_stack([Xe, lp_e]))[:, 1]
    N_, Y_ = len(y), y
    n2, i2 = nri_idi(y, po, pn)
    lo_n2, hi_n2 = boot_ci(lambda i: nri_idi(y[i], po[i], pn[i])[0], 1000)
    lo_i2, hi_i2 = boot_ci(lambda i: nri_idi(y[i], po[i], pn[i])[1], 1000)
    aold = roc_auc_score(y, po); anew = roc_auc_score(y, pn)
    _, _, z2, pv2 = delong(y, pn, po)
    print("  %s (n=%d)" % (tag, len(y)))
    print("    perinatal model AUC = %.4f  ->  + four genes AUC = %.4f   dAUC = %+.4f" % (aold, anew, anew - aold))
    print("    DeLong delta AUC: z = %.3f, P = %.4f" % (z2, pv2))
    print("    continuous NRI = %+.3f (95%%CI %+.3f to %+.3f)" % (n2, lo_n2, hi_n2))
    print("    IDI      = %+.3f (95%%CI %+.3f to %+.3f)" % (i2, lo_i2, hi_i2))
    print("")

# Save
pd.DataFrame(rows2, columns=["clinical_subset", "AUC_train", "n_vars"]).to_csv(
    P("data/derived/Clinical_subset_AUC_train.csv"), index=False, encoding="utf-8-sig")
pd.DataFrame(ext_rows, columns=["cohort", "model", "AUC", "CI_low", "CI_high", "n"]).to_csv(
    P("data/derived/Clinical_subset_AUC_external.csv"), index=False, encoding="utf-8-sig")
with open(P("data/derived/gene_vs_perinatal_variables_20260911.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print("\n[saved] data/derived/Clinical_subset_AUC_train.csv / _external.csv / gene_vs_perinatal_variables_20260911.txt")
