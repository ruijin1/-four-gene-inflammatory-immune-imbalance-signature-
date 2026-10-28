# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/external_clinical_comparison.py
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
External-cohort 'value' analysis (the fairest comparison design):
    * Model A: four-gene signature -- coefficients taken from the training set, fixed, applied directly to the external cohort
    * Model C: GA + BW + sex      -- coefficients likewise fitted on the TRAINING set, then applied to the external cohort
    * Model D: A + C combined     -- as above

All models share the same origin (training set), the same target (external cohort) and the
An optimism-corrected AUC (Harrell bootstrap) is also reported for a clinical model
refitted inside the external cohort, together with univariable AUCs and the sex distribution.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

RNG = np.random.default_rng(20260910)
B = 1000

TRAIN = P("data/derived/training_set_clinical_genotype.csv")
S9 = P("data/raw/supporting_information/S9 Data.xlsx")
CLIN = P("data/derived/gse220135_samples.csv")

COEF = {"CD3G": -0.820, "PRTN3": 0.610, "CD6": -1.402, "THBS1": 1.560}
INTERCEPT = -3.375


def auc_ci(y, s, n=B):
    y = np.asarray(y); s = np.asarray(s)
    a = roc_auc_score(y, s)
    bs = []
    for _ in range(n):
        idx = RNG.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        bs.append(roc_auc_score(y[idx], s[idx]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return a, float(lo), float(hi)


def fit_logit(X, y):
    lr = LogisticRegression(penalty=None, solver="lbfgs", max_iter=8000)
    lr.fit(X, y)
    return lr


def optimism_corrected_auc(X, y, labels=None, n=200):
    """Harrell optimism correction: fit on the original data, score on bootstrap samples, take the optimism."""
    lr = fit_logit(X, y)
    apparent = roc_auc_score(y, lr.predict_proba(X)[:, 1])
    opts = []
    for _ in range(n):
        idx = RNG.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        try:
            lb = fit_logit(X[idx], y[idx])
        except Exception:
            continue
        p_boot = lb.predict_proba(X[idx])[:, 1]
        p_orig = lb.predict_proba(X)[:, 1]
        if len(np.unique(y[idx])) < 2:
            continue
        opts.append(roc_auc_score(y[idx], p_boot) - roc_auc_score(y, p_orig))
    return apparent, apparent - float(np.mean(opts))


# ---------- Training set ----------
tr = pd.read_csv(TRAIN)
tr["sexM"] = tr["sex"].astype(str).str.lower().str.startswith("m").astype(int)
ytr = tr["outcome"].values
lp_tr = INTERCEPT + sum(COEF[g] * tr[g].values for g in COEF)

# Fit Model C (GA+BW+sex) and Model D (+ four genes) on the training set
mC_tr = fit_logit(tr[["ga", "bw", "sexM"]].values.astype(float), ytr)
Xd_tr = np.column_stack([tr[["ga", "bw", "sexM"]].values.astype(float), lp_tr])
mD_tr = fit_logit(Xd_tr, ytr)

print("=== Coefficients fitted on the training set (applied as-is to the external cohort) ===")
print("  Model C GA+BW+sex : intercept %+.3f | coef %s"
      % (mC_tr.intercept_[0], np.round(mC_tr.coef_[0], 4).tolist()))
print("  Model D GA+BW+sex+four genes : intercept %+.3f | coef %s"
      % (mD_tr.intercept_[0], np.round(mD_tr.coef_[0], 4).tolist()))
print()

# ---------- External cohort ----------
expr = {}
for sheet, tag in (("Day14d", "Day14"), ("Day28d", "Day28")):
    df = pd.read_excel(S9, sheet_name=sheet)
    df["gsm"] = df["Sample"].astype(str).str.split("_").str[0]
    expr[tag] = df

clin = pd.read_csv(CLIN)
clin["gsm"] = clin["gsm"].astype(str)
clin["case"] = (clin["disease"] == "BPD").astype(int)
clin["sexM"] = (clin["sex"] == "M").astype(int)

print("=== Sex distribution in the external cohort (sex availability to be verified) ===")
for tag in ("Day14", "Day28"):
    sub = clin[clin["time"] == tag]
    c1 = sub[sub.case == 1]; c0 = sub[sub.case == 0]
    print("  %s: cases male %d/female %d (missing %d) | controls male %d/female %d (missing %d)"
          % (tag, (c1.sex == "M").sum(), (c1.sex == "F").sum(), c1.sex.isna().sum(),
             (c0.sex == "M").sum(), (c0.sex == "F").sum(), c0.sex.isna().sum()))
print()

rows = []
for tag in ("Day14", "Day28"):
    m = expr[tag].merge(clin[["gsm", "case", "ga", "bw", "sexM", "time"]], on="gsm", how="inner")
    m = m[m["time"] == tag].copy()
    y = m["case"].values
    lp_ext = INTERCEPT + sum(COEF[g] * m[g + "_vst"].values for g in COEF)
    # Models C/D need complete clinical variables (sex has a few missing values) -> complete-case subset
    mc = m.dropna(subset=["ga", "bw", "sexM"]).copy()
    yc = mc["case"].values
    lp_c = INTERCEPT + sum(COEF[g] * mc[g + "_vst"].values for g in COEF)
    Xc = mc[["ga", "bw", "sexM"]].values.astype(float)
    Xd = np.column_stack([Xc, lp_c])

    print("=" * 80)
    print("### GSE220135 %s | n=%d (cases %d / controls %d)" % (tag, len(m), y.sum(), (1 - y).sum()))

    a, lo, hi = auc_ci(y, lp_ext)
    print("  A  Four-gene signature (fixed training-set coefficients)          AUC %.3f (95%% CI %.3f-%.3f)  [n=%d]"
          % (a, lo, hi, len(m)))
    rows.append(("GSE220135 %s" % tag, "A Four-gene signature (fixed coefficients)", a, lo, hi, len(m)))

    aC, loC, hiC = auc_ci(yc, mC_tr.predict_proba(Xc)[:, 1])
    print("  C  GA+BW+sex (training-set coefficients, extrapolated as-is)      AUC %.3f (95%% CI %.3f-%.3f)  [n=%d]"
          % (aC, loC, hiC, len(mc)))
    rows.append(("GSE220135 %s" % tag, "C GA+BW+sex (training-set coefficients, extrapolated)", aC, loC, hiC, len(mc)))

    aD, loD, hiD = auc_ci(yc, mD_tr.predict_proba(Xd)[:, 1])
    print("  D  GA+BW+sex+signature (training-set coefficients, extrapolated as-is) AUC %.3f (95%% CI %.3f-%.3f)  [n=%d]"
          % (aD, loD, hiD, len(mc)))
    rows.append(("GSE220135 %s" % tag, "D A+C (training-set coefficients, extrapolated)", aD, loD, hiD, len(mc)))
    print()

    # Clinical model refitted inside the external cohort (apparent vs optimism-corrected)
    app, corr = optimism_corrected_auc(Xc, yc)
    print("  C' GA+BW+sex (refitted inside the external cohort)       apparent AUC %.3f -> optimism-corrected %.3f  [n=%d]"
          % (app, corr, len(mc)))
    rows.append(("GSE220135 %s" % tag, "C' Refitted inside the external cohort (corrected)", corr, np.nan, np.nan, len(mc)))
    print()

    for c, lab, flip in (("ga", "GA only", True), ("bw", "BW only", True),
                         ("sexM", "male sex only", False)):
        s = m[c].dropna().values.astype(float)
        yy = m.loc[m[c].notna(), "case"].values
        if len(np.unique(s)) < 2:
            print("  -  %-14s not estimable (single value)" % lab)
            continue
        a3, lo3, hi3 = auc_ci(yy, s)
        if flip:                       # low GA / low BW = high risk -> orientation flipped
            a3, lo3, hi3 = 1 - a3, 1 - hi3, 1 - lo3
            lab = lab.replace("only", "only (risk-oriented)")
        print("  -  %-26s AUC %.3f (95%% CI %.3f-%.3f)" % (lab, a3, lo3, hi3))
        rows.append(("GSE220135 %s" % tag, lab, a3, lo3, hi3, len(s)))
    print()

# Training set, same definitions (reference)
print("=" * 80)
print("### Reference: day-14 training set GSE32472")
a, lo, hi = auc_ci(ytr, lp_tr)
print("  A  Four-gene signature (within the same set)              AUC %.3f (95%% CI %.3f-%.3f)" % (a, lo, hi))
Xc = tr[["ga", "bw", "sexM"]].values.astype(float)
aC, loC, hiC = auc_ci(ytr, mC_tr.predict_proba(Xc)[:, 1])
print("  C  GA+BW+sex (within the same set)              AUC %.3f (95%% CI %.3f-%.3f)" % (aC, loC, hiC))
Xd = np.column_stack([Xc, lp_tr])
aD, loD, hiD = auc_ci(ytr, mD_tr.predict_proba(Xd)[:, 1])
print("  D  Combined (within the same set)                   AUC %.3f (95%% CI %.3f-%.3f)" % (aD, loD, hiD))

out = pd.DataFrame(rows, columns=["cohort", "model", "AUC", "CI_low", "CI_high", "n"])
out.to_csv(P("data/derived/External_clinical_comparison.csv"),
           index=False, encoding="utf-8-sig")
print("\n[saved] data/derived/External_clinical_comparison.csv")
