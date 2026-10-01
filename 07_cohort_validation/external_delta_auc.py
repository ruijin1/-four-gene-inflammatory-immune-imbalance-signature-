# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/external_delta_auc.py
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

"""External cohort: dAUC from adding the four genes to clinical variables only (bootstrap percentile CI)."""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

RNG = np.random.default_rng(20260910)
B = 2000

TRAIN = P("data/derived/training_set_clinical_genotype.csv")
S9 = P("data/raw/supporting_information/S9 Data.xlsx")
CLIN = P("data/derived/gse220135_samples.csv")
COEF = {"CD3G": -0.820, "PRTN3": 0.610, "CD6": -1.402, "THBS1": 1.560}
INTERCEPT = -3.375

tr = pd.read_csv(TRAIN)
tr["sexM"] = tr["sex"].astype(str).str.lower().str.startswith("m").astype(int)
ytr = tr["outcome"].values
lp_tr = INTERCEPT + sum(COEF[g] * tr[g].values for g in COEF)
Xtr = tr[["ga", "bw", "sexM"]].values.astype(float)

mC = LogisticRegression(penalty=None, solver="lbfgs", max_iter=8000).fit(Xtr, ytr)
mD = LogisticRegression(penalty=None, solver="lbfgs", max_iter=8000).fit(
    np.column_stack([Xtr, lp_tr]), ytr)

clin = pd.read_csv(CLIN)
clin["gsm"] = clin["gsm"].astype(str)
clin["case"] = (clin["disease"] == "BPD").astype(int)
clin["sexM"] = (clin["sex"] == "M").astype(int)

for sheet, tag in (("Day14d", "Day14"), ("Day28d", "Day28")):
    e = pd.read_excel(S9, sheet_name=sheet)
    e["gsm"] = e["Sample"].astype(str).str.split("_").str[0]
    m = e.merge(clin[["gsm", "case", "ga", "bw", "sexM", "time"]], on="gsm", how="inner")
    m = m[m["time"] == tag].dropna(subset=["ga", "bw", "sexM"]).copy()
    y = m["case"].values
    lp = INTERCEPT + sum(COEF[g] * m[g + "_vst"].values for g in COEF)
    Xc = m[["ga", "bw", "sexM"]].values.astype(float)
    pc = mC.predict_proba(Xc)[:, 1]
    pd_ = mD.predict_proba(np.column_stack([Xc, lp]))[:, 1]
    obs = roc_auc_score(y, pd_) - roc_auc_score(y, pc)
    d = []
    for _ in range(B):
        idx = RNG.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        d.append(roc_auc_score(y[idx], pd_[idx]) - roc_auc_score(y[idx], pc[idx]))
    d = np.asarray(d)
    lo, hi = np.percentile(d, [2.5, 97.5])
    print("%s | n=%d (cases %d)  dAUC (combined - clinical only) = %+.3f (95%% CI %+.3f to %+.3f) | bootstrap <= 0 proportion %.3f"
          % (tag, len(y), y.sum(), obs, lo, hi, float((d <= 0).mean())))
