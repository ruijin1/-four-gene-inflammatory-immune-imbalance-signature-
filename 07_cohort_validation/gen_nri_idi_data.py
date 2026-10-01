# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/gen_nri_idi_data.py
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
Generate the summary and bootstrap detail for the incremental-value NRI/IDI analysis
Output:
  - NRI_IDI_summary.json          -- the numbers written to Table S18 and the manuscript
  - NRI_IDI_bootstrap_detail.csv  -- the NRI/IDI of every iteration (added as a new sheet to Data S6)
"""
import warnings
warnings.filterwarnings("ignore")
import sys, json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, P("lib"))
from delong_nri_idi import delong, nri_idi

RNG = np.random.default_rng(20260911)
B = 2000

TRAIN = P("data/derived/training_set_clinical_genotype.csv")
CLIN = P("data/derived/gse220135_samples.csv")
S9 = P("data/raw/supporting_information/S9 Data.xlsx")
OUT_JSON = P("data/derived/NRI_IDI_summary.json")
OUT_CSV = P("data/derived/NRI_IDI_bootstrap_detail.csv")

COEF = {"CD3G": -0.820, "PRTN3": 0.610, "CD6": -1.402, "THBS1": 1.560}
INTERCEPT = -3.375

# ---------- Training set ----------
tr = pd.read_csv(TRAIN)
tr["sexM"] = tr["sex"].astype(str).str.lower().str.startswith("m").astype(int)
ytr = tr["outcome"].values.astype(int)
lp_tr = INTERCEPT + sum(COEF[g] * tr[g].values for g in COEF)
Xc_tr = tr[["ga", "bw", "sexM"]].values.astype(float)

mC = LogisticRegression(C=np.inf, solver="lbfgs", max_iter=8000).fit(Xc_tr, ytr)
mD = LogisticRegression(C=np.inf, solver="lbfgs", max_iter=8000).fit(
    np.column_stack([Xc_tr, lp_tr]), ytr)
p_old_tr = mC.predict_proba(Xc_tr)[:, 1]
p_new_tr = mD.predict_proba(np.column_stack([Xc_tr, lp_tr]))[:, 1]

# ---------- External cohort ----------
clin = pd.read_csv(CLIN)
clin["gsm"] = clin["gsm"].astype(str)
clin["case"] = (clin["disease"] == "BPD").astype(int)
clin["sexM"] = (clin["sex"] == "M").astype(int)

ext = {}
for sheet, tag in (("Day14d", "day14"), ("Day28d", "day28")):
    e = pd.read_excel(S9, sheet_name=sheet)
    e["gsm"] = e["Sample"].astype(str).str.split("_").str[0]
    m = e.merge(clin[["gsm", "case", "ga", "bw", "sexM", "time"]], on="gsm", how="inner")
    m = m[m["time"] == tag.capitalize()].dropna(subset=["ga", "bw", "sexM"]).reset_index(drop=True)
    lp = INTERCEPT + sum(COEF[g] * m[g + "_vst"].values for g in COEF)
    Xe = m[["ga", "bw", "sexM"]].values.astype(float)
    ext[tag] = dict(y=m["case"].values.astype(int),
                    p_old=mC.predict_proba(Xe)[:, 1],
                    p_new=mD.predict_proba(np.column_stack([Xe, lp]))[:, 1])

# ---------- Computation + bootstrap ----------
datasets = [("training_day14", "Day-14 training set (GSE32472)", ytr, p_old_tr, p_new_tr)]
for tag in ("day14", "day28"):
    d = ext[tag]
    datasets.append(("external_%s" % tag, "GSE220135 %s" % tag, d["y"], d["p_old"], d["p_new"]))

summary = {}
detail_rows = []
for key, label, y, po, pn in datasets:
    y = np.asarray(y); po = np.asarray(po); pn = np.asarray(pn)
    n = len(y)
    auc_old = roc_auc_score(y, po); auc_new = roc_auc_score(y, pn)
    _, _, z, pv = delong(y, pn, po)
    nri, idi = nri_idi(y, po, pn)
    nris, idis = [], []
    for _ in range(B):
        idx = RNG.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        a, b = nri_idi(y[idx], po[idx], pn[idx])
        nris.append(a); idis.append(b)
        detail_rows.append(dict(dataset=label, iteration=len(nris), NRI=a, IDI=b))
    nris = np.asarray(nris); idis = np.asarray(idis)
    nri_lo, nri_hi = np.percentile(nris, [2.5, 97.5])
    idi_lo, idi_hi = np.percentile(idis, [2.5, 97.5])
    summary[key] = dict(
        label=label, n=int(n), n_case=int(y.sum()), n_ctrl=int((1 - y).sum()),
        auc_perinatal=round(float(auc_old), 4), auc_combined=round(float(auc_new), 4),
        delta_auc=round(float(auc_new - auc_old), 4),
        delong_z=round(float(z), 3), delong_p=round(float(pv), 4),
        nri=round(float(nri), 3), nri_lo=round(float(nri_lo), 3), nri_hi=round(float(nri_hi), 3),
        nri_frac_le0=round(float((nris <= 0).mean()), 3),
        idi=round(float(idi), 3), idi_lo=round(float(idi_lo), 3), idi_hi=round(float(idi_hi), 3),
        idi_frac_le0=round(float((idis <= 0).mean()), 3),
        n_boot=int(len(nris)))
    s = summary[key]
    print("%-42s n=%3d AUC %.4f -> %.4f (dAUC %+.4f, DeLong P=%.4f) | NRI %+.3f (%.3f, %.3f) | IDI %+.3f (%.3f, %.3f)"
          % (label, n, auc_old, auc_new, auc_new - auc_old, pv, nri, nri_lo, nri_hi, idi, idi_lo, idi_hi))

with open(OUT_JSON, "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)
pd.DataFrame(detail_rows).to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
print("\n[saved]", OUT_JSON)
print("[saved]", OUT_CSV, "(rows=%d)" % len(detail_rows))
