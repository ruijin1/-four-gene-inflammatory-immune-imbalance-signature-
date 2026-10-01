# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/day28_increment_repro.py
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

import warnings; warnings.filterwarnings('ignore')
import sys, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
sys.path.insert(0, P("lib"))
from delong_nri_idi import delong, nri_idi
RNG = np.random.default_rng(20260911); B = 2000
GENES = ['CD3G','PRTN3','CD6','THBS1']
COEF = {'CD3G':-0.820,'PRTN3':0.610,'CD6':-1.402,'THBS1':1.560}; INTERCEPT = -3.375

tr = pd.read_csv('training_set_clinical_genotype.csv')
tr['sexM'] = tr['sex'].astype(str).str.lower().str.startswith('m').astype(int)
ytr = tr['outcome'].values.astype(int)
lp_tr = INTERCEPT + sum(COEF[g]*tr[g].values for g in GENES)
Xc_tr = tr[['ga','bw','sexM']].values.astype(float)
mC = LogisticRegression(C=np.inf, solver='lbfgs', max_iter=8000).fit(Xc_tr, ytr)
mD = LogisticRegression(C=np.inf, solver='lbfgs', max_iter=8000).fit(np.column_stack([Xc_tr, lp_tr]), ytr)
po = mC.predict_proba(Xc_tr)[:,1]; pn = mD.predict_proba(np.column_stack([Xc_tr, lp_tr]))[:,1]
a_old, a_new, z, pv = delong(ytr, pn, po); nri, idi = nri_idi(ytr, po, pn)
nris=[]; idis=[]
for _ in range(B):
    idx = RNG.integers(0, len(ytr), len(ytr))
    a,b = nri_idi(ytr[idx], po[idx], pn[idx]); nris.append(a); idis.append(b)
nlo,nhi = np.percentile(nris,[2.5,97.5]); ilo,ihi = np.percentile(idis,[2.5,97.5])
print('=== Reproduction reference: day-14 training set (AF13 records +0.028 / NRI +0.801 (0.303, 1.252) / IDI +0.073 (0.016, 0.134)) ===')
print('  n=%d  AUC perinatal %.4f  combined %.4f  dAUC %+.4f  DeLong P=%.4f' % (len(ytr), a_old, a_new, a_new-a_old, pv))
print('  NRI %+.3f (%.3f, %.3f)   IDI %+.3f (%.3f, %.3f)' % (nri,nlo,nhi,idi,ilo,ihi))
