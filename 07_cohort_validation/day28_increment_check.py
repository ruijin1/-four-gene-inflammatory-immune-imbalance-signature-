# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/day28_increment_check.py
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
import sys, json, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
sys.path.insert(0, P("lib"))
from delong_nri_idi import delong, nri_idi
RNG = np.random.default_rng(20260911); B = 2000
GENES = ['CD3G','PRTN3','CD6','THBS1']
COEF = {'CD3G':-0.820,'PRTN3':0.610,'CD6':-1.402,'THBS1':1.560}; INTERCEPT = -3.375
ROOT = (P("data/raw/GSE32472/GSE32472_series_matrix.txt"))

tr = pd.read_csv('training_set_clinical_genotype.csv')
tr['sexM'] = tr['sex'].astype(str).str.lower().str.startswith('m').astype(int)
ytr = tr['outcome'].values.astype(int)
lp_tr = INTERCEPT + sum(COEF[g]*tr[g].values for g in GENES)
Xc_tr = tr[['ga','bw','sexM']].values.astype(float)
mC = LogisticRegression(C=np.inf, solver='lbfgs', max_iter=8000).fit(Xc_tr, ytr)
mD = LogisticRegression(C=np.inf, solver='lbfgs', max_iter=8000).fit(np.column_stack([Xc_tr, lp_tr]), ytr)

mat = pd.read_excel(ROOT + r'/正常组与中重度BPD的比较/22-27d/GSE32472_22-27d.xlsx',
                    sheet_name='GSE32472_22-27dcomplete', index_col=0)
mat.index = [str(x).strip() for x in mat.index]
wb = pd.read_excel(P("data/raw/GSE32472.xlsx"), header=None, sheet_name='GSE32472_28d')
grade = pd.Series([int(str(v).split('(')[0].strip()) for v in wb.iloc[1,1:] if pd.notna(v)],
                  index=[str(v).strip() for v in wb.iloc[0,1:] if pd.notna(v)])
clin = pd.read_csv('gse32472_samples.csv'); clin['gsm']=clin['gsm'].astype(str).str.strip(); clin=clin.set_index('gsm')
clin['ga']=pd.to_numeric(clin['gestational age (weeks)'],errors='coerce')
clin['bw']=pd.to_numeric(clin['birth weight (grams)'],errors='coerce')
clin['sexM']=clin['gender'].astype(str).str.strip().str.startswith('2').astype(int)
cols=[c for c in grade.index if c in mat.columns]
e = mat.loc[GENES, cols].astype(float).T.join(clin[['ga','bw','sexM']], how='inner')
e['grade'] = grade.loc[e.index].values

y = np.where(e.grade.values==0, 0, np.where(e.grade.values>=2, 1, np.nan))
msk = ~np.isnan(y); yy = y[msk].astype(int)
lp = INTERCEPT + e[GENES].values[msk] @ np.array([COEF[g] for g in GENES])
Xe = e[['ga','bw','sexM']].values[msk].astype(float)
po = mC.predict_proba(Xe)[:,1]; pn = mD.predict_proba(np.column_stack([Xe,lp]))[:,1]

A_per, A_gene, A_comb = roc_auc_score(yy,po), roc_auc_score(yy,lp), roc_auc_score(yy,pn)
_,_,_,pv = delong(yy, pn, po)     # combined vs perinatal
nri, idi = nri_idi(yy, po, pn)
nris=[];idis=[];das=[]
for _ in range(B):
    idx = RNG.integers(0,len(yy),len(yy))
    if len(np.unique(yy[idx]))<2: continue
    a,b = nri_idi(yy[idx],po[idx],pn[idx]); nris.append(a); idis.append(b)
    das.append(roc_auc_score(yy[idx],pn[idx]) - roc_auc_score(yy[idx],po[idx]))
nlo,nhi=np.percentile(nris,[2.5,97.5]); ilo,ihi=np.percentile(idis,[2.5,97.5]); dlo,dhi=np.percentile(das,[2.5,97.5])
print('GSE32472 day-28 window (all coefficients extrapolated from the day-14 training set, not re-estimated)')
print('  n=%d (cases %d / ctrl %d)' % (len(yy), yy.sum(), (yy==0).sum()))
print('  AUC  perinatal variables only %.4f | four-gene signature only %.4f | combined model %.4f' % (A_per, A_gene, A_comb))
print('  dAUC (combined - perinatal only) = %+.4f (95%% CI %+.3f to %+.3f)   DeLong P = %.4f' % (A_comb-A_per, dlo, dhi, pv))
print('  continuous NRI %+.3f (%.3f, %.3f)   P(NRI<=0)=%.3f' % (nri,nlo,nhi,(np.array(nris)<=0).mean()))
print('  IDI      %+.3f (%.3f, %.3f)   P(IDI<=0)=%.3f' % (idi,ilo,ihi,(np.array(idis)<=0).mean()))
json.dump(dict(n=int(len(yy)),ncase=int(yy.sum()),nctrl=int((yy==0).sum()),
  auc_perinatal=round(float(A_per),4),auc_gene=round(float(A_gene),4),auc_combined=round(float(A_comb),4),
  delta_auc=round(float(A_comb-A_per),4),delta_lo=round(float(dlo),3),delta_hi=round(float(dhi),3),
  delong_p=round(float(pv),4),nri=round(float(nri),3),nri_lo=round(float(nlo),3),nri_hi=round(float(nhi),3),
  nri_le0=round(float((np.array(nris)<=0).mean()),3),idi=round(float(idi),3),idi_lo=round(float(ilo),3),
  idi_hi=round(float(ihi),3),idi_le0=round(float((np.array(idis)<=0).mean()),3)),
  open('_day28_increment.json','w'), ensure_ascii=False, indent=1)
