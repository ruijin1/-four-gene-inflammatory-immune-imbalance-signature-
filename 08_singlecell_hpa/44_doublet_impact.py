# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 08_singlecell_hpa/44_doublet_impact.py
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
44_doublet_impact.py
Using the Scrublet predictions, re-checks on the correct log1p CP10K matrix whether the

(Note: the comparison block in script 42 mistakenly applied expm1 to the raw counts matrix,
 distorting the values; here the annotated log1p CP10K object is used with barcodes aligned to the Scrublet calls.)
"""
import os, warnings
import numpy as np, pandas as pd, scanpy as sc

warnings.filterwarnings("ignore")
DATA = P("results/scBPD_4gene")
GENES = ["CD3G", "PRTN3", "CD6", "THBS1"]
GM = ["ELANE", "AZU1", "MPO", "CTSG", "LTF", "DEFA4", "MMP8", "S100A9", "CSF3R"]

ann = sc.read_h5ad(os.path.join(DATA, "pbmc10k_annotated.h5ad"))
calls = pd.read_csv(os.path.join(DATA, "scrublet_doublet_calls.csv")).set_index("barcode")
calls = calls.reindex(ann.obs_names)
assert calls["predicted_doublet_r075"].notna().all(), "barcodes not fully aligned"
pred = calls["predicted_doublet_r075"].astype(bool).values
ct = ann.obs["celltype"].astype(str).values
print(f"total cells {ann.n_obs}| predicted doublets {pred.sum()} ({pred.mean()*100:.2f}%)")

X = ann[:, GENES].X
X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
lin = np.expm1(X)
print("maximum on the linear scale (should be < 10, CP10K order):", dict(zip(GENES, lin.max(0).round(3))))

print("\n=== cell type with the highest expression per gene: all vs doublet-removed (mean linear CP10K) ===")
keep = ~pred
for j, g in enumerate(GENES):
    a_ = pd.Series(lin[:, j]).groupby(ct).mean().sort_values(ascending=False)
    b_ = pd.Series(lin[keep, j]).groupby(ct[keep]).mean().sort_values(ascending=False)
    same = "same" if a_.index[0] == b_.index[0] else "**changed**"
    print(f"  {g:6s} all Top1 = {a_.index[0]:12s} ({a_.iloc[0]:6.3f}) | doublet-removed Top1 = "
          f"{b_.index[0]:12s} ({b_.iloc[0]:6.3f})  [{same}]")
    print(f"         all Top3 = {list(a_.index[:3])} | doublet-removed Top3 = {list(b_.index[:3])}")

print("\n=== THBS1 platelet fold enrichment (linear scale) ===")
j = GENES.index("THBS1")
fe_a = lin[:, j][ct == "Platelet"].mean() / lin[:, j].mean()
fe_b = lin[keep, j][ct[keep] == "Platelet"].mean() / lin[keep, j].mean()
print(f"  all {fe_a:.2f} | doublet-removed {fe_b:.2f}  (relative change {abs(fe_b-fe_a)/fe_a*100:.1f}%)")

print("\n=== granulocyte / myeloid marker detection % ===")
gm_ok = [g for g in GM if g in ann.var_names]
Xg = ann[:, gm_ok].X; Xg = Xg.toarray() if hasattr(Xg, "toarray") else np.asarray(Xg)
for k, g in enumerate(gm_ok):
    p1 = (Xg[:, k] > 0).mean() * 100
    p2 = (Xg[keep, k] > 0).mean() * 100
    d = p2 - p1
    print(f"  {g:8s} all {p1:6.2f}% | doublet-removed {p2:6.2f}%  (delta {d:+.2f})")

print(f"\nPRTN3 non-zero cells: all {int((X[:, GENES.index('PRTN3')]>0).sum())} | "
      f"doublet-removed {int((X[keep, GENES.index('PRTN3')]>0).sum())}")

print("\n=== change in cell number per cell type ===")
t = pd.DataFrame({"all": pd.Series(ct).value_counts(),
                  "doublet-removed": pd.Series(ct[keep]).value_counts()}).fillna(0).astype(int)
t["retention %"] = (t["doublet-removed"] / t["all"] * 100).round(1)
print(t.sort_values("all", ascending=False).to_string())
t.to_csv(os.path.join(DATA, "doublet_impact_by_celltype.csv"), encoding="utf-8-sig")
