# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 08_singlecell_hpa/42_scrublet_doublet_check.py
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
42_scrublet_doublet_check.py
Supplementary doublet detection (Scrublet) for the single-cell reference atlas, checking

Note: the h5ad does not store counts (layers empty, raw=False), so counts are read back from
    the original 10x filtered_feature_bc_matrix and the 10,667-cell set is reproduced with identical QC thresholds.
"""
import os, warnings
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
BASE = P("results")
DATA = os.path.join(BASE, "derived", "scBPD_4gene")
MTX = os.path.join(DATA, "reference", "filtered_feature_bc_matrix")
GENES = ["CD3G", "PRTN3", "CD6", "THBS1"]
GM = ["ELANE", "AZU1", "MPO", "CTSG", "LTF", "DEFA4", "MMP8", "S100A9", "CSF3R"]

# ---------- 1. Read raw counts + reproduce QC ----------
raw = sc.read_10x_mtx(MTX, var_names="gene_symbols", cache=False)
raw.var_names_make_unique()
print("Raw filtered matrix:", raw.shape)

raw.var["mt"] = raw.var_names.str.startswith(("MT-", "MT."))
sc.pp.calculate_qc_metrics(raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
n0 = raw.n_obs
raw = raw[(raw.obs.n_genes_by_counts >= 300) &
          (raw.obs.n_genes_by_counts <= 6000) &
          (raw.obs.pct_counts_mt <= 15)].copy()
print(f"QC: {n0} -> {raw.n_obs}")

ann = sc.read_h5ad(os.path.join(DATA, "pbmc10k_annotated.h5ad"))
common = raw.obs_names.intersection(ann.obs_names)
print(f"cells shared with the annotated object: {len(common)} (expected 10667)")
raw = raw[common].copy()
ann = ann[common].copy()
raw.obs["celltype"] = ann.obs["celltype"].astype(str).values

# ---------- 2. Scrublet ----------
res = {}
for rate in (0.06, 0.075):
    a = raw.copy()
    sc.pp.scrublet(a, expected_doublet_rate=rate, random_state=0, verbose=False)
    pred = a.obs["predicted_doublet"].astype(bool)
    res[rate] = pred
    print(f"\n=== expected_doublet_rate={rate} ===")
    print(f"  predicted doublets: {int(pred.sum())} / {len(pred)} = {pred.mean()*100:.2f}%")
    print("  Distribution by annotated cell type:")
    print(pd.crosstab(a.obs["celltype"], pred, normalize="index").round(3)
          .rename(columns={False: "singlet", True: "doublet"}).to_string())

pred = res[0.075]
print("\n" + "=" * 70)
print("=== Do the key conclusions change after removing predicted doublets (rate = 0.075)? ===")
b = raw[~pred.values].copy()
print(f"remaining cells: {b.n_obs}")

Xa = ann[:, GENES].X; Xa = Xa.toarray() if hasattr(Xa, "toarray") else np.asarray(Xa)
Xb_full = b[:, GENES].X; Xb_full = Xb_full.toarray() if hasattr(Xb_full, "toarray") else np.asarray(Xb_full)
la, lb = np.expm1(Xa), np.expm1(Xb_full)
ca = ann.obs["celltype"].astype(str).values
cb = b.obs["celltype"].astype(str).values

print("\n--- cell type with the highest expression per gene (all vs doublet-removed) ---")
for j, g in enumerate(GENES):
    ta = pd.Series(la[:, j]).groupby(ca).mean().sort_values(ascending=False)
    tb = pd.Series(lb[:, j]).groupby(cb).mean().sort_values(ascending=False)
    print(f"  {g}: all -> {ta.index[0]} ({ta.iloc[0]:.3f}) | doublet-removed -> {tb.index[0]} ({tb.iloc[0]:.3f})")

print("\n--- THBS1 platelet fold enrichment (linear scale) ---")
for tag, l, c in (("all", la, ca), ("doublet-removed", lb, cb)):
    m = pd.Series(l[:, GENES.index("THBS1")]).groupby(c).mean()
    fe = m["Platelet"] / l[:, GENES.index("THBS1")].mean()
    print(f"  {tag}: {fe:.2f}")

print("\n--- granulocyte marker detection % (all vs doublet-removed) ---")
gm_ok = [g for g in GM if g in b.var_names]
Xg_a = ann[:, gm_ok].X; Xg_a = Xg_a.toarray() if hasattr(Xg_a, "toarray") else np.asarray(Xg_a)
Xg_b = b[:, gm_ok].X;   Xg_b = Xg_b.toarray() if hasattr(Xg_b, "toarray") else np.asarray(Xg_b)
for j, g in enumerate(gm_ok):
    print(f"  {g:8s} all {(Xg_a[:,j]>0).mean()*100:6.2f}% | doublet-removed {(Xg_b[:,j]>0).mean()*100:6.2f}%")

print(f"\nPRTN3 non-zero cells (doublet-removed): {int((Xb_full[:, GENES.index('PRTN3')] > 0).sum())}")
out = pd.DataFrame({"barcode": raw.obs_names,
                    "celltype": raw.obs["celltype"].values,
                    "doublet_score": None,
                    "predicted_doublet_r06": res[0.06].values,
                    "predicted_doublet_r075": pred.values})
out.to_csv(os.path.join(DATA, "scrublet_doublet_calls.csv"), index=False, encoding="utf-8-sig")
print("\nSaved:", os.path.join(DATA, "scrublet_doublet_calls.csv"))
