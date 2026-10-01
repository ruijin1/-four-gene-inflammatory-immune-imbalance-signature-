# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 08_singlecell_hpa/10_ref_pbmc_localize.py
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
Approach A: localise the four manuscript genes using the 10x official PBMC 10k v3
reference dataset (CD3G / PRTN3 / CD6 / THBS1)

Pipeline: load -> QC -> normalise -> HVG -> PCA -> neighbourhood graph -> Leiden -> UMAP
      -> marker annotation -> four-gene localisation -> figures (600 dpi)

Run:
  python 10_ref_pbmc_localize.py
"""
import os, tarfile, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1
sc.settings.n_jobs = 4

BASE = P("results")
DATA = os.path.join(BASE, "derived", "scBPD_4gene")
REF = os.path.join(DATA, "reference")
FIG = os.path.join(BASE, "figures", "scBPD_4gene")
os.makedirs(FIG, exist_ok=True)

GENES = ["CD3G", "PRTN3", "CD6", "THBS1"]

# ---------------------------------------------------------------- load
tgz = os.path.join(REF, "pbmc_10k_v3.tar.gz")
mdir = os.path.join(REF, "filtered_feature_bc_matrix")
if not os.path.exists(mdir):
    with tarfile.open(tgz) as t:
        t.extractall(REF)
    print("extracted ->", mdir)
print("matrix dir:", mdir, os.listdir(mdir))

adata = sc.read_10x_mtx(mdir, var_names="gene_symbols", cache=True)
adata.var_names_make_unique()
adata.obs_names_make_unique()
print("raw:", adata.shape)

# ---------------------------------------------------------------- QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"],
                           percent_top=None, log1p=False, inplace=True)
print(adata.obs[["n_genes_by_counts", "total_counts",
                 "pct_counts_mt"]].describe().round(2))

n0 = adata.n_obs
adata = adata[(adata.obs.n_genes_by_counts >= 300) &
              (adata.obs.n_genes_by_counts <= 6000) &
              (adata.obs.pct_counts_mt <= 15)].copy()
print(f"QC: {n0} -> {adata.n_obs} cells")

for g in GENES:
    print(f"  {g} in var_names: {g in adata.var_names}")

# ---------------------------------------------------------------- norm / dim red
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1.0, key_added="leiden", flavor="igraph",
             n_iterations=2, directed=False)
print("clusters:", adata.obs["leiden"].nunique())

# ---------------------------------------------------------------- marker scores
full = adata.raw.to_adata()
full.obs = adata.obs.copy()
full.obsm = adata.obsm.copy()
full.uns = {}

markers = {
    "CD4 T":        ["IL7R", "CCR7", "LTB", "TCF7", "LEF1", "SELL", "CD4",
                     "LDHB", "NOSIP", "MAL"],
    "CD8 T":        ["CD8A", "CD8B", "GZMK", "CCL5", "GZMA", "KLRG1", "DUSP2"],
    "NK":           ["GNLY", "NKG7", "GZMB", "KLRD1", "KLRF1", "NCAM1",
                     "PRF1", "SPON2", "FGFBP2"],
    "B":            ["MS4A1", "CD79A", "CD79B", "IGHM", "TCL1A", "CD19",
                     "BANK1", "VPREB3"],
    "Plasma":       ["MZB1", "JCHAIN", "XBP1", "IGHG1", "IGHA1", "SSR4",
                     "DERL3"],
    "CD14 Mono":    ["LYZ", "CD14", "VCAN", "FCN1", "S100A8", "S100A9",
                     "AIF1", "FTL", "TYROBP"],
    "FCGR3A Mono":  ["FCGR3A", "MS4A7", "CDKN1C", "HES4", "LST1", "IFITM3",
                     "RHOC", "CFD"],
    "cDC":          ["FCER1A", "CST3", "CLEC10A", "CD1C", "SERPINF1", "ENHO",
                     "CLEC9A", "C1orf54"],
    "pDC":          ["IL3RA", "CLEC4C", "LILRA4", "IRF7", "IRF8", "TPM2",
                     "SERPINF1", "PLD4"],
    "Platelet":     ["PPBP", "PF4", "TUBB1", "GP9", "ITGA2B", "NRGN",
                     "CLU", "CAVIN2"],
    "Neutrophil":   ["FCGR3B", "CSF3R", "MPO", "ELANE", "PRTN3", "AZU1",
                     "DEFA4", "CTSG", "S100A12", "G0S2"],
}

avail = {k: [g for g in v if g in full.var_names] for k, v in markers.items()}
for k, v in avail.items():
    sc.tl.score_genes(full, v, score_name=f"sc_{k}")

sc_cols = [f"sc_{k}" for k in markers]
scores = full.obs.groupby("leiden")[sc_cols].mean()
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
print("\n===== cluster x lineage score =====")
print(scores.round(2))
print("\nassigned (argmax):")
print(scores.idxmax(axis=1).str.replace("sc_", "", regex=False))

lab = scores.idxmax(axis=1).str.replace("sc_", "", regex=False)
full.obs["celltype"] = full.obs["leiden"].map(lab).astype(str)
adata.obs["celltype"] = full.obs["celltype"].values
print("\ncell type counts:")
print(full.obs["celltype"].value_counts())

full.write(os.path.join(DATA, "pbmc10k_annotated.h5ad"))
print("saved pbmc10k_annotated.h5ad")
