# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 08_singlecell_hpa/01_load_qc_integrate.py
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
E-MTAB-14507  (infant PBMC, 5 BPD + 5 HC, 10x 3' v1)
Step 1: load the 10 samples -> QC -> integrate -> cluster -> save h5ad

Run:
#   python <script>.py
"""
import os, glob, gzip, warnings
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import scanpy as sc
import anndata as ad

warnings.filterwarnings("ignore")
sc.settings.verbosity = 2
sc.settings.n_jobs = 4

DATA = P("data/raw/E-MTAB-14507")
OUT = P("results/scBPD_4gene")
os.makedirs(OUT, exist_ok=True)

SAMPLES = ["BPD1", "BPD2", "BPD3", "BPD4", "BPD5",
           "HC1", "HC2", "HC3", "HC4", "HC5"]


def read_10x(sample):
    """Read matrix.mtx.gz / barcodes / features (features has two columns only: ENSG, symbol)"""
    mtx = os.path.join(DATA, f"matrix-{sample}.mtx.gz")
    bc = os.path.join(DATA, f"barcodes-{sample}.tsv.gz")
    ft = os.path.join(DATA, f"features-{sample}.tsv.gz")

    with gzip.open(mtx, "rt") as fh:
        M = sio.mmread(fh).tocsc()          # genes x cells
    M = M.T.tocsr()                          # cells x genes

    barcodes = pd.read_csv(bc, header=None)[0].astype(str).values
    feats = pd.read_csv(ft, header=None, sep="\t")
    feats.columns = ["gene_ids", "gene_symbols"][: feats.shape[1]]

    a = ad.AnnData(M.astype(np.float32))
    a.obs_names = [f"{sample}_{b.split('-')[0]}" for b in barcodes]
    a.var_names = feats["gene_symbols"].astype(str).values
    a.var["gene_ids"] = feats["gene_ids"].astype(str).values
    a.obs["sample"] = sample
    a.obs["group"] = "BPD" if sample.startswith("BPD") else "HC"
    a.var_names_make_unique()
    return a


# ---------------------------------------------------------------- load
adatas = []
for s in SAMPLES:
    a = read_10x(s)
    print(f"{s}: {a.n_obs} cells x {a.n_vars} genes")
    adatas.append(a)

adata = ad.concat(adatas, join="outer", label="batch", keys=SAMPLES,
                  index_unique=None)
adata.obs_names_make_unique()
print("merged:", adata.shape)

# ---------------------------------------------------------------- QC
adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.match(r"^HB[ABDGQZE]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], percent_top=None,
    log1p=False, inplace=True)

print(adata.obs[["n_genes_by_counts", "total_counts",
                 "pct_counts_mt"]].describe())

adata.layers["counts"] = adata.X.copy()

# --- Filtering thresholds (mitochondrial threshold relaxed for neonatal PBMC) ---
MIN_GENES = 300
MAX_GENES = 6000
MIN_COUNTS = 500
MAX_MT = 20

qc_before = adata.n_obs
adata = adata[(adata.obs.n_genes_by_counts >= MIN_GENES) &
              (adata.obs.n_genes_by_counts <= MAX_GENES) &
              (adata.obs.total_counts >= MIN_COUNTS) &
              (adata.obs.pct_counts_mt <= MAX_MT)].copy()
print(f"QC: {qc_before} -> {adata.n_obs} cells")

print(adata.obs.groupby("sample").size())

# ---------------------------------------------------------------- norm / HVG
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat",
                            batch_key="sample")
print("HVG:", int(adata.var.highly_variable.sum()))

adata = adata[:, adata.var.highly_variable].copy()

sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)

# ---------------------------------------------------------------- PCA + harmony
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
sc.external.pp.harmony_integrate(adata, key="sample", basis="X_pca",
                                 adjusted_basis="X_pca_harmony",
                                 max_iter_harmony=20)

sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca_harmony")
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1.0, key_added="leiden_1.0", flavor="igraph",
             n_iterations=2, directed=False)
for r in [0.4, 0.6, 1.2]:
    sc.tl.leiden(adata, resolution=r, key_added=f"leiden_{r}",
                 flavor="igraph", n_iterations=2, directed=False)

adata.obs["sample"] = adata.obs["sample"].astype("category")
adata.obs["group"] = adata.obs["group"].astype("category")

h5 = os.path.join(OUT, "E-MTAB-14507_PBMC_integrated.h5ad")
adata.write_h5ad(h5)
print("saved:", h5)
print(adata.obs["leiden_1.0"].value_counts())
