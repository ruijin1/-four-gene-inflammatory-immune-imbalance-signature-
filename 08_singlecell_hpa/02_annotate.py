# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 08_singlecell_hpa/02_annotate.py
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
Step 2: cell-type annotation (a marker table is written out first for manual review)
"""
import os, warnings
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

OUT = P("results/scBPD_4gene")
FIG = P("figures/scBPD_4gene")
os.makedirs(FIG, exist_ok=True)

adata = sc.read_h5ad(os.path.join(OUT, "E-MTAB-14507_PBMC_integrated.h5ad"))
print("cells:", adata.n_obs, "genes(kept):", adata.n_vars)
print("raw genes:", adata.raw.shape if adata.raw else None)

# use the full gene set for marker assessment and plotting
full = adata.raw.to_adata()
full.obs = adata.obs.copy()
full.obsm = adata.obsm.copy()
full.obsl = {}
if "X_pca_harmony" in adata.obsm:
    full.obsm["X_pca_harmony"] = adata.obsm["X_pca_harmony"]
full.uns = {}
print("full:", full.shape)

CLU = "leiden_1.0"
full.obs[CLU] = adata.obs[CLU].astype(str)

markers = {
    "T cell":      ["CD3D", "CD3E", "CD3G", "CD2", "IL7R", "TRAC", "LAT"],
    "CD4 T":       ["CD4", "IL7R", "LTB", "CCR7", "SELL", "S100A4", "TCF7"],
    "CD8 T":       ["CD8A", "CD8B", "GZMK", "CCL5", "NKG7"],
    "NK":          ["GNLY", "NKG7", "GZMB", "KLRD1", "KLRF1", "NCAM1", "PRF1"],
    "B cell":      ["MS4A1", "CD79A", "CD79B", "IGHM", "CD19", "TCL1A"],
    "Plasma":      ["MZB1", "JCHAIN", "XBP1", "IGHG1", "IGHA1", "SDC1"],
    "CD14 Mono":   ["LYZ", "CD14", "VCAN", "S100A8", "S100A9", "FCN1", "AIF1"],
    "FCGR3A Mono": ["FCGR3A", "MS4A7", "CDKN1C", "HES4", "LYZ"],
    "cDC":         ["FCER1A", "CST3", "CLEC10A", "CD1C", "SERPINF1", "HLA-DRA"],
    "pDC":         ["IL3RA", "CLEC4C", "LILRA4", "SERPINF1", "GZMB"],
    "Neutrophil":  ["FCGR3B", "CSF3R", "S100A8", "S100A9", "MPO", "ELANE",
                    "PRTN3", "AZU1", "DEFA4", "LYZ"],
    "Platelet":    ["PPBP", "PF4", "TUBB1", "GP9", "ITGA2B"],
    "Erythrocyte": ["HBB", "HBA1", "HBA2", "HBD", "ALAS2"],
    "HSPC":        ["CD34", "GATA2", "SOX4", "SPINK2", "AVP"],
}

# ---- mean marker expression per cluster (log1p scale; unscaled log-normalised values)
rows = []
for ct, genes in markers.items():
    genes = [g for g in genes if g in full.var_names]
    if not genes:
        continue
    df = sc.get.obs_df(full, keys=genes)
    m = df.groupby(full.obs[CLU].values).mean()
    rows.append(m.mean(axis=1).rename(ct))
tab = pd.concat(rows, axis=1)
tab.insert(0, "n_cells", full.obs[CLU].value_counts())
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 50)
print("\n===== cluster x lineage mean expression (log1p CP10K) =====")
print(tab.round(2))

tab.to_csv(os.path.join(OUT, "cluster_marker_means.csv"))

# ---- dot plot of canonical markers (inspect clusters before naming)
sc.pl.dotplot(full, markers, groupby=CLU, use_raw=False, standard_scale="var",
              show=False)
import matplotlib.pyplot as plt
plt.savefig(os.path.join(FIG, "S0_cluster_marker_dotplot.png"), dpi=300,
            bbox_inches="tight")
plt.close()

full.write(os.path.join(OUT, "full_genes.h5ad"))
print("saved full_genes.h5ad")
