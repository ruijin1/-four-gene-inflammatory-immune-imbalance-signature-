# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 08_singlecell_hpa/02b_marker_table.py
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
Write a cluster x key-marker mean-expression table for manual annotation review
"""
import os, warnings
import numpy as np, pandas as pd, scanpy as sc

warnings.filterwarnings("ignore")
OUT = P("results/scBPD_4gene")

a = sc.read_h5ad(os.path.join(OUT, "E-MTAB-14507_PBMC_integrated.h5ad"))
full = a.raw.to_adata()
full.obs = a.obs.copy()

for clu in ["leiden_0.4", "leiden_0.6", "leiden_1.0"]:
    if clu in a.obs:
        print(clu, "n_clusters =", a.obs[clu].nunique())

CLU = "leiden_0.6"
full.obs[CLU] = a.obs[CLU].astype(str)

key = ["CD3D", "CD3E", "CD4", "CD8A", "CD8B", "IL7R", "CCR7", "SELL",
       "S100A4", "GZMK", "NKG7", "GNLY", "KLRD1", "GZMB",
       "MS4A1", "CD79A", "CD79B", "MZB1", "XBP1",
       "LYZ", "CD14", "VCAN", "FCN1", "S100A8", "S100A9",
       "FCGR3A", "MS4A7", "HLA-DRA", "CST3", "FCER1A", "CD1C",
       "IL3RA", "CLEC4C", "LILRA4",
       "MPO", "ELANE", "PRTN3", "AZU1", "DEFA4", "FCGR3B", "CSF3R",
       "PPBP", "PF4", "HBB", "HBA1", "ALAS2", "CD34", "SOX4", "SPINK2",
       "THBS1", "CD6", "CD3G"]
key = [g for g in key if g in full.var_names]

df = sc.get.obs_df(full, keys=key)
tab = df.groupby(full.obs[CLU].values).mean()
tab.insert(0, "n", full.obs[CLU].value_counts())

# sort by size and transpose for readability
tab = tab.sort_values("n", ascending=False)
pd.set_option("display.width", 400)
pd.set_option("display.max_columns", 100)
pd.set_option("display.max_rows", 200)
print("\n===== cluster x marker (mean log1p CP10K) =====")
print(tab.round(2).T.to_string())
tab.round(4).to_csv(os.path.join(OUT, f"marker_table_{CLU}.csv"))
print("\nsaved:", f"marker_table_{CLU}.csv")
