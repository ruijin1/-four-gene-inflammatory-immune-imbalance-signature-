# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 08_singlecell_hpa/12_composite_figure.py
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
Composite figure for submission: exactly 17 cm wide, 600 dpi TIFF
A UMAP (cell type, on-data legend) | B four genes x cell type dot plot
C-F four-gene UMAP feature plots
"""
import os, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec

warnings.filterwarnings("ignore")

DATA = P("results/scBPD_4gene")
FIG = P("figures/scBPD_4gene")

rcParams["font.family"] = "Arial"
rcParams["font.size"] = 7
rcParams["axes.labelsize"] = 7
rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42

CM = 1 / 2.54
W, H = 17 * CM, 14 * CM
DPI = 600

GENES = ["CD3G", "PRTN3", "CD6", "THBS1"]

a = sc.read_h5ad(os.path.join(DATA, "pbmc10k_annotated.h5ad"))
order = ["CD4 T", "CD8 T", "NK", "B", "Plasma", "CD14 Mono",
         "FCGR3A Mono", "cDC", "pDC", "Platelet"]
a.obs["celltype"] = pd.Categorical(a.obs["celltype"], categories=order,
                                   ordered=True)

fig = plt.figure(figsize=(W, H))
gs = GridSpec(3, 4, figure=fig,
              left=0.055, right=0.965, top=0.93, bottom=0.075,
              height_ratios=[1.25, 0.78, 0.78], wspace=0.35, hspace=0.55)

# ---- A: UMAP celltype (left half, on-data labels)
axA = fig.add_subplot(gs[0, :2])
sc.pl.umap(a, color="celltype", ax=axA, show=False, frameon=False,
           size=4.5, title="", legend_loc="on data", legend_fontsize=5.5,
           legend_fontoutline=2)
axA.set_xlabel("UMAP1", fontsize=7)
axA.set_ylabel("UMAP2", fontsize=7)
axA.text(-0.10, 1.05, "A", transform=axA.transAxes,
         fontsize=11, fontweight="bold", va="top")

# ---- B: dotplot (right half)
axB = fig.add_subplot(gs[0, 2:])
sc.pl.dotplot(a, GENES, groupby="celltype", ax=axB, show=False,
              cmap="Reds", colorbar_title="Mean expr.\n(log1p CP10K)",
              size_title="Cells\n(%)")
axB.set_xlabel("")
axB.set_ylabel("")
axB.tick_params(axis="y", labelsize=6.5)
axB.text(-0.22, 1.08, "B", transform=axB.transAxes,
         fontsize=11, fontweight="bold", va="top")

# ---- C-F: four gene feature plots
for i, g in enumerate(GENES):
    vmax = max(0.5, float(np.asarray(a[:, g].X.todense()).ravel().max()))
    ax = fig.add_subplot(gs[1 + i // 2, 2 * (i % 2):2 * (i % 2) + 2])
    sc.pl.umap(a, color=g, ax=ax, show=False, frameon=False, title=g,
               size=4, cmap="viridis", vmin=0, vmax=vmax)
    ax.set_xlabel("UMAP1", fontsize=6.5)
    ax.set_ylabel("UMAP2", fontsize=6.5)
    ax.set_title(g, fontsize=8, fontweight="bold")
    ax.text(-0.10, 1.08, chr(ord("C") + i), transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top")

out_tif = os.path.join(FIG, "Fig10_four_gene_localization.tif")
out_png = os.path.join(FIG, "Fig10_four_gene_localization.png")
fig.savefig(out_tif, dpi=DPI)
fig.savefig(out_png, dpi=DPI)
plt.close(fig)

from PIL import Image
im = Image.open(out_tif)
print("saved:", out_tif, os.path.getsize(out_tif), "bytes")
print("size(px):", im.size, "| dpi:", im.info.get("dpi"),
      "| width(cm):", round(im.size[0] / DPI * 2.54, 2),
      "| height(cm):", round(im.size[1] / DPI * 2.54, 2))
