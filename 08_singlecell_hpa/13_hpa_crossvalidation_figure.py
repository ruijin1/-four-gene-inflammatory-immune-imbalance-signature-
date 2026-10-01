# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 08_singlecell_hpa/13_hpa_crossvalidation_figure.py
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
Supplementary Figure -- cross-validation against HPA immune-cell bulk RNA-seq
A: heat map of 14 genes x 19 immune cell types (log10(nTPM+1), row-normalised to 0-1)
B: nTPM of the four genes across key immune cell types (with the FCGR3B positive control)
exactly 17 cm wide x 600 dpi
"""
import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap

warnings.filterwarnings("ignore")

FIG = P("figures/scBPD_4gene")
DATA = P("results/scBPD_4gene")

rcParams["font.family"] = "Arial"
rcParams["font.size"] = 7
rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42

CM = 1 / 2.54
DPI = 600

csv = os.path.join(DATA, "hpa_immune_cell_crossvalidation.csv")
df = pd.read_csv(csv, index_col=0)

ORDER = ['neutrophil', 'basophil', 'eosinophil',
         'classical monocyte', 'intermediate monocyte', 'non-classical monocyte',
         'myeloid DC', 'plasmacytoid DC',
         'naive CD4 T-cell', 'memory CD4 T-cell', 'naive CD8 T-cell',
         'memory CD8 T-cell', 'T-reg', 'MAIT T-cell', 'gdT-cell',
         'NK-cell', 'naive B-cell', 'memory B-cell', 'total PBMC']

GENE_ORDER = ['FCGR3B', 'S100A9', 'CSF3R', 'ELANE', 'DEFA4', 'LTF', 'MMP8',
              'CTSG', 'MPO', 'AZU1', 'PRTN3', 'THBS1', 'CD6', 'CD3G']

SHORT = {'neutrophil': 'Neutrophil', 'basophil': 'Basophil',
         'eosinophil': 'Eosinophil', 'classical monocyte': 'cMono',
         'intermediate monocyte': 'iMono', 'non-classical monocyte': 'ncMono',
         'myeloid DC': 'mDC', 'plasmacytoid DC': 'pDC',
         'naive CD4 T-cell': 'nCD4', 'memory CD4 T-cell': 'mCD4',
         'naive CD8 T-cell': 'nCD8', 'memory CD8 T-cell': 'mCD8',
         'T-reg': 'Treg', 'MAIT T-cell': 'MAIT', 'gdT-cell': 'gdT',
         'NK-cell': 'NK', 'naive B-cell': 'nB', 'memory B-cell': 'mB',
         'total PBMC': 'total\nPBMC'}

df = df[ORDER]
df = df.loc[[g for g in GENE_ORDER if g in df.index]]

# log transform + row-normalisation to 0-1
L = np.log10(df.values + 1.0)
rowmin = L.min(axis=1, keepdims=True)
rowmax = L.max(axis=1, keepdims=True)
N = np.divide(L - rowmin, np.where(rowmax - rowmin == 0, 1, rowmax - rowmin))

W, H = 17 * CM, 11.5 * CM
fig = plt.figure(figsize=(W, H))
gs = GridSpec(2, 1, figure=fig, left=0.10, right=0.985, top=0.90, bottom=0.20,
              height_ratios=[1.35, 1.0], hspace=0.62)

# ---------- A: heatmap ----------
axA = fig.add_subplot(gs[0])
cmap = LinearSegmentedColormap.from_list("w2r", ["#FFFFFF", "#FDE0D0", "#E8642B", "#8B1A0A"])
im = axA.imshow(N, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
axA.set_xticks(range(len(ORDER)))
axA.set_xticklabels([SHORT[c] for c in ORDER], rotation=60, ha="right", fontsize=5.6)
axA.set_yticks(range(len(df.index)))
axA.set_yticklabels(df.index, fontsize=6.6)
for g in ['PRTN3', 'CD6', 'CD3G', 'THBS1']:
    if g in df.index:
        axA.get_yticklabels()[list(df.index).index(g)].set_fontweight("bold")

# underscore separating cell-type groups
groups = [("Granulo", 0, 3), ("Mono", 3, 6), ("DC", 6, 8),
          ("T", 8, 15), ("NK/B", 15, 18), ("PBMC", 18, 19)]
for name, i0, i1 in groups:
    x = (i0 + i1 - 1) / 2
    axA.text(x, len(df.index) + 0.35, name, ha="center", va="bottom",
             fontsize=6.2, fontweight="bold")
    axA.plot([i0 - 0.5, i1 - 0.5], [len(df.index) - 0.45, len(df.index) - 0.45],
             color="black", lw=0.8, clip_on=False)

axA.set_title("HPA immune-cell bulk RNA-seq (nTPM, row-scaled $\\log_{10}$)",
              fontsize=7.5, pad=4)
axA.text(-0.115, 1.03, "A", transform=axA.transAxes, fontsize=11,
         fontweight="bold", va="top")

cb = fig.colorbar(im, ax=axA, fraction=0.012, pad=0.008)
cb.ax.tick_params(labelsize=5.5)
cb.set_label("row-scaled", fontsize=5.8)

# ---------- B: four-gene comparison ----------
axB = fig.add_subplot(gs[1])
show = ['neutrophil', 'eosinophil', 'classical monocyte',
        'memory CD4 T-cell', 'NK-cell', 'memory B-cell', 'total PBMC']
labels = [SHORT[c].replace("\n", " ") for c in show]
bg = ['CD3G', 'PRTN3', 'CD6', 'THBS1']
colors = {"CD3G": "#2C6FA8", "PRTN3": "#C0392B", "CD6": "#17A589", "THBS1": "#D68910"}

x = np.arange(len(show))
w = 0.20
for i, g in enumerate(bg):
    vals = df.loc[g, show].values.astype(float)
    axB.bar(x + (i - 1.5) * w, vals, w, label=g, color=colors[g],
            edgecolor="black", linewidth=0.35)

axB.set_yscale("symlog", linthresh=1.0, linscale=0.6)
axB.set_ylabel("nTPM (symlog)", fontsize=7)
axB.set_xticks(x)
axB.set_xticklabels(labels, rotation=30, ha="right", fontsize=6.2)
axB.legend(fontsize=6.2, frameon=False, ncol=4, loc="upper center",
           bbox_to_anchor=(0.5, 1.20))
axB.axhline(0, color="black", lw=0.5)
axB.tick_params(axis="y", labelsize=6)
axB.text(-0.115, 1.10, "B", transform=axB.transAxes, fontsize=11,
         fontweight="bold", va="top")
axB.grid(axis="y", ls=":", lw=0.3, alpha=0.5)
axB.set_axisbelow(True)

out_tif = os.path.join(FIG, "FigS2_HPA_crossvalidation.tif")
out_png = os.path.join(FIG, "FigS2_HPA_crossvalidation.png")
fig.savefig(out_tif, dpi=DPI)
fig.savefig(out_png, dpi=DPI)
plt.close(fig)

from PIL import Image
im2 = Image.open(out_tif)
print("saved:", out_tif, os.path.getsize(out_tif), "bytes")
print("size(px):", im2.size, "| dpi:", im2.info.get("dpi"),
      "| width(cm):", round(im2.size[0] / DPI * 2.54, 2),
      "| height(cm):", round(im2.size[1] / DPI * 2.54, 2))
