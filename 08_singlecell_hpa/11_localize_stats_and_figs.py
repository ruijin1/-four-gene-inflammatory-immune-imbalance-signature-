# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 08_singlecell_hpa/11_localize_stats_and_figs.py
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
Four-gene (CD3G / PRTN3 / CD6 / THBS1) localisation: statistics + publication-grade figures (600 dpi)
"""
import os, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

warnings.filterwarnings("ignore")

BASE = P("results")
DATA = os.path.join(BASE, "derived", "scBPD_4gene")
FIG = os.path.join(BASE, "figures", "scBPD_4gene")
os.makedirs(FIG, exist_ok=True)

rcParams["font.family"] = "Arial"
rcParams["font.size"] = 7
rcParams["pdf.fonttype"] = 42

GENES = ["CD3G", "PRTN3", "CD6", "THBS1"]
COEF = {"PRTN3": 0.610, "CD3G": -0.820, "CD6": -1.402, "THBS1": 1.560}
INTERCEPT = -3.375

a = sc.read_h5ad(os.path.join(DATA, "pbmc10k_annotated.h5ad"))
print("loaded:", a.shape)

ORDER = ["CD4 T", "CD8 T", "NK", "B", "Plasma", "CD14 Mono",
         "FCGR3A Mono", "cDC", "pDC", "Platelet"]
present = [c for c in ORDER if c in a.obs["celltype"].unique()]
a.obs["celltype"] = pd.Categorical(a.obs["celltype"],
                                   categories=present, ordered=True)
a = a[~a.obs["celltype"].isna()].copy()
print("\ncell types:", a.obs["celltype"].value_counts().to_dict())

# ============================================================ 1. Localisation statistics
rows = []
for ct in present:
    sub = a[a.obs["celltype"] == ct]
    r = {"cell_type": ct, "n_cells": sub.n_obs}
    for g in GENES:
        x = np.asarray(sub[:, g].X.todense()).ravel()
        # log1p CP10K > 0 is equivalent to at least one UMI detected in that cell
        r[f"{g}_mean"] = round(float(x.mean()), 3)
        r[f"{g}_pct"] = round(100 * float((x > 0).mean()), 1)
    rows.append(r)
tab = pd.DataFrame(rows)
pd.set_option("display.width", 300)
pd.set_option("display.max_columns", 50)
print("\n===== four-gene expression per cell type (mean log1p CP10K / % cells detected / mean UMI) =====")
for g in GENES:
    print(f"\n--- {g}")
    print(tab[["cell_type", "n_cells", f"{g}_mean", f"{g}_pct"]]
          .sort_values(f"{g}_mean", ascending=False).to_string(index=False))

tab.to_csv(os.path.join(DATA, "four_gene_localization_stats.csv"),
           index=False, encoding="utf-8-sig")

# ---- fold enrichment: highest cell type / overall mean
print("\n===== fold enrichment (cell-type mean / overall mean) =====")
for g in GENES:
    overall = float(np.asarray(a[:, g].X.todense()).ravel().mean())
    s = (tab.set_index("cell_type")[f"{g}_mean"] / (overall + 1e-9)).sort_values(ascending=False)
    print(f"\n{g} (overall mean {overall:.3f}):")
    print(s.round(2).to_string())

# ---- signature score
X = np.column_stack([np.asarray(a[:, g].X.todense()).ravel() for g in GENES])
a.obs["sig_score"] = INTERCEPT + X @ np.array([COEF[g] for g in GENES])
print("\n===== four-gene signature score (manuscript coefficients), mean per cell type =====")
print(a.obs.groupby("celltype", observed=True)["sig_score"].mean()
      .sort_values(ascending=False).round(3).to_string())

# ============================================================ 2. Figures
# ---- Fig A: UMAP celltype
fig, ax = plt.subplots(figsize=(4.6, 4.0))
sc.pl.umap(a, color="celltype", ax=ax, show=False, frameon=False,
           title="", size=6, legend_loc="right margin",
           legend_fontsize=6, palette="tab20")
ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
fig.savefig(os.path.join(FIG, "FigA_umap_celltype.tif"), dpi=600,
            bbox_inches="tight")
fig.savefig(os.path.join(FIG, "FigA_umap_celltype.png"), dpi=600,
            bbox_inches="tight")
plt.close(fig)

# ---- Fig B: UMAP features (vmin = 0 throughout, to avoid odd colour scales for all-zero genes)
fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))
for ax, g in zip(axes.ravel(), GENES):
    vmax = max(0.5, float(np.asarray(a[:, g].X.todense()).ravel().max()))
    sc.pl.umap(a, color=g, ax=ax, show=False, frameon=False, title=g,
               size=5, cmap="viridis", vmin=0, vmax=vmax)
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title(g, fontsize=9, fontweight="bold")
fig.suptitle("Cellular localization of the four hub genes (PBMC, 10x 10k v3)",
             fontsize=9)
fig.savefig(os.path.join(FIG, "FigB_umap_four_genes.tif"), dpi=600,
            bbox_inches="tight")
fig.savefig(os.path.join(FIG, "FigB_umap_four_genes.png"), dpi=600,
            bbox_inches="tight")
plt.close(fig)

# ---- Fig C: dotplot
fig, ax = plt.subplots(figsize=(6.4, 2.6))
sc.pl.dotplot(a, GENES, groupby="celltype", ax=ax, show=False,
              cmap="Reds",
              colorbar_title="Mean expression\n(log1p CP10K)",
              size_title="Fraction of\ncells (%)")
ax.set_xlabel(""); ax.set_ylabel("")
fig.savefig(os.path.join(FIG, "FigC_dotplot_four_genes.tif"), dpi=600,
            bbox_inches="tight")
fig.savefig(os.path.join(FIG, "FigC_dotplot_four_genes.png"), dpi=600,
            bbox_inches="tight")
plt.close(fig)

# ---- Fig D: violin
fig, axes = plt.subplots(1, 4, figsize=(11.2, 3.2))
for ax, g in zip(axes.ravel(), GENES):
    sc.pl.violin(a, g, groupby="celltype", ax=ax, show=False,
                 rotation=90, stripplot=False, inner="box", linewidth=0.4)
    ax.set_title(g, fontsize=9, fontweight="bold")
    ax.set_xlabel(""); ax.set_ylabel("Expression (log1p CP10K)")
    ax.tick_params(labelsize=6)
fig.savefig(os.path.join(FIG, "FigD_violin_four_genes.tif"), dpi=600,
            bbox_inches="tight")
fig.savefig(os.path.join(FIG, "FigD_violin_four_genes.png"), dpi=600,
            bbox_inches="tight")
plt.close(fig)

# ---- Fig E: signature score by celltype (bar)
d = a.obs.groupby("celltype", observed=True)["sig_score"].mean().reindex(present)
fig, ax = plt.subplots(figsize=(5.2, 3.0))
colors = ["#c0392b" if v > 0 else "#2980b9" for v in d.values]
ax.bar(range(len(d)), d.values, color=colors)
ax.set_xticks(range(len(d)))
ax.set_xticklabels(d.index, rotation=60, ha="right", fontsize=7)
ax.axhline(0, color="k", lw=0.6)
ax.set_ylabel("Four-gene signature score\n(manuscript coefficients)")
ax.set_title("Signature score by cell type", fontsize=9)
fig.savefig(os.path.join(FIG, "FigE_signature_score.tif"), dpi=600,
            bbox_inches="tight")
fig.savefig(os.path.join(FIG, "FigE_signature_score.png"), dpi=600,
            bbox_inches="tight")
plt.close(fig)

print("\nfigures saved to:", FIG)
for f in sorted(os.listdir(FIG)):
    if f.startswith("Fig"):
        print("  ", f, os.path.getsize(os.path.join(FIG, f)), "bytes")
