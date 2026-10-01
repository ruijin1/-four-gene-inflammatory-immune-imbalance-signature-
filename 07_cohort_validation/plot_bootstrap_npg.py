# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/plot_bootstrap_npg.py
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

"""Fig. S1. Bootstrap internal-validation figure -- colours match Figure 5 (ggsci npg palette)"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve

XLSX = P("data/raw/GSE32472/GSE32472_series_matrix.txt/正常组与中重度BPD的比较/12-17d/炎症并集/Nomogram/预测值-系数.xlsx")
CSV = P("results/bootstrap_auc_iterations_1000.csv")

# ggsci npg palette (matches Figure 5)
RED = '#E64B35'        # npg 1 -- Figure 5 main red
BLUE = '#4DBBD5'       # npg 2 -- Figure 5 main blue
GREEN = '#00A087'      # npg 3
DARKBLUE = '#3C5488'   # npg 4
ORANGE = '#F39B7F'     # npg 5
GREYBLUE = '#8491B4'   # npg 6
LIGHTBLUE = '#A5DBE8'  # light blue used in Figure 5
LIGHTRED = '#F1A399'   # light red used in Figure 5

df = pd.read_excel(XLSX, sheet_name='glm模型')
y = df['outcome'].astype(int).values
lp = df['linear_predictors'].values
auc_app = roc_auc_score(y, lp)

b = pd.read_csv(CSV)
corr = b['corrected_auc'].values
opt = b['optimism'].values
auc_corr = auc_app - opt.mean()
ci = np.percentile(corr, [2.5, 97.5])
lo, hi = 0.8657, 0.9911  # DeLong apparent CI

plt.rcParams.update({'font.family': 'Arial', 'font.size': 8,
                     'axes.linewidth': 0.8, 'xtick.major.width': 0.8,
                     'ytick.major.width': 0.8})
fig, axes = plt.subplots(1, 2, figsize=(6.69, 2.95), dpi=600)
fig.subplots_adjust(left=0.09, right=0.985, top=0.86, bottom=0.17, wspace=0.32)

# --- Panel A: ROC ---
ax = axes[0]
fpr, tpr, _ = roc_curve(y, lp)
ax.plot(fpr, tpr, color=RED, lw=1.6,
        label=f'Training set (Day 14)\nAUC = {auc_app:.3f} (95% CI {lo:.3f}-{hi:.3f})')
ax.plot([0, 1], [0, 1], color=GREYBLUE, lw=0.9, ls='--', label='Chance')
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
ax.set_xlabel('1 - Specificity (False positive rate)', fontsize=8)
ax.set_ylabel('Sensitivity (True positive rate)', fontsize=8)
ax.legend(loc='lower right', frameon=False, fontsize=6.6, handlelength=1.6)
ax.set_title('A   Apparent ROC curve', fontsize=8.5, loc='left', pad=6, fontweight='bold')
ax.tick_params(labelsize=7.2)

# --- Panel B: bootstrap distribution ---
ax2 = axes[1]
ax2.hist(corr, bins=40, color=LIGHTBLUE, edgecolor=BLUE, linewidth=0.4)
ax2.axvline(auc_app, color=RED, lw=1.4, ls='-', label=f'Apparent AUC = {auc_app:.3f}')
ax2.axvline(auc_corr, color=DARKBLUE, lw=1.4, ls='--',
            label=f'Optimism-corrected AUC = {auc_corr:.3f}')
ax2.axvspan(ci[0], ci[1], color=DARKBLUE, alpha=0.10,
            label=f'95% CI {ci[0]:.3f}-{ci[1]:.3f}')
ax2.axvline(0.5, color=GREYBLUE, lw=0.9, ls=':')
ax2.set_xlim(0.60, 1.00)
ax2.set_xlabel('Optimism-corrected AUC (1000 bootstrap replicates)', fontsize=8)
ax2.set_ylabel('Frequency', fontsize=8)
ax2.legend(loc='upper left', frameon=False, fontsize=6.6, handlelength=1.8)
ax2.set_title('B   Bootstrap internal validation', fontsize=8.5, loc='left', pad=6,
              fontweight='bold')
ax2.tick_params(labelsize=7.2)

for a in axes:
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)

out = P("figures/S1_Fig_Bootstrap_internal_validation_AUC_600dpi.tiff")
fig.savefig(out, dpi=600, format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
fig.savefig(out.replace('.tiff', '.png'), dpi=600)
print('saved:', out)
print(f'apparent={auc_app:.4f}  corrected={auc_corr:.4f}  CI={ci[0]:.4f}-{ci[1]:.4f}  '
      f'optimism={opt.mean():.4f}')
