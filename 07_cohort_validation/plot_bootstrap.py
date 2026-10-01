# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/plot_bootstrap.py
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

import numpy as np, pandas as pd, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve

XL = P("data/raw/GSE32472/GSE32472_series_matrix.txt/正常组与中重度BPD的比较/12-17d/炎症并集/Nomogram/预测值-系数.xlsx")
GENES = ['CD3G', 'PRTN3', 'CD6', 'THBS1']
df = pd.read_excel(XL, sheet_name='glm模型')
y = df['outcome'].astype(int).values; X = df[GENES].values; lp = df['linear_predictors'].values
n = len(y); auc_app = roc_auc_score(y, lp)

rng = np.random.default_rng(20260909)
opt = []
for _ in range(1000):
    idx = rng.integers(0, n, n)
    Xb, yb = X[idx], y[idx]
    if len(np.unique(yb)) < 2: continue
    m = LogisticRegression(C=np.inf, max_iter=5000, solver='lbfgs').fit(Xb, yb)
    opt.append(roc_auc_score(yb, m.decision_function(Xb)) - roc_auc_score(y, m.decision_function(X)))
opt = np.array(opt); auc_corr = auc_app - opt.mean()
ci = np.percentile(auc_app - opt, [2.5, 97.5])
lo, hi = 0.8657, 0.9911  # DeLong

plt.rcParams.update({'font.family': 'Arial', 'font.size': 8,
                     'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8})
fig, axes = plt.subplots(1, 2, figsize=(6.69, 2.95), dpi=600)
fig.subplots_adjust(left=0.09, right=0.985, top=0.86, bottom=0.17, wspace=0.32)

# --- Panel A: ROC ---
ax = axes[0]
fpr, tpr, _ = roc_curve(y, lp)
ax.plot(fpr, tpr, color='#B03030', lw=1.6,
        label=f'Training set (Day 14)\nAUC = {auc_app:.3f} (95% CI {lo:.3f}-{hi:.3f})')
ax.plot([0, 1], [0, 1], color='#9AA0AA', lw=0.9, ls='--', label='Chance')
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
ax.set_xlabel('1 - Specificity (False positive rate)', fontsize=8)
ax.set_ylabel('Sensitivity (True positive rate)', fontsize=8)
ax.legend(loc='lower right', frameon=False, fontsize=6.6, handlelength=1.6)
ax.set_title('A   Apparent ROC curve', fontsize=8.5, loc='left', pad=6, fontweight='bold')
ax.tick_params(labelsize=7.2)

# --- Panel B: bootstrap distribution ---
ax2 = axes[1]
corr = auc_app - opt
ax2.hist(corr, bins=40, color='#D9E2F3', edgecolor='#5B7FB5', linewidth=0.4)
ax2.axvline(auc_app, color='#B03030', lw=1.4, ls='-', label=f'Apparent AUC = {auc_app:.3f}')
ax2.axvline(auc_corr, color='#1F5FA8', lw=1.4, ls='--',
            label=f'Optimism-corrected AUC = {auc_corr:.3f}')
ax2.axvspan(ci[0], ci[1], color='#1F5FA8', alpha=0.10,
            label=f'95% CI {ci[0]:.3f}-{ci[1]:.3f}')
ax2.axvline(0.5, color='#9AA0AA', lw=0.9, ls=':')
ax2.set_xlim(0.60, 1.00)
ax2.set_xlabel('Optimism-corrected AUC (1000 bootstrap replicates)', fontsize=8)
ax2.set_ylabel('Frequency', fontsize=8)
ax2.legend(loc='upper left', frameon=False, fontsize=6.6, handlelength=1.8)
ax2.set_title('B   Bootstrap internal validation', fontsize=8.5, loc='left', pad=6, fontweight='bold')
ax2.tick_params(labelsize=7.2)

for a in axes:
    a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)

out = P("figures/FigS_Bootstrap_internal_validation_AUC_600dpi.tiff")
fig.savefig(out, dpi=600, format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
fig.savefig(out.replace('.tiff', '.png'), dpi=600)
print('saved:', out)
print(f'apparent={auc_app:.4f}  corrected={auc_corr:.4f}  CI={ci[0]:.4f}-{ci[1]:.4f}  optimism={opt.mean():.4f}')
