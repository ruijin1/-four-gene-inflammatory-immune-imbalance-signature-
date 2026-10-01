# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/external_validation_gse32472.py
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

"""GSE220135 external validation: using the four-gene logistic model with coefficients fixed
- fixed coefficients: logit(P) = -3.375 + 0.610*PRTN3 - 0.820*CD3G - 1.402*CD6 + 1.560*THBS1
- validation features: DESeq2 vst (same log2 scale as the training-set expression values)
- output: external-validation ROC figures (PNG + PDF, with the training set as reference) and a per-sample prediction Excel
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import roc_curve, roc_auc_score
from pydeseq2.dds import DeseqDataSet

GENES = ["CD3G", "PRTN3", "CD6", "THBS1"]   # independent of the row order in the coefficient file; used explicitly below
COEF = {"(Intercept)": -3.375128, "CD3G": -0.819751,
        "PRTN3": 0.610007, "CD6": -1.402259, "THBS1": 1.560158}
TRAIN_XL = P("data/raw/GSE32472/GSE32472_series_matrix.txt/正常组与中重度BPD的比较/12-17d/炎症并集/logestic回归/预测值-系数.xlsx")
CSVs = {
    "14d": P("data/raw/GSE220135/GSE220135_14d.csv"),
    "28d": P("data/raw/GSE220135/GSE220135_28d.csv"),
}
OUT_IMG = P("figures/GSE220135")
OUT_XL = P("results/GSE220135/external_validation_GSE32472_model_GSE220135.xlsx")
os.makedirs(OUT_IMG, exist_ok=True)


def get_vst(csv_path):
    counts = pd.read_csv(csv_path, index_col=0).astype(int).T
    grp = {c: ("nonBPD" if "NONBPD" in c.upper() else "BPD") for c in counts.index}
    meta = pd.DataFrame({"condition": [grp[c] for c in counts.index]},
                        index=counts.index)
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition")
    dds.deseq2(); dds.vst_fit()
    v = pd.DataFrame(dds.vst_transform(), index=counts.index, columns=counts.columns)
    return v, grp


def predict(vst_row):
    lp = COEF["(Intercept)"]
    for g in GENES:
        lp += COEF[g] * vst_row[g]
    return 1.0 / (1.0 + np.exp(-lp))


def delong_ci(y, proba):
    """AUC 95% CI via the DeLong variance."""
    y = np.asarray(y); p = np.asarray(proba)
    pos = y == 1; neg = y == 0
    m = int(pos.sum()); n = int(neg.sum())
    xp, xn = p[pos], p[neg]
    v10 = np.array([(np.sum(xp[i] > xn) + 0.5 * np.sum(xp[i] == xn)) / n for i in range(m)])
    v01 = np.array([(np.sum(xp > xn[j]) + 0.5 * np.sum(xp == xn[j])) / m for j in range(n)])
    auc = v10.mean()
    var = np.var(v10, ddof=1) / m + np.var(v01, ddof=1) / n
    return auc, auc - 1.96 * np.sqrt(var), auc + 1.96 * np.sqrt(var)


def youden(y, proba):
    fpr, tpr, thr = roc_curve(y, proba)
    j = tpr - fpr
    i = np.argmax(j)
    tn = ((proba <= thr[i]) & (y == 0)).sum(); fn = ((proba <= thr[i]) & (y == 1)).sum()
    fp = ((proba > thr[i]) & (y == 0)).sum(); tp = ((proba > thr[i]) & (y == 1)).sum()
    return thr[i], tp / (tp + fn), tn / (tn + fp)


def plot_roc(ax, y, proba, auc, title, color="#185FA5", ref=None):
    fpr, tpr, _ = roc_curve(y, proba)
    ax.plot(fpr, tpr, color=color, lw=2.4, label=f"GSE220135 {title} (AUC={auc:.3f})")
    if ref is not None:
        ax.plot(ref[0], ref[1], color="#888888", lw=1.6, ls="--",
                label=f"Training set GSE32472 (AUC={ref[2]:.3f})")
    ax.plot([0, 1], [0, 1], ls=":", color="#BBBBBB", lw=1)
    ax.set_xlabel("False Positive Rate", fontsize=11)
    ax.set_ylabel("True Positive Rate", fontsize=11)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=10)


def main():
    # training-set ROC (reference)
    tr = pd.read_excel(TRAIN_XL, sheet_name="glm模型")
    tr_fpr, tr_tpr, _ = roc_curve(tr.outcome, tr.predicted_probability)
    tr_auc = roc_auc_score(tr.outcome, tr.predicted_probability)
    ref = (tr_fpr, tr_tpr, tr_auc)

    rows = []
    detail = []
    figs = {}
    for tag in ["14d", "28d"]:
        v, grp = get_vst(CSVs[tag])
        y = np.array([1 if grp[s] == "BPD" else 0 for s in v.index])
        proba = np.array([predict(v.loc[s]) for s in v.index])
        auc, lo, hi = delong_ci(y, proba)
        thr, sens, spec = youden(y, proba)
        print(f"{tag}: n={len(y)} (BPD={y.sum()}), external validation AUC={auc:.3f} "
              f"(95%CI {lo:.3f}-{hi:.3f}), youden_cutoff={thr:.3f}, Se={sens:.3f}, Sp={spec:.3f}")
        rows.append(dict(timepoint=tag, n_validation_samples=len(y), BPD=y.sum(), control=len(y)-y.sum(),
                         external_validation_auc=round(auc, 4), AUC_95CI=f"{lo:.3f}-{hi:.3f}",
                         youden_cutoff=round(thr, 4), sensitivity=round(sens, 4), specificity=round(spec, 4),
                         training_auc_reference=round(tr_auc, 4)))
        # per sample
        for s in v.index:
            detail.append(dict(timepoint=tag, sample=s, true_group=grp[s],
                               predicted_probability=round(proba[v.index.tolist().index(s)], 4)))
        fig, ax = plt.subplots(figsize=(6.4, 5.6))
        fig.patch.set_facecolor("white")
        plot_roc(ax, y, proba, auc, tag, ref=ref)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_IMG, f"external_validation_ROC_{tag}.png"), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(OUT_IMG, f"external_validation_ROC_{tag}.pdf"), bbox_inches="tight")
        plt.close(fig)
        figs[tag] = (y, proba, auc)

    # combined figure
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6))
    fig.patch.set_facecolor("white")
    for ax, tag in zip(axes, ["14d", "28d"]):
        y, proba, auc = figs[tag]
        plot_roc(ax, y, proba, auc, tag, ref=ref)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_IMG, "external_validation_ROC_14d_28d.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUT_IMG, "external_validation_ROC_14d_28d.pdf"), bbox_inches="tight")
    plt.close(fig)

    with pd.ExcelWriter(OUT_XL) as w:
        pd.DataFrame(rows).to_excel(w, sheet_name="Validation performance summary", index=False)
        pd.DataFrame(detail).to_excel(w, sheet_name="Per-sample predictions", index=False)
    print(f"\nExcel: {OUT_XL}")
    print("Figures: external_validation_ROC_{14d,28d,14d_28d}.{png,pdf}")


if __name__ == "__main__":
    main()
