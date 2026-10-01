# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/roc_diagnosis_4v3.py
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

"""14d / 28d ROC curves: four-gene vs three-gene diagnostic models
- features: DESeq2 vst values (same source as the group-comparison figure)
- model: logistic regression (StandardScaler fitted inside the training fold, to avoid leakage)
- evaluation: single-gene AUC / four-gene model / four three-gene combinations; training AUC and LOOCV AUC
- comparison: DeLong test (four-gene vs best three-gene model)
- output: ROC figures (PNG + PDF) and a summary Excel workbook
"""
import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, roc_auc_score
from pydeseq2.dds import DeseqDataSet

GENES = ["PRTN3", "CD3G", "CD6", "THBS1"]
CSVs = {
    "14d": P("data/raw/GSE220135/GSE220135_14d.csv"),
    "28d": P("data/raw/GSE220135/GSE220135_28d.csv"),
}
OUT_IMG = P("figures/GSE220135")
OUT_XL = P("results/GSE220135/ROC_diagnostic_model_comparison_14d_28d.xlsx")
os.makedirs(OUT_IMG, exist_ok=True)
os.makedirs(os.path.dirname(OUT_XL), exist_ok=True)

C4 = "#E15759"   # four genes (red, consistent with the BPD red used earlier)
C3 = "#4FA8C5"   # three genes (blue, consistent with the Normal blue used earlier)


def parse_group(col):
    return "nonBPD" if "NONBPD" in col.upper() else "BPD"


def get_vst(csv_path):
    counts = pd.read_csv(csv_path, index_col=0).astype(int).T
    grp = {c: parse_group(c) for c in counts.index}
    meta = pd.DataFrame({"condition": [grp[c] for c in counts.index]},
                        index=counts.index)
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition")
    dds.deseq2()
    dds.vst_fit()
    arr = dds.vst_transform()
    vst = pd.DataFrame(arr, index=counts.index, columns=counts.columns)
    return vst, grp   # samples (rows) x genes (columns)


def loocv_proba(X, y):
    """Leave-one-out cross-validation: fit StandardScaler + LR within each fold; return per-sample probabilities."""
    n = len(y)
    proba = np.zeros(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        sc = StandardScaler().fit(X[mask])
        clf = LogisticRegression(max_iter=1000).fit(sc.transform(X[mask]), y[mask])
        proba[i] = clf.predict_proba(sc.transform(X[i:i+1]))[:, 1][0]
    return proba


def delong_test(y, p1, p2):
    """DeLong test for the AUC difference between two correlated ROCs (placements method, ties handled)."""
    y = np.asarray(y); p1 = np.asarray(p1); p2 = np.asarray(p2)
    pos = y == 1; neg = y == 0
    m = int(pos.sum()); n = int(neg.sum())

    def placement(pred):
        v10 = np.zeros(m); v01 = np.zeros(n)
        xp, xn = pred[pos], pred[neg]
        for i in range(m):
            v10[i] = (np.sum(xp[i] > xn) + 0.5 * np.sum(xp[i] == xn)) / n
        for j in range(n):
            v01[j] = (np.sum(xp > xn[j]) + 0.5 * np.sum(xp == xn[j])) / m
        auc = v10.mean()
        var = np.var(v10, ddof=1) / m + np.var(v01, ddof=1) / n
        return auc, v10, v01, var

    a1, v10_1, v01_1, var1 = placement(p1)
    a2, v10_2, v01_2, var2 = placement(p2)
    cov = np.cov(v10_1, v10_2, ddof=1)[0, 1] / m + np.cov(v01_1, v01_2, ddof=1)[0, 1] / n
    se = max(var1 + var2 - 2 * cov, 1e-12)
    z = (a1 - a2) / np.sqrt(se)
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return p


def plot_roc(ax, y, p4, p3, auc4, auc3, title, tag):
    for proba, c, lbl in [(p4, C4, "4-gene"), (p3, C3, "3-gene")]:
        fpr, tpr, _ = roc_curve(y, proba)
        ax.plot(fpr, tpr, color=c, lw=2.2,
                label=f"{lbl} (AUC={roc_auc_score(y, proba):.3f})")
    ax.plot([0, 1], [0, 1], ls="--", color="#888888", lw=1)
    ax.set_xlabel("False Positive Rate", fontsize=11)
    ax.set_ylabel("True Positive Rate", fontsize=11)
    ax.set_title(f"{title}", fontsize=13, fontweight="bold", loc="left")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=10)


def main():
    rows = []
    for tag in ["14d", "28d"]:
        vst, grp = get_vst(CSVs[tag])
        X = vst[GENES].values
        y = np.array([1 if grp[s] == "BPD" else 0 for s in vst.index])

        # single-gene AUC (using vst values directly)
        print(f"\n===== {tag}: BPD={y.sum()} nonBPD={len(y)-y.sum()} =====")
        single = {}
        for g in GENES:
            single[g] = roc_auc_score(y, X[:, GENES.index(g)])
            print(f"  single gene {g:6s} AUC={single[g]:.3f}")

        # combined models: four genes + four three-gene combinations
        combos = {"4-gene": GENES}
        for c3 in itertools.combinations(GENES, 3):
            combos["+".join(c3)] = list(c3)

        res = {}
        for name, feats in combos.items():
            Xc = X[:, [GENES.index(f) for f in feats]]
            p_loocv = loocv_proba(Xc, y)
            sc = StandardScaler().fit(Xc)
            clf = LogisticRegression(max_iter=1000).fit(sc.transform(Xc), y)
            p_train = clf.predict_proba(sc.transform(Xc))[:, 1]
            auc_tr = roc_auc_score(y, p_train)
            auc_cv = roc_auc_score(y, p_loocv)
            res[name] = dict(train_auc=auc_tr, cv_auc=auc_cv, proba=p_loocv)
            print(f"  {name:22s} training_auc={auc_tr:.3f}  LOOCV AUC={auc_cv:.3f}")

        best3 = max([k for k in res if k != "4-gene"], key=lambda k: res[k]["cv_auc"])
        p_delong = delong_test(y, res["4-gene"]["proba"], res[best3]["proba"])
        print(f"  best three genes: {best3} (LOOCV AUC={res[best3]['cv_auc']:.3f})")
        print(f"  four genes vs {best3}: DeLong p={p_delong:.4f}")

        # plotting: single figure
        fig, ax = plt.subplots(figsize=(6.2, 5.6))
        fig.patch.set_facecolor("white")
        plot_roc(ax, y, res["4-gene"]["proba"], res[best3]["proba"],
                 res["4-gene"]["cv_auc"], res[best3]["cv_auc"],
                 f"{tag}", tag)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_IMG, f"ROC_{tag}.png"), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(OUT_IMG, f"ROC_{tag}.pdf"), bbox_inches="tight")
        plt.close(fig)

        # summary table
        for g in GENES:
            rows.append(dict(timepoint=tag, model=f"single gene-{g}", n_genes=1,
                             training_auc=single[g], LOOCV_AUC=single[g],
                             note="single genes use vst values directly; training = LOOCV"))
        for name, r in res.items():
            rows.append(dict(timepoint=tag, model=name, n_genes=4 if name == "4-gene" else 3,
                             training_auc=r["train_auc"], LOOCV_AUC=r["cv_auc"],
                             note="best three genes" if name == best3 else ("four genes" if name == "4-gene" else "")))
        rows.append(dict(timepoint=tag, model=f"DeLong: four genes vs {best3}", n_genes=np.nan,
                         training_auc=np.nan, LOOCV_AUC=np.nan, note=f"p={p_delong:.4f}"))

        # cached data for the combined figure (overwritten on the last round; used for the 14d_28d figure)
        if tag == "14d":
            cache14 = (y, res["4-gene"]["proba"], res[best3]["proba"], best3)
        else:
            cache28 = (y, res["4-gene"]["proba"], res[best3]["proba"], best3)

    # combined figure 14d+28d
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6))
    fig.patch.set_facecolor("white")
    for ax, (tag, cache) in zip(axes, [("14d", cache14), ("28d", cache28)]):
        y, p4, p3, best3 = cache
        plot_roc(ax, y, p4, p3, None, None, f"{tag}", tag)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_IMG, "ROC_14d_28d.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUT_IMG, "ROC_14d_28d.pdf"), bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(rows).to_excel(OUT_XL, index=False)
    print(f"\nExcel: {OUT_XL}")
    print("Figures: ROC_14d.png/pdf, ROC_28d.png/pdf, ROC_14d_28d.png/pdf")


if __name__ == "__main__":
    main()
