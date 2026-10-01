# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/auc_broad_grouping.py
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

"""AUC / ROC comparison of the four-gene signature under a broad case definition of
GSE32472 (0 = no BPD as control, 1/2/3 = any BPD as case).

Coefficients are fixed at the day-14 MS-BPD training-set values (never re-estimated).
The narrow definition (MS-BPD vs control) and the mild-only definition (mild vs
control) are also reported, to locate the source of the drop in AUC when the case
definition is broadened.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score

# ---------------- Paths ----------------
ROOT = P("data/raw/GSE32472/GSE32472_series_matrix.txt")
FULL = os.path.join(ROOT, "正常组与BPD比较")
XLSX = P("data/raw/GSE32472.xlsx")
OUT_IMG = P("figures/GSE32472_broad_grouping")
OUT_DAT = P("results/GSE32472_broad_grouping")
os.makedirs(OUT_IMG, exist_ok=True)
os.makedirs(OUT_DAT, exist_ok=True)

GENES = ["CD3G", "PRTN3", "CD6", "THBS1"]
COEF = {"(Intercept)": -3.375128, "CD3G": -0.819751,
        "PRTN3": 0.610007, "CD6": -1.402259, "THBS1": 1.560158}
# The two data sources specified by the user:
#   day-14 = 正常组与中重度BPD的比较/12-17d/GSE32472_12-17d.csv    (MS subset, n = 63)
#   day-28 = 正常组与中重度BPD的比较/22-27d/GSE32472_22-27d.xlsx   (sheet GSE32472_22-27dcomplete, n = 100)
USER_14D = P("data/raw/GSE32472/GSE32472_series_matrix.txt/正常组与中重度BPD的比较/12-17d/GSE32472_12-17d.csv")
USER_28D = P("data/raw/GSE32472/GSE32472_series_matrix.txt/正常组与中重度BPD的比较/22-27d/GSE32472_22-27d.xlsx")


def load_matrix(tag):
    """day-14: the full 97-sample matrix (same time point; values identical to the
       63 samples of the user-supplied csv);
       day-28: the complete sheet of the user-specified xlsx (100 samples)."""
    if tag == "14d":
        p = os.path.join(FULL, "12-17d", "GSE32472_12-17d.csv")
        df = pd.read_csv(p, index_col=0)
    else:
        df = pd.read_excel(USER_28D, sheet_name="GSE32472_22-27dcomplete", index_col=0)
    df.index = [str(x).strip() for x in df.index]
    return df

# Direction: CD3G/CD6 down in cases, PRTN3/THBS1 up in cases
# Single genes must be oriented so that a higher value means "more case-like"
GENE_SIGN = {"CD3G": -1, "CD6": -1, "PRTN3": 1, "THBS1": 1}


# ---------------- DeLong ----------------
def _midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T
    return T2


def _fast_delong(preds, m):
    """preds: (k, n) with the first m rows as the positive class. Returns aucs, cov."""
    n = preds.shape[1]
    k = preds.shape[0]
    tx = np.empty((k, m)); ty = np.empty((k, n - m))
    tz = np.empty((k, n))
    for r in range(k):
        tx[r] = _midrank(preds[r, :m])
        ty[r] = _midrank(preds[r, m:])
        tz[r] = _midrank(preds[r, :])
    aucs = tz[:, :m].mean(axis=1) / (n - m) - (m + 1.0) / (2.0 * (n - m))
    v01 = (tz[:, :m] - tx) / (n - m)
    v10 = 1.0 - (tz[:, m:] - ty) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    cov = sx / m + sy / (n - m)
    return aucs, np.atleast_2d(cov)


def delong_ci(y, proba):
    y = np.asarray(y); p = np.asarray(proba)
    order = np.argsort(-y)          # Positive class first
    preds = p[order][None, :]
    m = int((y == 1).sum())
    aucs, cov = _fast_delong(preds, m)
    se = np.sqrt(cov[0, 0])
    return float(aucs[0]), float(aucs[0] - 1.96 * se), float(aucs[0] + 1.96 * se), float(se)


def delong_test(y, p1, p2):
    """Paired DeLong test between two ROCs on the same set of samples."""
    y = np.asarray(y)
    order = np.argsort(-y)
    m = int((y == 1).sum())
    preds = np.vstack([np.asarray(p1)[order], np.asarray(p2)[order]])
    aucs, cov = _fast_delong(preds, m)
    var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    if var <= 0:
        return float(aucs[1] - aucs[0]), np.nan, np.nan
    z = (aucs[1] - aucs[0]) / np.sqrt(var)
    from scipy import stats
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return float(aucs[1] - aucs[0]), float(z), float(p)


def boot_ci(y, proba, n_boot=2000, seed=42):
    rng = np.random.default_rng(seed)
    y = np.asarray(y); p = np.asarray(proba)
    idx = np.arange(len(y))
    vals = []
    for _ in range(n_boot):
        b = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[b])) < 2:
            continue
        vals.append(roc_auc_score(y[b], p[b]))
    vals = np.array(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def youden(y, proba):
    fpr, tpr, thr = roc_curve(y, proba)
    j = tpr - fpr
    i = int(np.argmax(j))
    tn = int(((proba <= thr[i]) & (y == 0)).sum()); fn = int(((proba <= thr[i]) & (y == 1)).sum())
    fp = int(((proba > thr[i]) & (y == 0)).sum()); tp = int(((proba > thr[i]) & (y == 1)).sum())
    return float(thr[i]), tp / (tp + fn), tn / (tn + fp)


# ---------------- Data ----------------
def load_group():
    wb = pd.read_excel(XLSX, header=None, sheet_name=None)
    out = {}
    for sn, df in wb.items():
        tag = sn.replace("GSE32472_", "")
        ids = [str(v).strip() for v in df.iloc[0, 1:] if pd.notna(v)]
        grp = [str(v).strip() for v in df.iloc[1, 1:] if pd.notna(v)]
        grade = [int(g.split("(")[0].strip()) for g in grp]
        out[tag] = pd.DataFrame({"grade": grade}, index=ids)
    return out


def predict(vec):
    lp = COEF["(Intercept)"]
    for g in GENES:
        lp += COEF[g] * float(vec[g])
    return 1.0 / (1.0 + np.exp(-lp))


def plot_roc(ax, y, proba, label, color, ls="-", lw=2.2):
    fpr, tpr, _ = roc_curve(y, proba)
    ax.plot(fpr, tpr, color=color, lw=lw, ls=ls, label=label)


def main():
    groups = load_group()
    summary = []
    persample = []
    store = {}

    for tag in ["14d", "28d"]:
        mat = load_matrix(tag)
        grp = groups[tag]
        common = [s for s in grp.index if s in mat.columns]
        g = grp.loc[common]
        expr = mat.loc[GENES, common].astype(float).T          # samples x genes

        # Four-gene signature probability (fixed coefficients)
        proba = np.array([predict(expr.loc[s]) for s in common])
        grades = g["grade"].values

        store[tag] = dict(common=common, expr=expr, proba=proba, grades=grades)

        # ---- Three case definitions ----
        scenarios = {
            "wide_any_BPD":      (np.where(grades >= 1, 1, 0), "Any BPD vs no BPD"),
            "narrow_MS_BPD":     (np.where(grades == 0, 0, np.where(grades >= 2, 1, np.nan)), "MS-BPD vs no BPD"),
            "mild_only":         (np.where(grades == 0, 0, np.where(grades == 1, 1, np.nan)), "Mild BPD vs no BPD"),
            "modsev_only":       (np.where(grades == 0, 0, np.where(grades >= 2, 1, np.nan)), "Moderate/severe BPD vs no BPD"),
        }
        for key, (y, nice) in scenarios.items():
            mask = ~np.isnan(y)
            yy = y[mask].astype(int); pp = proba[mask]
            if len(np.unique(yy)) < 2:
                continue
            auc, lo, hi, se = delong_ci(yy, pp)
            blo, bhi = boot_ci(yy, pp)
            thr, sens, spec = youden(yy, pp)
            summary.append(dict(
                timepoint=tag, definition=key, definition_note=nice,
                n_cases=int(yy.sum()), n_controls=int((yy == 0).sum()), total=len(yy),
                AUC=round(auc, 4), AUC_95CI_DeLong="%.3f-%.3f" % (lo, hi),
                AUC_95CI_bootstrap="%.3f-%.3f" % (blo, bhi),
                youden_cutoff=round(thr, 4), sensitivity=round(sens, 4), specificity=round(spec, 4)))
            # Per-sample predictions only for the two main definitions
            if key in ("wide_any_BPD", "narrow_MS_BPD"):
                for s, yv, pv, gv in zip(np.array(common)[mask], yy, pp, grades[mask]):
                    persample.append(dict(timepoint=tag, definition=key, sample=s, bpd_grade=int(gv),
                                          true_label=int(yv), predicted_probability=round(float(pv), 4)))

        # ---- Single-gene AUC (broad definition, orientation corrected) ----
        yw = np.where(grades >= 1, 1, 0)
        for gene in GENES:
            score = GENE_SIGN[gene] * expr[gene].values
            auc = roc_auc_score(yw, score)
            _, lo, hi, _ = delong_ci(yw, score)
            summary.append(dict(
                timepoint=tag, definition="single_gene_" + gene, definition_note="%s alone (any BPD vs no BPD)" % gene,
                n_cases=int(yw.sum()), n_controls=int((yw == 0).sum()), total=len(yw),
                AUC=round(auc, 4), AUC_95CI_DeLong="%.3f-%.3f" % (lo, hi),
                AUC_95CI_bootstrap="", youden_cutoff="", sensitivity="", specificity=""))

        # ---- Control: apparent AUC from refitting on the same samples (not validation) ----
        import statsmodels.api as sm
        for key in ("narrow_MS_BPD", "wide_any_BPD"):
            y = scenarios_apply(grades, key)
            mask = ~np.isnan(y); yy = y[mask].astype(int)
            X = sm.add_constant(expr.loc[np.array(common)[mask], GENES].astype(float).values)
            res = sm.Logit(yy, X).fit(disp=0)
            pr = res.predict(X)
            summary.append(dict(
                timepoint=tag, definition="REFIT_apparent_" + key,
                definition_note="Refitted on the same samples (apparent AUC, not validation)",
                n_cases=int(yy.sum()), n_controls=int((yy == 0).sum()), total=len(yy),
                AUC=round(float(roc_auc_score(yy, pr)), 4), AUC_95CI_DeLong="", AUC_95CI_bootstrap="",
                youden_cutoff="intercept=%.2f" % res.params[0],
                sensitivity="CD3G=%.2f" % res.params[1], specificity="THBS1=%.2f" % res.params[4]))

        # ---- Paired DeLong: signature vs each single gene (same samples, valid) ----
        yw_ = np.where(grades >= 1, 1, 0)
        for gene in GENES:
            score = GENE_SIGN[gene] * expr[gene].values
            d, z, p = delong_test(yw_, score, proba)
            summary.append(dict(
                timepoint=tag, definition="delong_sig_vs_" + gene,
                definition_note="AUC(signature)-AUC(%s), paired DeLong" % gene,
                n_cases=int(yw_.sum()), n_controls=int((yw_ == 0).sum()), total=len(yw_),
                AUC=round(d, 4), AUC_95CI_DeLong="", AUC_95CI_bootstrap="",
                youden_cutoff="z=%.3f" % z if not np.isnan(z) else "z=NA",
                sensitivity="P=%.4f" % p if not np.isnan(p) else "P=NA", specificity=""))

    # ---------------- Figures ----------------
    COL_W = "#C0392B"   # Broad definition (red)
    COL_N = "#1F4E79"   # Narrow definition (blue)
    COL_M = "#27AE60"   # Mild only (green)

    # Figure 1: broad definition, 14d vs 28d
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.8))
    fig.patch.set_facecolor("white")
    for ax, tag in zip(axes, ["14d", "28d"]):
        st = store[tag]
        y = np.where(st["grades"] >= 1, 1, 0)
        row = [r for r in summary if r["timepoint"] == tag and r["definition"] == "wide_any_BPD"][0]
        plot_roc(ax, y, st["proba"], "Four-gene signature (AUC=%.3f, 95%%CI %s)" % (row["AUC"], row["AUC_95CI_DeLong"]), COL_W)
        ax.plot([0, 1], [0, 1], ls=":", color="#BBBBBB", lw=1)
        ax.set_xlabel("1 - Specificity (False Positive Rate)", fontsize=11)
        ax.set_ylabel("Sensitivity (True Positive Rate)", fontsize=11)
        ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.set_title("GSE32472 %s  (cases=%d, controls=%d)" % (tag, int(y.sum()), int((y == 0).sum())), fontsize=12)
        ax.legend(loc="lower right", frameon=False, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_IMG, "ROC_broad_anyBPD_14d_28d.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUT_IMG, "ROC_broad_anyBPD_14d_28d.pdf"), bbox_inches="tight")
    plt.close(fig)

    # Figure 2: the three definitions side by side (one panel per time point)
    for tag in ["14d", "28d"]:
        st = store[tag]; grades = st["grades"]
        fig, ax = plt.subplots(figsize=(6.8, 6.0)); fig.patch.set_facecolor("white")
        for key, col, nice, ls in [("wide_any_BPD", COL_W, "Any BPD", "-"),
                                   ("modsev_only", COL_N, "Moderate/severe BPD", "-"),
                                   ("mild_only", COL_M, "Mild BPD", "--")]:
            y = scenarios_apply(grades, key)
            mask = ~np.isnan(y); yy = y[mask].astype(int); pp = st["proba"][mask]
            r = [x for x in summary if x["timepoint"] == tag and x["definition"] == key][0]
            plot_roc(ax, yy, pp, "%s (n=%d vs %d, AUC=%.3f)" % (nice, int(yy.sum()), int((yy == 0).sum()), r["AUC"]), col, ls=ls)
        ax.plot([0, 1], [0, 1], ls=":", color="#BBBBBB", lw=1)
        ax.set_xlabel("1 - Specificity", fontsize=11); ax.set_ylabel("Sensitivity", fontsize=11)
        ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.set_title("GSE32472 %s: effect of case definition" % tag, fontsize=12)
        ax.legend(loc="lower right", frameon=False, fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_IMG, "ROC_definition_comparison_%s.png" % tag), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(OUT_IMG, "ROC_definition_comparison_%s.pdf" % tag), bbox_inches="tight")
        plt.close(fig)

    # Figure 3: ROC of each single gene (broad definition)
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.8)); fig.patch.set_facecolor("white")
    palette = {"CD3G": "#1F4E79", "PRTN3": "#C0392B", "CD6": "#27AE60", "THBS1": "#8E44AD"}
    for ax, tag in zip(axes, ["14d", "28d"]):
        st = store[tag]
        y = np.where(st["grades"] >= 1, 1, 0)
        plot_roc(ax, y, st["proba"], "Four-gene signature (AUC=%.3f)" % roc_auc_score(y, st["proba"]), "#111111", lw=2.8)
        for gene in GENES:
            score = GENE_SIGN[gene] * st["expr"][gene].values
            plot_roc(ax, y, score, "%s (AUC=%.3f)" % (gene, roc_auc_score(y, score)), palette[gene], lw=1.8)
        ax.plot([0, 1], [0, 1], ls=":", color="#BBBBBB", lw=1)
        ax.set_xlabel("1 - Specificity", fontsize=11); ax.set_ylabel("Sensitivity", fontsize=11)
        ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.set_title("GSE32472 %s: signature vs single genes (any BPD)" % tag, fontsize=12)
        ax.legend(loc="lower right", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_IMG, "ROC_single_gene_broad.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUT_IMG, "ROC_single_gene_broad.pdf"), bbox_inches="tight")
    plt.close(fig)

    # ---------------- Export ----------------
    sdf = pd.DataFrame(summary)
    with pd.ExcelWriter(os.path.join(OUT_DAT, "GSE32472_broad_grouping_four_gene_AUC.xlsx")) as w:
        sdf.to_excel(w, sheet_name="Summary", index=False)
        pd.DataFrame(persample).to_excel(w, sheet_name="Per-sample predictions", index=False)

    sdf.to_csv(os.path.join(OUT_DAT, "GSE32472_broad_grouping_summary.csv"), index=False, encoding="utf-8-sig")
    print(sdf.to_string(index=False))
    print()
    print("Figures:", sorted(os.listdir(OUT_IMG)))
    print("Data:", sorted(os.listdir(OUT_DAT)))
    return sdf


def scenarios_apply(grades, key):
    if key == "wide_any_BPD":
        return np.where(grades >= 1, 1, 0).astype(float)
    if key == "narrow_MS_BPD":
        return np.where(grades == 0, 0, np.where(grades >= 2, 1, np.nan))
    if key == "mild_only":
        return np.where(grades == 0, 0, np.where(grades == 1, 1, np.nan))
    if key == "modsev_only":
        return np.where(grades == 0, 0, np.where(grades >= 2, 1, np.nan))
    raise KeyError(key)


if __name__ == "__main__":
    main()
