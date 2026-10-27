# -*- coding: utf-8 -*-
"""
Supplementary check run on the authors' own processed GSE32472 day-14 matrix
(Supplementary Data S1) + per-sample labels (Supplementary Data S14).

Purpose (peer-review response):
  R1 major #4 : do lung-derived cells contribute to the whole-blood signal?
  R2 major #3 : is the deconvolution-inferred myeloid/T-cell pattern supported
                by an orthogonal, label-free marker-based analysis?

Outputs: marker_check_results.csv, marker_check_report.md
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))
SUPPL = r"C:\Users\44937\Desktop\投稿内容\送审包_20260911\Supplementary_FrontPediatr"

expr = pd.read_excel(os.path.join(SUPPL, "Data S1.xlsx"), sheet_name="S1")
expr = expr.set_index("ID")
expr.index = expr.index.astype(str)

labels = pd.read_excel(os.path.join(SUPPL, "Data S14.xlsx"),
                       sheet_name="Per-sample predictions")
labels = labels[["sample", "group", "fixed_score"]].dropna()

# ---- sanity check: reproduce the published fixed-coefficient score ----
COEF = {"PRTN3": 0.610, "CD3G": -0.820, "CD6": -1.402, "THBS1": 1.560}
INTERCEPT = -3.375
common = [s for s in labels["sample"] if s in expr.columns]
lab = labels.set_index("sample").loc[common]
X = expr[common]
score = INTERCEPT
for g, c in COEF.items():
    if g in expr.index:
        score = score + c * X.loc[g].values
    else:
        raise SystemExit("gene missing from matrix: " + g)
r = np.corrcoef(score, lab["fixed_score"].values)[0, 1]

y = (lab["group"].str.lower().str.contains("bpd")).astype(int).values
case_mask = y == 1
ctrl_mask = y == 0
n_case, n_ctrl = int(case_mask.sum()), int(ctrl_mask.sum())

# background level: median expression of all genes across all samples
bg_median = float(np.median(X.values))
bg_q75 = float(np.percentile(X.values, 75))

SETS = {
    "Lung epithelium (alveolar)": ["SFTPB", "SFTPC", "SFTPA1", "SFTPA2", "AGER", "AQP5"],
    "Lung epithelium (airway)": ["SCGB1A1", "SCGB3A2", "FOXJ1", "TP63"],
    "Epithelium (pan)": ["EPCAM", "KRT8", "KRT18", "CDH1"],
    "Endothelium": ["PECAM1", "CDH5", "CLDN5", "VWF", "TEK"],
    "Mesenchyme / smooth muscle": ["ACTA2", "MYH11", "COL1A1", "PDGFRB"],
    "Immature granulocyte / emergency granulopoiesis":
        ["PRTN3", "ELANE", "AZU1", "CTSG", "DEFA4", "CAMP", "LCN2", "LTF", "CEACAM8", "BPI"],
    "Mature neutrophil": ["FCGR3B", "CXCR2", "CSF3R", "SELL", "ITGAM", "MPO"],
    "T cell": ["CD3D", "CD3E", "CD3G", "CD6", "CD2", "IL7R", "LCK", "TRAC"],
    "Platelet": ["PPBP", "PF4", "ITGA2B", "GP9"],
    "Signature genes": ["CD3G", "PRTN3", "CD6", "THBS1"],
}


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


rows = []
for setname, genes in SETS.items():
    present = [g for g in genes if g in expr.index]
    absent = [g for g in genes if g not in expr.index]
    # per-gene test
    for g in present:
        a = X.loc[g].values[case_mask].astype(float)
        b = X.loc[g].values[ctrl_mask].astype(float)
        t, p = stats.ttest_ind(a, b, equal_var=False)
        rows.append({
            "set": setname, "gene": g, "level": "gene",
            "mean_all": float(X.loc[g].values.mean()),
            "mean_MSBPD": float(a.mean()), "mean_control": float(b.mean()),
            "log2FC_MSBPD_vs_control": float(a.mean() - b.mean()),
            "t": float(t), "P": float(p),
            "above_background_median": bool(X.loc[g].values.mean() > bg_median),
        })
    # set-level signature score (z-scored across samples, then averaged)
    if present:
        Z = X.loc[present].apply(lambda r_: (r_ - r_.mean()) / r_.std(ddof=1), axis=1)
        s = Z.mean(axis=0).values
        t, p = stats.ttest_ind(s[case_mask], s[ctrl_mask], equal_var=False)
        rows.append({
            "set": setname, "gene": "|".join(present), "level": "signature score",
            "mean_all": np.nan, "mean_MSBPD": float(s[case_mask].mean()),
            "mean_control": float(s[ctrl_mask].mean()),
            "log2FC_MSBPD_vs_control": float(s[case_mask].mean() - s[ctrl_mask].mean()),
            "t": float(t), "P": float(p), "above_background_median": np.nan,
        })
    if absent:
        rows.append({"set": setname, "gene": "|".join(absent), "level": "absent from platform",
                     "mean_all": np.nan, "mean_MSBPD": np.nan, "mean_control": np.nan,
                     "log2FC_MSBPD_vs_control": np.nan, "t": np.nan, "P": np.nan,
                     "above_background_median": np.nan})

res = pd.DataFrame(rows)
tested = res["P"].notna()
res.loc[tested, "FDR_BH"] = bh(res.loc[tested, "P"].values)
res = res.sort_values(["set", "level", "P"], na_position="last")
res.to_csv(os.path.join(BASE, "marker_check_results.csv"), index=False, encoding="utf-8-sig")

sig = res[res["level"] == "signature score"].copy()

md = []
md.append("# Whole-blood marker check on the authors' own GSE32472 day-14 matrix\n")
md.append("Data: Supplementary Data S1 (processed expression) + Supplementary Data S14 (labels).")
md.append(f"Samples: {n_case} MS-BPD vs {n_ctrl} controls (n = {n_case + n_ctrl}).\n")
md.append(f"**Validation of data handling:** score recomputed with the published fixed "
          f"coefficients correlates with the deposited fixed_score at r = {r:.4f} "
          f"(n = {len(common)}), confirming that the matrix and labels were read correctly.\n")
md.append(f"Array-wide background: median log2 expression = {bg_median:.2f}, "
          f"75th percentile = {bg_q75:.2f}.\n")
md.append("## Set-level signature scores (z-scored genes, averaged per sample)\n")
md.append("| Gene set | Mean z (MS-BPD) | Mean z (control) | Delta | t | P | FDR (BH) |")
md.append("|---|---|---|---|---|---|---|")
for _, row in sig.iterrows():
    md.append(f"| {row['set']} | {row['mean_MSBPD']:.3f} | {row['mean_control']:.3f} | "
              f"{row['log2FC_MSBPD_vs_control']:+.3f} | {row['t']:.2f} | {row['P']:.2e} | "
              f"{row['FDR_BH']:.3f} |")
md.append("\n## Per-gene detail\n")
md.append("| Gene set | Gene | Mean log2 (all) | MS-BPD | Control | log2FC | P | FDR (BH) | Above array median |")
md.append("|---|---|---|---|---|---|---|---|---|")
for _, row in res[res["level"] == "gene"].iterrows():
    md.append(f"| {row['set']} | {row['gene']} | {row['mean_all']:.2f} | {row['mean_MSBPD']:.2f} | "
              f"{row['mean_control']:.2f} | {row['log2FC_MSBPD_vs_control']:+.3f} | {row['P']:.2e} | "
              f"{row['FDR_BH']:.3f} | {row['above_background_median']} |")
absent_rows = res[res["level"] == "absent from platform"]
if len(absent_rows):
    md.append("\n## Genes not represented on GPL6244\n")
    for _, row in absent_rows.iterrows():
        md.append(f"- {row['set']}: {row['gene']}")
md.append("\n## Interpretation (layered, as it should be reported)\n")
md.append("### 1. Deconvolution-free confirmation of the immune pattern (strongest result)")
md.append("Every T-cell-restricted transcript present on the platform moved in the same direction "
          "(CD2, CD3D, CD3E, CD3G, CD6, LCK; all FDR < 1e-5), and the granulocyte set moved in the "
          "opposite direction (CSF3R FDR = 4.5e-4, ITGAM FDR = 1.6e-3, MPO FDR = 0.034). Because "
          "this analysis uses only raw marker-gene levels and does **not** run CIBERSORT or any "
          "deconvolution, it is an orthogonal check: the myeloid-versus-T-cell pattern is visible "
          "in the marker genes themselves and is not an artefact of one algorithm.")
md.append("### 2. Platelet contribution")
md.append("ITGA2B rose by +0.51 log2 (FDR = 0.012), consistent with the single-cell localisation of "
          "THBS1 to platelets. This also provides the most parsimonious explanation for the mild "
          "rise in endothelial markers (see 3).")
md.append("### 3. Lung-derived cells — what this dataset can and cannot answer")
absent_lung = res[(res["set"].str.contains("Lung epithelium")) &
                  (res["level"] == "absent from platform")]["gene"].tolist()
md.append(f"**Important limitation:** the most lung-restricted transcripts "
          f"({' / '.join(absent_lung)}) are **not represented on the GPL6244 array**, so their "
          "contribution cannot be measured in this cohort at all. The lung-related probes that do "
          "exist (AQP5, SCGB1A1, FOXJ1) sit below the array-wide median and showed no "
          "between-group difference (all P > 0.27, all FDR > 0.35).")
md.append("Endothelial markers rose slightly (VWF +0.14, P = 0.044; CDH5 +0.09, P = 0.056; "
          "PECAM1 +0.11, P = 0.27), but PECAM1, CDH5 and VWF are also expressed by platelets, "
          "monocytes and endothelial progenitor cells, and the platelet set rose in parallel. The "
          "parsimonious reading is therefore a circulating platelet / progenitor signal rather "
          "than shed pulmonary endothelium.")
md.append("**Honest bottom line:** these data do not exclude a lung-derived contribution, and we "
          "should not claim that they do. What they show is that the whole-blood signal is "
          "dominated by circulating immune and platelet populations, and that the four-gene "
          "signature and the deconvolution pattern are reproduced by marker genes alone. "
          "Definitive exclusion of a lung-derived contribution requires paired whole-blood and "
          "cell-fraction (or cell-free RNA) samples, which we state as a limitation.")
md.append("### 4. PRTN3 and immature granulocytes")
md.append("PRTN3 itself rose (+0.61 log2, FDR = 0.015), with ELANE (+0.62, FDR = 0.063) and AZU1 "
          "(+0.44, FDR = 0.067) in the same direction. The set-level score did not pass FDR < 0.05 "
          "(P = 0.034, FDR = 0.063). The direction is consistent with granulopoiesis-driven "
          "release, but the evidence is suggestive rather than conclusive.")
md.append("\n### Statistics")
md.append(f"Welch t-test per gene; Benjamini-Hochberg FDR within the tested pool "
          f"({int(tested.sum())} tests). Set-level scores are the mean of per-gene z-scores across "
          "the 63 samples.")

with open(os.path.join(BASE, "marker_check_report.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(md))

print("sanity r vs deposited fixed_score:", round(r, 4))
print("n case / ctrl:", n_case, n_ctrl)
print("array median:", round(bg_median, 2))
print(sig[["set", "mean_MSBPD", "mean_control", "log2FC_MSBPD_vs_control", "P", "FDR_BH"]].to_string(index=False))
