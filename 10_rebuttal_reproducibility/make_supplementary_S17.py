# -*- coding: utf-8 -*-
"""Assemble Supplementary Data S17 for the revision (R1 major #3 response)."""
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ANA = os.path.join(os.path.dirname(HERE), "02_补充分析")
OUT = os.path.join(HERE, "Supplementary Data S17.xlsx")

sheets = {
    "S17-1 Pipeline by definition": "irg_sensitivity_v3_summary.csv",
    "S17-2 Bootstrap selection freq": "irg_sensitivity_v3_bootstrap.csv",
    "S17-3 CytoHubba cut-off sweep": "irg_sensitivity_v3_module_sweep.csv",
    "S17-4 Prediction stability": "irg_sensitivity_v3_prediction_stability.csv",
}

NOTE = [
    ("Title",
     "Sensitivity of the four-gene signature to the definition of the "
     "inflammatory-response gene set (response to Reviewer 1, major comment 3)"),
    ("Data", "GSE32472, postnatal day 14 training set; 24 MS-BPD vs 39 controls "
             "(Supplementary Data S1, labels in Supplementary Data S14)"),
    ("DEG pool",
     "limma moderated t (R, limma 3.54); P < 0.05 and |log2FC| > 0.58. "
     "This reproduces the numbers reported in the manuscript exactly: "
     "139 DEGs (84 up / 55 down), 57 DE-IRGs."),
    ("Screening chain",
     "DEG x IRG -> STRING PPI (>0.4, Supplementary Data S9) -> CytoHubba "
     "MCC/EPC/Degree module -> LASSO (10-fold CV) intersected with SVM-RFE (8 features)"),
    ("Implementation",
     "Python re-implementation (igraph + scikit-learn) of the manuscript's "
     "Cytoscape + R (glmnet/e1071) pipeline. The intermediate quantities match "
     "the published ones (57 DE-IRGs; 48 nodes / 172 interactions); the final "
     "intersection is obtained for CD3G, CD6 and THBS1, while PRTN3 - the most "
     "marginal of the four at the DEG stage (log2FC 0.61, FDR 0.039, rank "
     "1178/9323) - is not recovered by the re-implementation at any cut-off."),
    ("Definitions",
     "Def-1 GeneCards 'inflammatory response' union MSigDB human collection (as used in the "
     "manuscript); Def-2 GO:0006954 with descendant terms; Def-3 MSigDB "
     "HALLMARK_INFLAMMATORY_RESPONSE; Def-4 union of GO:0006954/0006955/0045087/0002526; "
     "Def-5 Def-4 plus eight inflammation-related Hallmark sets"),
    ("Sheets",
     "S17-1 full-data run per definition; S17-2 selection frequency in 100 bootstrap "
     "resamples with the DEG pool held fixed; S17-3 sensitivity of the result to the "
     "CytoHubba top-N cut-off; S17-4 bootstrap confidence intervals for the AUC of the "
     "locked-coefficient score and of the cross-validated prediction"),
    ("Caveat",
     "The definitions differ in curation philosophy (manual vs experimental vs "
     "computational), so differences between them are not attributable to the search "
     "term alone. Selection frequencies quantify resampling variance at n=63; they are "
     "not a test of the signature's predictive validity, which rests on the "
     "locked-coefficient transport validation reported in the manuscript."),
]

with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    pd.DataFrame(NOTE, columns=["Item", "Description"]).to_excel(
        w, sheet_name="S17-0 Notes", index=False)
    for sheet, fn in sheets.items():
        p = os.path.join(ANA, fn)
        if os.path.exists(p):
            pd.read_csv(p, encoding="utf-8-sig").to_excel(w, sheet_name=sheet, index=False)
            print("  +", sheet, "<-", fn)
        else:
            print("  ! missing", fn)
    for ws in w.book.worksheets:
        for col in ws.columns:
            width = max(len(str(c.value)) for c in col if c.value is not None)
            ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 10), 70)

print("written ->", OUT)
