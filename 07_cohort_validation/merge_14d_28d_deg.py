# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/merge_14d_28d_deg.py
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

"""Merge the differential-expression results for all genes at both GSE220135 time
points (Day14 and Day28) into a single Excel workbook.
Output: merged all-gene table (time points side by side) + significant DEGs per time
point + shared DEGs + genes with a direction change
"""
import pandas as pd
import os

BASE = P("results/GSE220135")
F14 = os.path.join(BASE, "GSE220135_14d_DEG_results.xlsx")
F28 = os.path.join(BASE, "GSE220135_28d_DEG_results.xlsx")
OUT = os.path.join(BASE, "GSE220135_14d_28d_merged_DEG.xlsx")

stat_cols = ["baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj_FDR", "gene_type"]

def load(tag):
    df = pd.read_excel((F14 if tag == "14d" else F28), sheet_name="All genes")
    df = df[["Gene"] + stat_cols].copy()
    df.columns = ["Gene"] + [f"{c}_{tag}" for c in stat_cols]
    return df

d14 = load("14d")
d28 = load("28d")

# Merge (outer join, in case the two gene sets are not fully identical)
m = pd.merge(d14, d28, on="Gene", how="outer")

# Significance flag (padj < 0.05 and |log2FC| > 1)
m["sig_14d"] = (m["padj_FDR_14d"] < 0.05) & (m["log2FC_14d"].abs() > 1)
m["sig_28d"] = (m["padj_FDR_28d"] < 0.05) & (m["log2FC_28d"].abs() > 1)

# Direction concordance (only for genes significant at both time points)
def direction(fc):
    if pd.isna(fc):
        return "NA"
    return "up" if fc > 0 else "down"

m["dir_14d"] = m["log2FC_14d"].apply(direction)
m["dir_28d"] = m["log2FC_28d"].apply(direction)

both_sig = m["sig_14d"] & m["sig_28d"]
m["Concordant across time points"] = ""
m.loc[both_sig, "Concordant across time points"] = m.loc[both_sig].apply(
    lambda r: "concordant" if r["dir_14d"] == r["dir_28d"] else "reversed", axis=1)

# Sort: significant genes first, then by the smaller padj of the two time points
m["_minpadj"] = m[["padj_FDR_14d", "padj_FDR_28d"]].min(axis=1)
m = m.sort_values("_minpadj", na_position="last").drop(columns="_minpadj")

# ---- Sub-tables ----
sig14 = m[m["sig_14d"]].copy().sort_values("padj_FDR_14d")
sig28 = m[m["sig_28d"]].copy().sort_values("padj_FDR_28d")
shared = m[both_sig].copy()
shared = shared.sort_values(["Concordant across time points", "padj_FDR_14d"], ascending=[True, True])
flip = shared[shared["Concordant across time points"] == "reversed"].copy()

# ---- Summary sheet ----
n14 = int(m["sig_14d"].sum())
n28 = int(m["sig_28d"].sum())
summary = pd.DataFrame({
    "Item": [
        "Data source", "Total genes (rows after merge)",
        "Day14 samples (BPD/nonBPD)", "Day28 samples (BPD/nonBPD)",
        "Day14 significant DEGs (padj<0.05 & |log2FC|>1)",
        "  - up-regulated (log2FC>0)", "  - down-regulated (log2FC<0)",
        "Day28 significant DEGs (padj<0.05 & |log2FC|>1)",
        "  - up-regulated (log2FC>0)", "  - down-regulated (log2FC<0)",
        "Significant at both time points", "  - concordant", "  - reversed",
        "Day14 only", "Day28 only",
        "DE method", "Significance threshold",
    ],
    "Value": [
        "GSE220135 (BPD lung tissue vs nonBPD)", f"{len(m)}",
        "14 / 54", "13 / 55",
        f"{n14}", f"{int((sig14['log2FC_14d']>0).sum())}", f"{int((sig14['log2FC_14d']<0).sum())}",
        f"{n28}", f"{int((sig28['log2FC_28d']>0).sum())}", f"{int((sig28['log2FC_28d']<0).sum())}",
        f"{len(shared)}", f"{int((shared['Concordant across time points']=='concordant').sum())}", f"{int((shared['Concordant across time points']=='reversed').sum())}",
        f"{int((m['sig_14d'] & ~m['sig_28d']).sum())}", f"{int((~m['sig_14d'] & m['sig_28d']).sum())}",
        "DESeq2 (pydeseq2 0.5.4, negative-binomial Wald test + BH correction)",
        "padj_FDR < 0.05 and |log2FC| > 1",
    ],
})

with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
    m.to_excel(xw, sheet_name="All genes (14d+28d)", index=False)
    sig14.to_excel(xw, sheet_name="Day14 significant DEGs", index=False)
    sig28.to_excel(xw, sheet_name="Day28 significant DEGs", index=False)
    shared.to_excel(xw, sheet_name="Shared significant DEGs", index=False)
    flip.to_excel(xw, sheet_name="Direction-reversed genes", index=False)
    summary.to_excel(xw, sheet_name="Analysis summary", index=False)

# Auto-fit column widths
from openpyxl import load_workbook
wb = load_workbook(OUT)
for ws in wb.worksheets:
    for col in ws.columns:
        width = max((len(str(c.value)) for c in col if c.value is not None), default=10)
        ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 10), 40)
    ws.freeze_panes = "A2"
wb.save(OUT)

print("Output:", OUT)
print(f"All genes: {len(m)} rows")
print(f"Day14 significant: {n14} | Day28 significant: {n28} | shared: {len(shared)} (reversed {len(flip)})")
