# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/rnaseq_deseq2.py
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
GSE220135 Generic RNA-seq differential-expression script (DESeq2 / pydeseq2)
Usage:
  python rnaseq_deseq2.py <input_csv> <output_xlsx> <timepoint_label>
Example:
  python rnaseq_deseq2.py GSE220135_28d.csv GSE220135_28d_DEG_results.xlsx Day28
Grouping rule: column names containing "_BPD" are cases, "_nonBPD" are controls (log2FC = BPD/nonBPD)
"""
import sys
import re
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

if len(sys.argv) < 4:
    print("Usage: python rnaseq_deseq2.py <input_csv> <output_xlsx> <timepoint_label>")
    sys.exit(1)

CSV_PATH = sys.argv[1]
OUT_XLSX = sys.argv[2]
TIME_LABEL = sys.argv[3]

# ---------- 1. Read data ----------
print("=" * 60)
print(f"STEP 1: Read data ({TIME_LABEL})")
df = pd.read_csv(CSV_PATH, index_col=0)
print(f"raw matrix: {df.shape[0]} genes x {df.shape[1]} samples")

# ---------- 2. Parse groups ----------
sample_cols = list(df.columns)
group = []
for c in sample_cols:
    if "_BPD" in c:
        group.append("BPD")
    elif "_nonBPD" in c:
        group.append("nonBPD")
    else:
        group.append("unknown")

df = df.loc[:, [g in ("BPD", "nonBPD") for g in group]]
sample_cols = list(df.columns)
group = [g for g in group if g in ("BPD", "nonBPD")]
df.columns = [f"{g}_{i}" for i, g in enumerate(group)]
metadata = pd.DataFrame({"condition": group}, index=df.columns)
n_bpd = group.count("BPD")
n_non = group.count("nonBPD")
print(f"sample groups: BPD={n_bpd}, nonBPD={n_non}, total={len(group)}")
print(f"comparison: BPD vs nonBPD (log2FC = BPD/nonBPD)")

# ---------- 3. Data cleaning ----------
df = df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
if df.index.duplicated().any():
    df = df.groupby(df.index).sum()
print(f"after filtering: {df.shape[0]} genes")

# ---------- 4. DESeq2 differential expression ----------
print("=" * 60)
print("STEP 2: DESeq2 differential expression (pydeseq2)")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

dds = DeseqDataSet(
    counts=df.T,
    metadata=metadata,
    design_factors="condition",
    refit_cooks=True,
    n_cpus=4,
)
dds.deseq2()

stat_res = DeseqStats(dds, contrast=["condition", "BPD", "nonBPD"], alpha=0.05)
stat_res.summary()

res = stat_res.results_df.copy()
res.index.name = "Gene"
print(f"genes tested: {res.shape[0]}")

# ---------- 5. Assemble results ----------
print("=" * 60)
print("STEP 3: Assemble results")
full = res[["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]].copy()
full.columns = ["baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj_FDR"]

sig = full.dropna(subset=["padj_FDR"])
sig = sig[(sig["padj_FDR"] < 0.05) & (sig["log2FC"].abs() > 1)].copy()
sig = sig.sort_values("padj_FDR")
up = sig[sig["log2FC"] > 0].sort_values("log2FC", ascending=False)
down = sig[sig["log2FC"] < 0].sort_values("log2FC")
print(f"significant DEGs (padj<0.05 and |log2FC|>1): {len(sig)}")
print(f"  up: {len(up)}   down: {len(down)}")

def annotate(g):
    return "ENSEMBL" if re.match(r"^ENSG\d+", str(g)) else "GeneSymbol"

full["gene_type"] = full.index.map(annotate)
sig["gene_type"] = sig.index.map(annotate)

# ---------- 6. Write Excel ----------
print("=" * 60)
print("STEP 4: Write Excel")
with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
    sig.to_excel(writer, sheet_name="Significant DEGs (up+down)")
    up.to_excel(writer, sheet_name="Up-regulated")
    down.to_excel(writer, sheet_name="Down-regulated")
    full.sort_values("padj_FDR", na_position="last").to_excel(writer, sheet_name="All genes")
    summary = pd.DataFrame({
        "Item": ["Dataset", "Time point", "Comparison", "n (BPD)", "n (nonBPD)",
                 "Test", "Significance threshold", "Total significant DEGs", "Up-regulated", "Down-regulated"],
        "Content": ["GSE220135", TIME_LABEL, "BPD vs nonBPD", n_bpd, n_non,
                 "DESeq2 (negative-binomial model, pydeseq2)", "padj<0.05 & |log2FC|>1",
                 len(sig), len(up), len(down)],
    })
    summary.to_excel(writer, sheet_name="Analysis summary", index=False)

print(f"saved: {OUT_XLSX}")
print("=" * 60)
print("Analysis complete!")
