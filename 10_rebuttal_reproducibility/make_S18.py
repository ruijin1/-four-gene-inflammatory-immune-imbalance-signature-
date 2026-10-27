# -*- coding: utf-8 -*-
"""Build Supplementary Data S18: marker-set analysis (R1 major 4) + GSE151974 tables."""
import os

import pandas as pd

SRC = r"C:\Users\44937\Desktop\投稿内容\送审包_20260911\答复_审稿意见_20261026\02_补充分析"
SI = r"C:\Users\44937\Desktop\投稿内容\送审包_20260911\Supplementary_FrontPediatr"
OUT = os.path.join(SI, "Data S18.xlsx")

sheets = [
    ("Marker_set_analysis", "marker_check_results.csv"),
    ("GSE151974_composition", "composition_shift.csv"),
    ("GSE151974_within_celltype", "de_within_celltypes.csv"),
    ("GSE151974_decomposition", "decomposition_contribution.csv"),
    ("GSE151974_pseudobulk", "pseudobulk_log2FC.csv"),
]

with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
    for name, fn in sheets:
        p = os.path.join(SRC, fn)
        df = pd.read_csv(p, encoding="utf-8-sig")
        df.to_excel(xw, sheet_name=name[:31], index=False)
        print("  %-28s %-32s %d rows x %d cols" % (name, fn, len(df), df.shape[1]))

print()
print("written:", OUT, os.path.getsize(OUT), "bytes")

# read back
xl = pd.ExcelFile(OUT)
print("sheets:", xl.sheet_names)
