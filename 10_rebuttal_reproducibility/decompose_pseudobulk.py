# -*- coding: utf-8 -*-
"""
伪bulk 分解 (GSE151974)

目的: 回答 Reviewer 2 的核心质疑 ——
  "CD3G/CD6 在 bulk 中下降, 到底是 T 细胞变少了, 还是每个 T 细胞表达变少了?"

做法:
  1. 按 Sample(6 个 MULTI-seq pool) 聚合 UMI, 得到伪bulk 的 CP10K 表达
     —— 这一步复现了"鼠肺 bulk RNA-seq"看到的现象
  2. 把伪bulk 表达拆成两部分:
       contribution(cell type) = 该类型细胞占比 x 该类型平均表达
     比较常氧 vs 高氧, 判断整体变化来自"比例"还是"表达强度"

输出:
  pseudobulk_by_sample.csv
  decomposition_contribution.csv
  分解结论 (打印 + 写入 decompose_report.md)
"""

import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

from scBPD_4gene_lung_mapping import (
    SIGNATURE, MARKERS, ALL_TARGETS,
    load_metadata, load_matrix_subset, normalize, log,
)

HERE = os.path.dirname(os.path.abspath(__file__))


def pseudobulk(mat, md):
    """按 Sample 聚合 UMI -> CP10K"""
    samp = md["Sample"].astype(str)
    rows = {}
    for s, idx in samp.groupby(samp).groups.items():
        idx = list(idx)
        tot = float(mat.loc[idx].values.sum())
        if tot <= 0:
            continue
        rows[s] = (mat.loc[idx].values.sum(axis=0) / tot * 1e4)
    df = pd.DataFrame.from_dict(rows, orient="index",
                                columns=mat.columns)
    df.index.name = "Sample"
    return df


def contribution(norm, md):
    """contribution = 细胞占比 x 该类型平均(未log)表达近似

    这里用 log 归一化的均值作为"表达强度"的可比量度, 并在报告中明确说明
    它是相对量而非绝对 UMI 拆解。
    """
    ct = md["CellType"].astype(str)
    ox = md["Oxygen"].astype(str)
    rows = []
    for ctype, idx in ct.groupby(ct).groups.items():
        idx = list(idx)
        n_idx = [i for i in idx if ox.get(i) == "Normoxia"]
        h_idx = [i for i in idx if ox.get(i) == "Hyperoxia"]
        for g in SIGNATURE:
            if g not in norm.columns:
                continue
            mn = norm[g].reindex(n_idx).dropna().mean() if n_idx else np.nan
            mh = norm[g].reindex(h_idx).dropna().mean() if h_idx else np.nan
            pn = 100.0 * len(n_idx) / int((ox == "Normoxia").sum())
            ph = 100.0 * len(h_idx) / int((ox == "Hyperoxia").sum())
            rows.append({
                "cell_type": ctype, "gene": g,
                "pct_normoxia": pn, "pct_hyperoxia": ph,
                "mean_expr_normoxia": mn, "mean_expr_hyperoxia": mh,
                "contrib_normoxia": pn * mn if mn == mn else np.nan,
                "contrib_hyperoxia": ph * mh if mh == mh else np.nan,
            })
    df = pd.DataFrame(rows)
    df["contrib_diff"] = df["contrib_hyperoxia"] - df["contrib_normoxia"]
    return df


def main():
    data = os.path.join(HERE, "data")
    md = load_metadata(data)
    mat = load_matrix_subset(data, ALL_TARGETS)
    norm = normalize(mat)

    pb = pseudulk = pseudobulk(mat, md)
    log("")
    log("=== 伪bulk CP10K (按 Sample) ===")
    cols = [c for c in SIGNATURE if c in pb.columns]
    log(pb[cols].round(2).to_string())

    ox = md["Oxygen"].astype(str)
    samp_ox = md.groupby(md["Sample"].astype(str))["Oxygen"].first().astype(str)
    log("")
    log("=== 高氧 vs 常氧 (伪bulk, 配对到 Age) ===")
    res = []
    for g in cols:
        by_age = {}
        for s in pb.index:
            a = s.split("_")[0]
            by_age.setdefault(a, {})[samp_ox.get(s, "")] = pb.loc[s, g]
        for a, d in sorted(by_age.items()):
            if "Normoxia" in d and "Hyperoxia" in d:
                n, h = d["Normoxia"], d["Hyperoxia"]
                res.append({"gene": g, "age": a, "normoxia": n, "hyperoxia": h,
                            "log2FC": np.log2((h + 1e-4) / (n + 1e-4))})
    rdf = pd.DataFrame(res)
    if not rdf.empty:
        log(rdf.round(3).to_string(index=False))
        log("")
        log("--- 每个基因的 log2FC 均值 ---")
        log(rdf.groupby("gene")["log2FC"].mean().round(3).to_string())

    contrib = contribution(norm, md)
    log("")
    log("=== 贡献度分解: 比例 vs 表达强度 (取每个基因贡献最大的 8 类细胞) ===")
    for g in cols:
        s = contrib[contrib["gene"] == g].sort_values("contrib_normoxia",
                                                      ascending=False)
        if s.empty:
            continue
        log(f"\n--- {g} ---")
        log(s.head(8)[["cell_type", "pct_normoxia", "pct_hyperoxia",
                       "mean_expr_normoxia", "mean_expr_hyperoxia",
                       "contrib_normoxia", "contrib_hyperoxia",
                       "contrib_diff"]].round(3).to_string(index=False))

    os.makedirs(HERE, exist_ok=True)
    pseudulk.to_csv(os.path.join(HERE, "pseudobulk_by_sample.csv"),
                    encoding="utf-8-sig")
    contrib.to_csv(os.path.join(HERE, "decomposition_contribution.csv"),
                   index=False, encoding="utf-8-sig")
    rdf.to_csv(os.path.join(HERE, "pseudobulk_log2FC.csv"),
               index=False, encoding="utf-8-sig")
    log("")
    log("[saved] pseudobulk_by_sample.csv, decomposition_contribution.csv, "
        "pseudobulk_log2FC.csv")


if __name__ == "__main__":
    main()
