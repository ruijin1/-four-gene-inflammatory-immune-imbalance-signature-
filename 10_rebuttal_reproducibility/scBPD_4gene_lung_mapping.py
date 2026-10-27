# -*- coding: utf-8 -*-
"""
GSE151974 (Hurskainen 2021, Nat Commun) 新生鼠高氧肺单细胞
        -> 人血四基因签名 (CD3G / PRTN3 / CD6 / THBS1) 的肺内细胞定位

直接回应 Reviewer 2 的三条 major:
  R2#1  鼠肺 bulk 只能看方向, 不能定细胞来源
        -> 这里用真实高氧肺 scRNA-seq (61,839 cells, 34 种注释细胞类型)
           给出四基因在哪类细胞表达
  R2#2  A549 不能验证免疫签名
        -> 这里给出体内免疫细胞的真实组成变化 (直接计数, 非去卷积推断)
  R2#3  "中性粒细胞-T 细胞免疫失衡" 证据不足 / CIBERSORT 只是推断
        -> 这里用细胞类型标签直接计数 + Fisher 检验 + BH 校正

数据: GSE151974
  - GSE151974_cell_metadata_postfilter.csv.gz  (61,839 cells, 含 CellType/Oxygen/Age)
  - GSE151974_raw_umi_matrix_postfilter.csv    (22,100 genes x 61,839 cells, 解压后 2.55 GB)

用法:
    python scBPD_4gene_lung_mapping.py --data <目录>
"""

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------- 基因集 ----
# 人 -> 小鼠同源符号 (小鼠用首字母大写其余小写)
SIGNATURE = ["Cd3g", "Prtn3", "Cd6", "Thbs1"]

MARKERS = {
    "T cell":      ["Cd3d", "Cd3e", "Cd3g", "Cd2", "Lck", "Cd4", "Cd8a", "Cd8b1", "Cd6", "Thy1"],
    "Neutrophil":  ["Prtn3", "Elane", "Ctsg", "Mpo", "Itgam", "Csf3r", "Camp", "Ltf",
                    "Ngp", "Lcn2", "Mmp8", "Mmp9", "S100a8", "S100a9", "Fcnb"],
    "Macrophage":  ["Cd68", "Adgre1", "Lyz2", "Csf1r", "Cd14", "Fcgr3", "Mrc1", "Siglecf"],
    "Monocyte":    ["Lyz2", "Cd14", "Fcgr1", "Ccr2", "Ly6c2", "Cd68"],
    "B cell":      ["Cd19", "Ms4a1", "Cd79a", "Cd79b", "H2-Ob"],
    "NK/ILC":      ["Nkg7", "Klrb1c", "Gzma", "Ncr1", "Eomes", "Il7r", "Rora"],
    "DC":          ["Cd209a", "Flt3", "Ccr7", "H2-Ab1", "Itgax"],
    "Platelet":    ["Pf4", "Ppbp", "Itga2b", "Tubb1", "Gp9"],
    "Endothelial": ["Pecam1", "Cdh5", "Cldn5", "Ager", "Cav1"],
    "Epithelial":  ["Sftpc", "Sftpb", "Sftpa1", "Scgb1a1", "Foxj1", "Aqp5", "Epcam", "Cdh1"],
    "Fibroblast":  ["Col1a1", "Col13a1", "Col14a1", "Pdgfra", "Acta2"],
}

ALL_TARGETS = sorted(set(SIGNATURE) |
                     {g for v in MARKERS.values() for g in v})

# 免疫细胞类型关键词(用于组成分析的分组着色)
IMMUNE_KW = ["T cell", "Neut", "Mf", "Mono", "B cell", "NK", "DC", "ILC",
             "Mast", "gd T"]


def log(msg=""):
    print(msg, flush=True)


def load_metadata(data_dir):
    p = os.path.join(data_dir, "GSE151974_cell_metadata_postfilter.csv.gz")
    if not os.path.exists(p):
        sys.exit(f"[FATAL] metadata not found: {p}")
    md = pd.read_csv(p, index_col=0)
    md.columns = [c.strip() for c in md.columns]
    need = ["CellType", "Oxygen", "Age", "Sample"]
    miss = [c for c in need if c not in md.columns]
    if miss:
        sys.exit(f"[FATAL] metadata missing columns: {miss}")
    return md


def load_matrix_subset(data_dir, targets):
    """流式扫描 2.55 GB 的 csv, 只保留目标基因的行。"""
    p = os.path.join(data_dir, "GSE151974_raw_umi_matrix_postfilter.csv")
    if not os.path.exists(p):
        sys.exit(f"[FATAL] matrix not found: {p}\n"
                 f"        (注意: 若文件名带 .csv 却是目录, 真实文件在其内部)")
    tset = {t.lower(): t for t in targets}
    found = {}
    with open(p, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split(",")
        cellnames = header[1:]
        for line in f:
            g = line.split(",", 1)[0].strip().strip('"')
            k = g.lower()
            if k in tset:
                vals = np.asarray(line.rstrip("\n").split(",")[1:],
                                  dtype=np.float32)
                found[g] = vals
    log(f"[matrix] header cells = {len(cellnames)}")
    log(f"[matrix] genes scanned = 22100, targets kept = {len(found)}/{len(targets)}")
    missing = [t for t in targets if t not in found]
    if missing:
        log(f"[matrix] WARNING absent: {missing}")
    return pd.DataFrame(found, index=cellnames)


def normalize(mat):
    """CP10K + log1p (标准的 scRNA-seq 归一化)"""
    # pandas 3.x 下 .values 可能是只读视图, 先 np.array(...) 拷贝一份
    tot = np.array(mat.sum(axis=1).values, dtype=np.float64).reshape(-1, 1)
    tot[tot == 0] = 1.0
    norm = np.log1p(np.array(mat.values, dtype=np.float64) / tot * 1e4)
    return pd.DataFrame(norm, index=mat.index, columns=mat.columns)


def bh(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    if n == 0:
        return np.array([])
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1))
        q[o[i]] = prev
    return np.minimum(q, 1.0)


def analyse_localisation(norm, md):
    """四基因在各注释细胞类型中的平均表达 + 表达细胞比例"""
    ct = md["CellType"].astype(str)
    out = []
    for g in SIGNATURE:
        if g not in norm.columns:
            continue
        v = norm[g]
        for ctype, idx in ct.groupby(ct).groups.items():
            x = v.reindex(idx).dropna()
            out.append({
                "gene": g,
                "cell_type": ctype,
                "n_cells": len(x),
                "mean_log_expr": float(x.mean()),
                "pct_expressing": float(100 * (x > 0).mean()),
            })
    df = pd.DataFrame(out)
    if df.empty:
        return df
    df = df.sort_values(["gene", "mean_log_expr"], ascending=[True, False])
    return df


def analyse_composition(md):
    """高氧 vs 常氧 各细胞类型比例变化 + Fisher 检验 + BH 校正"""
    ct = md["CellType"].astype(str)
    ox = md["Oxygen"].astype(str)
    tot = {"Normoxia": int((ox == "Normoxia").sum()),
           "Hyperoxia": int((ox == "Hyperoxia").sum())}
    tab = pd.crosstab(ct, ox)
    for c in ("Normoxia", "Hyperoxia"):
        if c not in tab.columns:
            tab[c] = 0

    rows = []
    for ctype in tab.index:
        n = int(tab.loc[ctype, "Normoxia"])
        h = int(tab.loc[ctype, "Hyperoxia"])
        # 2x2: [this type norm, this type hyper] vs [other norm, other hyper]
        a, b = n, h
        c = tot["Normoxia"] - n
        d = tot["Hyperoxia"] - h
        odds, p = stats.fisher_exact([[a, b], [c, d]])
        pn = 100.0 * n / tot["Normoxia"] if tot["Normoxia"] else np.nan
        ph = 100.0 * h / tot["Hyperoxia"] if tot["Hyperoxia"] else np.nan
        rows.append({
            "cell_type": ctype,
            "n_normoxia": n,
            "n_hyperoxia": h,
            "pct_normoxia": pn,
            "pct_hyperoxia": ph,
            "pct_diff": ph - pn,
            "log2FC": np.log2((ph / 100 + 1e-4) / (pn / 100 + 1e-4)),
            "odds_ratio": odds,
            "p_value": p,
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["FDR"] = bh(df["p_value"].values)
        df = df.sort_values("p_value")
    return df


def analyse_composition_by_age(md):
    """分 P3 / P7 / P14 三个时间点看免疫组成变化"""
    ct = md["CellType"].astype(str)
    age = md["Age"].astype(str)
    ox = md["Oxygen"].astype(str)
    rows = []
    for a in ["P3", "P7", "P14"]:
        sel = (age == a)
        sub = md[sel]
        if sub.empty:
            continue
        t_n = int((sub["Oxygen"].astype(str) == "Normoxia").sum())
        t_h = int((sub["Oxygen"].astype(str) == "Hyperoxia").sum())
        cts = ct[sel]
        for ctype in cts.value_counts().index:
            n = int(((cts == ctype) & (sub["Oxygen"].astype(str) == "Normoxia")).sum())
            h = int(((cts == ctype) & (sub["Oxygen"].astype(str) == "Hyperoxia")).sum())
            pn = 100.0 * n / t_n if t_n else np.nan
            ph = 100.0 * h / t_h if t_h else np.nan
            rows.append({"age": a, "cell_type": ctype,
                         "n_normoxia": n, "n_hyperoxia": h,
                         "pct_normoxia": pn, "pct_hyperoxia": ph,
                         "pct_diff": ph - pn})
    return pd.DataFrame(rows)


def analyse_de_within_celltypes(norm, md, min_cells=30):
    """在每种细胞类型内部, 比较高氧 vs 常氧的四基因表达 (Wilcoxon)"""
    ct = md["CellType"].astype(str)
    ox = md["Oxygen"].astype(str)
    rows = []
    for ctype, idx in ct.groupby(ct).groups.items():
        i_n = [i for i in idx if ox.get(i) == "Normoxia"]
        i_h = [i for i in idx if ox.get(i) == "Hyperoxia"]
        if len(i_n) < min_cells or len(i_h) < min_cells:
            continue
        for g in SIGNATURE:
            if g not in norm.columns:
                continue
            a = norm[g].reindex(i_n).dropna().values
            b = norm[g].reindex(i_h).dropna().values
            if len(a) < min_cells or len(b) < min_cells:
                continue
            try:
                u, p = stats.mannwhitneyu(b, a, alternative="two-sided")
            except ValueError:
                continue
            rows.append({
                "cell_type": ctype,
                "gene": g,
                "n_normoxia": len(a),
                "n_hyperoxia": len(b),
                "mean_normoxia": float(a.mean()),
                "mean_hyperoxia": float(b.mean()),
                "log2FC_hyper_vs_norm": float(np.log2((b.mean() + 1e-4) /
                                                      (a.mean() + 1e-4))),
                "p_value": float(p),
            })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["FDR"] = bh(df["p_value"].values)
    return df


# ------------------------------------------------------------------ 绘图 ----
def make_figures(loc, comp, comp_age, out_dir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception as e:
        log(f"[plot] skipped: {e}")
        return []

    plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    paths = []

    # Fig A: 四基因在各细胞类型的平均表达
    if not loc.empty:
        piv = loc.pivot(index="cell_type", columns="gene",
                        values="mean_log_expr")
        piv = piv.reindex(comp["cell_type"].tolist()) if not comp.empty else piv
        piv = piv.dropna(how="all")
        fig, ax = plt.subplots(figsize=(7.2, max(4, 0.28 * len(piv) + 1.6)))
        sns.heatmap(piv, cmap="viridis", ax=ax, cbar_kws={"label": "mean log-normalised expression"})
        ax.set_title("Four-gene signature across annotated lung cell types (GSE151974)")
        ax.set_xlabel(""); ax.set_ylabel("")
        plt.tight_layout()
        p = os.path.join(out_dir, "figA_4gene_localisation_heatmap.png")
        fig.savefig(p, dpi=300); plt.close(fig); paths.append(p)

    # Fig B: 细胞组成变化 (比例差)
    if not comp.empty:
        d = comp.sort_values("pct_diff")
        fig, ax = plt.subplots(figsize=(7.2, max(4, 0.28 * len(d) + 1.6)))
        colors = ["#c0392b" if v > 0 else "#2980b9" for v in d["pct_diff"]]
        ax.barh(d["cell_type"], d["pct_diff"], color=colors)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel("Hyperoxia % - Normoxia % (percentage points)")
        ax.set_title("Shift in lung cell composition under hyperoxia")
        plt.tight_layout()
        p = os.path.join(out_dir, "figB_composition_shift.png")
        fig.savefig(p, dpi=300); plt.close(fig); paths.append(p)

    # Fig C: 关键免疫类群的组成对比
    if not comp.empty:
        keys = ["Neut 1", "Neut 2", "CD4 T cell 1", "CD4 T cell 2",
                "CD8 T cell 1", "CD8 T cell 2", "gd T cell", "B cell 1",
                "NK cell", "Alv Mf", "Int Mf", "Mono"]
        sub = comp[comp["cell_type"].isin(keys)].copy()
        if not sub.empty:
            sub["__o"] = sub["cell_type"].apply(lambda x: keys.index(x))
            sub = sub.sort_values("__o")
            x = np.arange(len(sub)); w = 0.38
            fig, ax = plt.subplots(figsize=(8, 4.2))
            ax.bar(x - w/2, sub["pct_normoxia"], w, label="Normoxia", color="#7fb3d5")
            ax.bar(x + w/2, sub["pct_hyperoxia"], w, label="Hyperoxia", color="#e07b58")
            ax.set_xticks(x); ax.set_xticklabels(sub["cell_type"], rotation=35, ha="right")
            ax.set_ylabel("% of all cells in condition")
            ax.set_title("Immune cell composition: normoxia vs hyperoxia")
            ax.legend()
            plt.tight_layout()
            p = os.path.join(out_dir, "figC_immune_composition.png")
            fig.savefig(p, dpi=300); plt.close(fig); paths.append(p)

    # Fig D: 分时间点的中性粒 / T 细胞比例
    if not comp_age.empty:
        keys = ["Neut 1", "Neut 2", "CD4 T cell 1", "CD8 T cell 1"]
        sub = comp_age[comp_age["cell_type"].isin(keys)]
        if not sub.empty:
            fig, ax = plt.subplots(figsize=(6.6, 4))
            sns.barplot(data=sub, x="age", y="pct_hyperoxia", hue="cell_type", ax=ax)
            ax.set_ylabel("% of cells (hyperoxia)")
            ax.set_title("Neutrophil / T-cell fractions across postnatal days")
            ax.legend(fontsize=8)
            plt.tight_layout()
            p = os.path.join(out_dir, "figD_age_stratified.png")
            fig.savefig(p, dpi=300); plt.close(fig); paths.append(p)

    return paths


# ------------------------------------------------------------------ 报告 ----
def write_report(md, loc, comp, comp_age, de, figs, out_dir):
    L = []
    A = L.append
    A("# 四基因签名在新生鼠高氧肺单细胞中的定位")
    A("")
    A(f"**数据集**：GSE151974（Hurskainen M et al., *Nat Commun* 2021）  ")
    A("**设计**：新生鼠常氧 vs 高氧，P3 / P7 / P14 三个时间点，MULTI-seq 混样  ")
    A(f"**细胞数**：{len(md):,}  ")
    A(f"**注释细胞类型**：{md['CellType'].nunique()} 种")
    A("")
    A("## 1. 样本构成")
    A("")
    A("| 分组 | 细胞数 |")
    A("|---|---|")
    for k, v in md["Sample"].astype(str).value_counts().sort_index().items():
        A(f"| {k} | {v:,} |")
    A("")

    A("## 2. 四基因的细胞定位")
    A("")
    if not loc.empty:
        for g in SIGNATURE:
            s = loc[loc["gene"] == g].sort_values("mean_log_expr",
                                                  ascending=False)
            if s.empty:
                continue
            A(f"**{g}** — 表达最高的 6 类细胞：")
            A("")
            A("| 细胞类型 | 细胞数 | 平均表达(log) | 表达细胞比例(%) |")
            A("|---|---|---|---|")
            for _, r in s.head(6).iterrows():
                A(f"| {r['cell_type']} | {int(r['n_cells']):,} | "
                  f"{r['mean_log_expr']:.3f} | {r['pct_expressing']:.1f} |")
            A("")

    A("## 3. 高氧对肺内细胞组成的影响")
    A("")
    if not comp.empty:
        A("| 细胞类型 | 常氧 % | 高氧 % | 差值(百分点) | OR | P | FDR |")
        A("|---|---|---|---|---|---|---|")
        for _, r in comp.head(20).iterrows():
            A(f"| {r['cell_type']} | {r['pct_normoxia']:.2f} | "
              f"{r['pct_hyperoxia']:.2f} | {r['pct_diff']:+.2f} | "
              f"{r['odds_ratio']:.2f} | {r['p_value']:.2e} | {r['FDR']:.2e} |")
        A("")

    A("## 4. 四基因在细胞类型内部的表达差异（高氧 vs 常氧）")
    A("")
    if not de.empty:
        sig = de[(de["FDR"] < 0.05)]
        if not sig.empty:
            A(f"FDR < 0.05 的组合共 **{len(sig)}** 个：")
            A("")
            A("| 细胞类型 | 基因 | 常氧 | 高氧 | log2FC | FDR |")
            A("|---|---|---|---|---|---|")
            for _, r in sig.sort_values("FDR").head(25).iterrows():
                A(f"| {r['cell_type']} | {r['gene']} | {r['mean_normoxia']:.3f} | "
                  f"{r['mean_hyperoxia']:.3f} | {r['log2FC_hyper_vs_norm']:+.3f} | "
                  f"{r['FDR']:.2e} |")
            A("")
        else:
            A("无 FDR < 0.05 的组合。")
            A("")

    if figs:
        A("## 5. 图件")
        A("")
        for f in figs:
            A(f"![{os.path.basename(f)}]({os.path.basename(f)})")
        A("")

    A("## 6. 使用限制（务必在稿件/rebuttal 中如实说明）")
    A("")
    A("- 本数据为**小鼠肺组织**，非人外周血；仅用于说明四个基因在肺内的**细胞表达谱**，"
      "不能直接外推为人血中该基因的细胞来源。")
    A("- 高氧模型为新生鼠早期（P3–P14）暴露，与人 BPD 的慢性演进并不完全等价。")
    A("- 细胞类型标签来自原作者注释（deMULTIplex 拆分 + 人工注释），本文未重新分群。")
    A("- 组成变化的比较以**比例**而非绝对细胞数为依据，因两组总细胞数不同"
      f"（常氧 {int((md['Oxygen'].astype(str)=='Normoxia').sum()):,} vs "
      f"高氧 {int((md['Oxygen'].astype(str)=='Hyperoxia').sum()):,}）。")

    p = os.path.join(out_dir, "四基因肺单细胞定位报告.md")
    open(p, "w", encoding="utf-8").write("\n".join(L))
    return p


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--data", default=os.path.join(here, "data"))
    ap.add_argument("--out", default=here)
    a = ap.parse_args()

    t0 = time.time()
    log(f"[data] {a.data}")

    md = load_metadata(a.data)
    log(f"[meta] cells = {len(md):,}, cell types = {md['CellType'].nunique()}")
    log(f"[meta] oxygen: {dict(md['Oxygen'].astype(str).value_counts())}")
    log(f"[meta] age:    {dict(md['Age'].astype(str).value_counts())}")

    mat = load_matrix_subset(a.data, ALL_TARGETS)
    log(f"[mat ] subset = {mat.shape}")

    norm = normalize(mat)
    log(f"[norm] done ({time.time()-t0:.0f}s)")

    loc = analyse_localisation(norm, md)
    comp = analyse_composition(md)
    comp_age = analyse_composition_by_age(md)
    de = analyse_de_within_celltypes(norm, md)

    os.makedirs(a.out, exist_ok=True)
    loc.to_csv(os.path.join(a.out, "4gene_localisation.csv"),
               index=False, encoding="utf-8-sig")
    comp.to_csv(os.path.join(a.out, "composition_shift.csv"),
                index=False, encoding="utf-8-sig")
    comp_age.to_csv(os.path.join(a.out, "composition_by_age.csv"),
                    index=False, encoding="utf-8-sig")
    de.to_csv(os.path.join(a.out, "de_within_celltypes.csv"),
              index=False, encoding="utf-8-sig")

    figs = make_figures(loc, comp, comp_age, a.out)
    rp = write_report(md, loc, comp, comp_age, de, figs, a.out)

    log("")
    log(f"[done] {time.time()-t0:.0f}s")
    log(f"[report] {rp}")
    log("[figs] " + ", ".join(os.path.basename(f) for f in figs))


if __name__ == "__main__":
    main()
