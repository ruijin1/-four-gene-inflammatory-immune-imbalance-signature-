# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 01_data_prep/02_deg_limma_volcano_heatmap_venn.R
#  Status : RECONSTRUCTED from the Methods section - NOT validated in the original environment
#  Paths  : the authors' working copy used absolute paths on a local Windows
#           machine. They have been replaced here by P("<relative path>")
#           calls. Set BPD_REPO to the repository root (or run the script from
#           the repository root) and place inputs as described in docs/DATA.md.
#  Note   : the statistical code itself is unchanged.
# ==========================================================================
REPO_ROOT <- Sys.getenv("BPD_REPO", unset = "")
if (!nzchar(REPO_ROOT)) REPO_ROOT <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
P <- function(rel) file.path(REPO_ROOT, rel)

# =============================================================================
# 02  Differential expression analysis (limma) + volcano plot + heatmap + DEG/IRG Venn diagram
# Manuscript: Methods, "Differential expression analysis"; Fig 2A-C; Supplementary Data S4
# [RECONSTRUCTED] Rewritten from the written description in the Methods section, for code release (Zenodo/GitHub) and reviewer reproducibility.
# NOTE: not validated in the original authors' analysis environment. Verify locally against the existing result files before use.
# =============================================================================
library(limma); library(ggplot2); library(pheatmap); library(VennDiagram)

expr <- as.matrix(read.csv("data/GSE32472_day14_expr.csv", row.names = 1))
grp  <- read.csv("data/GSE32472_day14_group.csv")
group <- factor(grp$group, levels = c("control", "MS-BPD"))

## ---- limma differential expression: adj.P < 0.05 and |log2FC| > 0.58 (FC > 1.5) ----
design <- model.matrix(~ group)
fit    <- lmFit(expr, design)
fit    <- eBayes(fit)
deg_all <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
deg_all$symbol <- rownames(deg_all)
deg_sig <- subset(deg_all, adj.P.Val < 0.05 & abs(logFC) > 0.58)
write.csv(deg_all, "results/DEG_all.csv"); write.csv(deg_sig, "results/DEG_significant.csv")

## ---- Volcano plot (ggplot2; all genes and the significant DEGs labelled) ----
deg_all$sig <- ifelse(deg_all$adj.P.Val < 0.05 & deg_all$logFC >  0.58, "up",
               ifelse(deg_all$adj.P.Val < 0.05 & deg_all$logFC < -0.58, "down", "ns"))
p <- ggplot(deg_all, aes(logFC, -log10(adj.P.Val), colour = sig)) +
     geom_point(size = 0.8) + theme_bw() +
     scale_colour_manual(values = c(down = "blue", ns = "grey70", up = "red"))
ggsave("results/volcano_DEG.pdf", p, width = 5, height = 4)

## ---- Heatmap (Z-scores; Euclidean distance with complete linkage for rows and columns) ----
top <- head(deg_sig$symbol[order(deg_sig$adj.P.Val)], 50)
ph  <- pheatmap(expr[top, ], scale = "row", clustering_distance_rows = "euclidean",
                clustering_distance_cols = "euclidean", clustering_method = "complete",
                annotation_col = data.frame(Group = group, row.names = colnames(expr)))
ggsave("results/heatmap_DEG.pdf", ph$gtable, width = 7, height = 8)

## ---- DEG & IRG = DE-IRG (n = 57, Supplementary Data S4), drawn with VennDiagram ----
irg  <- read.csv("data/IRG_union_3336.csv")$symbol
deirg <- intersect(deg_sig$symbol, irg)
write.csv(data.frame(DE_IRG = deirg), "results/DE_IRG_57.csv", row.names = FALSE)
venn.diagram(list(DEG = deg_sig$symbol, IRG = irg), filename = "results/venn.tiff",
             imagetype = "tiff", height = 1800, width = 1800, resolution = 300)
