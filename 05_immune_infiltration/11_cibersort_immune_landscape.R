# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 05_immune_infiltration/11_cibersort_immune_landscape.R
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
# 11  CIBERSORT (1.06): 22 immune cell fractions + between-group comparison + correlation heatmaps
# Manuscript: Methods, "Immune infiltration analysis" / "Statistical analysis"; Fig 6A-D
# [RECONSTRUCTED] Rewritten from the written description in the Methods section, for code release (Zenodo/GitHub) and reviewer reproducibility.
# NOTE: not validated in the original authors' analysis environment. Verify locally against the existing result files before use.
# =============================================================================
library(CIBERSORT); library(pheatmap); library(ggplot2)

expr <- as.matrix(read.csv("data/GSE32472_day14_expr.csv", row.names = 1))
grp  <- read.csv("data/GSE32472_day14_group.csv")
lm22 <- read.table("data/LM22.txt", header = TRUE, sep = "\t", row.names = 1)   # CIBERSORT official signature matrix

res <- CIBERSORT(lm22, expr, perm = 1000, QN = FALSE)
res <- as.data.frame(res)
res <- res[res$`P-value` < 0.05, ]              # keep only samples that passed deconvolution
write.csv(res, "results/CIBERSORT_22cells.csv")

## Between-group comparison: Mann-Whitney U test (BH adjustment)
cells <- colnames(res)[1:22]
pvals <- sapply(cells, function(c) wilcox.test(res[[c]] ~ grp$group)$p.value)
out <- data.frame(cell = cells, p = pvals, padj = p.adjust(pvals, "BH"))
write.csv(out, "results/CIBERSORT_group_diff.csv", row.names = FALSE)

## Correlation heatmaps (Spearman): cell-cell and cell-hub gene
hub <- c("CD3G", "PRTN3", "CD6", "THBS1")
cm  <- cor(res[, cells], t(expr[hub, rownames(res)]), method = "spearman")
pheatmap(cm, cluster_rows = TRUE, cluster_cols = TRUE,
         filename = "results/heatmap_immune_hub.pdf", width = 6, height = 8)
