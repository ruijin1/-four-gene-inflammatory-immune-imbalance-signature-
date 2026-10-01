# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 02_enrichment/05_single_gene_GSEA_KEGG.R
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
# 05  Single-gene GSEA (median split -> limma -> ranking by log2FC -> clusterProfiler pre-ranked GSEA)
# Manuscript: Methods, "Single-gene gene set enrichment analysis"; Fig 7
# criterion: |NES| > 1.5 and FDR < 0.25
# [RECONSTRUCTED] Rewritten from the written description in the Methods section, for code release (Zenodo/GitHub) and reviewer reproducibility.
# NOTE: not validated in the original authors' analysis environment. Verify locally against the existing result files before use.
# =============================================================================
library(limma); library(clusterProfiler); library(org.Hs.eg.db)

expr <- as.matrix(read.csv("data/GSE32472_day14_expr.csv", row.names = 1))
genes <- c("CD3G", "PRTN3", "CD6", "THBS1")

for (g in genes) {
  hi <- expr[g, ] >= median(expr[g, ])          # high-expression group n = 32 / low-expression group n = 31
  group <- factor(ifelse(hi, "high", "low"), levels = c("low", "high"))
  fit <- eBayes(lmFit(expr, model.matrix(~ group)))
  tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "logFC")
  ranks <- setNames(tt$logFC, rownames(tt))
  ranks <- sort(ranks, decreasing = TRUE)

  eg <- bitr(names(ranks), "SYMBOL", "ENTREZID", org.Hs.eg.db)
  ranks <- ranks[eg$SYMBOL]; names(ranks) <- eg$ENTREZID

  gsea <- gseKEGG(ranks, organism = "hsa", minGSSize = 10, maxGSSize = 500,
                  pvalueCutoff = 1, verbose = FALSE)
  out <- as.data.frame(gsea)
  sig <- subset(out, abs(NES) > 1.5 & p.adjust < 0.25)
  write.csv(out, paste0("results/GSEA_KEGG_", g, ".csv"))
  write.csv(sig, paste0("results/GSEA_KEGG_", g, "_significant.csv"))
}
