# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 02_enrichment/06_GSVA_hallmark.R
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
# 06  GSVA (GSVA 1.44.3 + Hallmark) + limma comparison between high- and low-expression groups
# Manuscript: Methods, "Gene set variation analysis"; Fig 8
# criterion: |t| > 2 and BH-adjusted p < 0.05
# [RECONSTRUCTED] Rewritten from the written description in the Methods section, for code release (Zenodo/GitHub) and reviewer reproducibility.
# NOTE: not validated in the original authors' analysis environment. Verify locally against the existing result files before use.
# =============================================================================
library(GSVA); library(msigdbr); library(limma)

expr  <- as.matrix(read.csv("data/GSE32472_day14_expr.csv", row.names = 1))
hall  <- msigdbr(species = "Homo sapiens", category = "H")
sets  <- split(hall$gene_symbol, hall$gs_name)
sets  <- lapply(sets, function(x) intersect(x, rownames(expr)))

gsva_par <- gsvaParam(expr, sets)             # alternatively gsva(expr, sets, method = "gsva")
scores   <- gsva(gsva_par)

for (g in c("CD3G", "PRTN3", "CD6", "THBS1")) {
  hi <- expr[g, ] >= median(expr[g, ])
  group <- factor(ifelse(hi, "high", "low"), levels = c("low", "high"))
  fit <- eBayes(lmFit(scores, model.matrix(~ group)))
  tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "t")
  tt$pathway <- rownames(tt)
  tt$sig <- abs(tt$t) > 2 & tt$adj.P.Val < 0.05   # positive t = enriched in the high-expression group
  write.csv(tt, paste0("results/GSVA_", g, ".csv"))
}
