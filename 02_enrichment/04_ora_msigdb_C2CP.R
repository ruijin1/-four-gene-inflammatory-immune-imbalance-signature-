# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 02_enrichment/04_ora_msigdb_C2CP.R
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
# 04  ORA (clusterProfiler::enricher; MSigDB C2.CP v2025.1; background = genes detected in the training set)
# Manuscript: Methods, "Over-representation analysis (ORA)"; supplementary to Fig 3
# [RECONSTRUCTED] Rewritten from the written description in the Methods section, for code release (Zenodo/GitHub) and reviewer reproducibility.
# NOTE: not validated in the original authors' analysis environment. Verify locally against the existing result files before use.
# =============================================================================
library(clusterProfiler); library(msigdbr); library(dplyr)

irg_1498 <- read.csv("data/IRG_detected_1498.csv")$symbol      # Supplementary Data S5
expr     <- as.matrix(read.csv("data/GSE32472_day14_expr.csv", row.names = 1))
bg       <- rownames(expr)

c2cp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP") %>%
        dplyr::select(gs_name, gene_symbol)
ora  <- enricher(irg_1498, TERM2GENE = c2cp, universe = bg,
                 pAdjustMethod = "BH", pvalueCutoff = 0.05)     # Fisher's exact test + BH adjustment
write.csv(as.data.frame(ora), "results/ORA_C2CP.csv")
