# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 02_enrichment/03_go_kegg_circle_plot.R
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
# 03  GO / KEGG enrichment (clusterProfiler 4.4.4) + logFC circle plot
# Manuscript: Methods, "GO and KEGG analyses"; Fig 3
# [RECONSTRUCTED] Rewritten from the written description in the Methods section, for code release (Zenodo/GitHub) and reviewer reproducibility.
# NOTE: not validated in the original authors' analysis environment. Verify locally against the existing result files before use.
# =============================================================================
library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2); library(dplyr)

deirg  <- read.csv("results/DE_IRG_57.csv")$DE_IRG
deg    <- read.csv("results/DEG_significant.csv")
expr   <- as.matrix(read.csv("data/GSE32472_day14_expr.csv", row.names = 1))
bg     <- rownames(expr)                      # genes detected on the platform, used as background
eg     <- bitr(deirg, "SYMBOL", "ENTREZID", org.Hs.eg.db)
bg_eg  <- bitr(bg,    "SYMBOL", "ENTREZID", org.Hs.eg.db)

go <- enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, universe = bg_eg$ENTREZID,
               ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
kk <- enrichKEGG(eg$ENTREZID, organism = "hsa", universe = bg_eg$ENTREZID,
                 pAdjustMethod = "BH", pvalueCutoff = 0.05)
write.csv(as.data.frame(go), "results/GO_enrichment.csv")
write.csv(as.data.frame(kk), "results/KEGG_enrichment.csv")

## Circle plot: each DE-IRG mapped to its enriched terms, coloured by logFC (red > 0 / blue < 0)
circ <- data.frame(term = rep(go$Description, lengths(strsplit(go$geneID, "/"))),
                   gene = unlist(strsplit(go$geneID, "/"))) %>%
        left_join(deg[, c("symbol", "logFC")], by = c("gene" = "symbol"))
p <- ggplot(circ, aes(x = term, y = gene, fill = logFC)) +
     geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
     theme_minimal() + coord_polar()
ggsave("results/GO_KEGG_circle.pdf", p, width = 8, height = 8)
