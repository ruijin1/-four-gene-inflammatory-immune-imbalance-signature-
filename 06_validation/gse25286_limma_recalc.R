# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 06_validation/gse25286_limma_recalc.R
#  Status : ORIGINAL (as used in the analysis)
#  Paths  : the authors' working copy used absolute paths on a local Windows
#           machine. They have been replaced here by P("<relative path>")
#           calls. Set BPD_REPO to the repository root (or run the script from
#           the repository root) and place inputs as described in docs/DATA.md.
#  Note   : the statistical code itself is unchanged.
# ==========================================================================
REPO_ROOT <- Sys.getenv("BPD_REPO", unset = "")
if (!nzchar(REPO_ROOT)) REPO_ROOT <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
P <- function(rel) file.path(REPO_ROOT, rel)

# ----------------------------------------------------------------------------
# GSE25286 (hyperoxia-induced BPD mouse model, lung tissue, GPL1261)
# Cross-species validation of the four hub genes (THBS1, CD6, CD3G, PRTN3)
# Reproduction of the manuscript "S9 Data" table from the raw expression matrix
# using limma (exact pooled P / adjusted P + time-point logFC).
#
# Methods mirrored from the manuscript:
#   - pooled analysis  : all P14 + P29 samples, hyperoxia ("test") vs normoxia ("ref")
#   - moderation       : limma empirical Bayes (eBayes) across the full array
#   - multiple testing : Benjamini-Hochberg over ALL probes on the array
#   - time-course      : P14 and P29 subgroups, logFC only (no tests, n_ref=2)
#
# Usage:
#   1) Install dependencies once:
#        install.packages("BiocManager")
#        BiocManager::install("limma")
#        install.packages("readxl")
#   2) Run:  Rscript gse25286_limma_recalc.R
#      (R >= 4.2 recommended; run inside RStudio if the console shows garbled
#       characters when printing the hard-coded data path.)
#   3) The console prints computed vs. manuscript-S9 values side by side;
#      CSV copies are written to OUT_DIR.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(limma)
  library(readxl)
})

# ----------------------------- paths ---------------------------------------
DATA_DIR <- P("data/raw/GSE32472/GSE25286_mouse_validation/complete")
EXPR_FILE <- file.path(DATA_DIR, "expression_matrix.csv")
META_FILE <- file.path(DATA_DIR, "sample_metadata.xlsx")
OUT_DIR   <- P("results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- read inputs ---------------------------------
if (!file.exists(EXPR_FILE)) {
  message("Expression file not found at:\n  ", EXPR_FILE)
  DATA_DIR <- readline("Paste the full path of the 'complete' data folder and press Enter:\n")
  EXPR_FILE <- file.path(DATA_DIR, "expression_matrix.csv")
  META_FILE <- file.path(DATA_DIR, "sample_metadata.xlsx")
}
expr_all <- read.csv(EXPR_FILE, check.names = FALSE, fileEncoding = "UTF-8-BOM",
                     stringsAsFactors = FALSE)
stopifnot(all(c("id", "anno") %in% colnames(expr_all)))
gsm_cols <- setdiff(colnames(expr_all), c("id", "anno"))

meta <- tryCatch(as.data.frame(read_excel(META_FILE)), error = function(e) NULL)
if (!is.null(meta)) {
  meta$time  <- sub("^(P14|P29).*", "\\1", meta$title)
  meta$group <- tolower(as.character(meta$group))     # "ref" (normoxia) / "test" (hyperoxia)
  rownames(meta) <- as.character(meta$gsm)
  meta <- meta[gsm_cols, ]                            # align to matrix column order
} else {
  meta <- data.frame(gsm   = gsm_cols,
                     time  = c("P14","P14","P29","P29","P14","P14","P14","P29","P29","P29"),
                     group = c("ref","ref","ref","ref","test","test","test","test","test","test"),
                     stringsAsFactors = FALSE)
}
stopifnot(all(meta$group %in% c("ref", "test")),
          all(c("time", "group") %in% colnames(meta)))

X <- as.matrix(expr_all[gsm_cols])
rownames(X) <- expr_all$id
message(sprintf("Loaded %d probes x %d samples (%d ref, %d test)",
                nrow(X), ncol(X), sum(meta$group == "ref"), sum(meta$group == "test")))

# ------------------------- target probes / genes ----------------------------
# THBS1: manuscript S9 values correspond to probe 1460302_at (both probes reported).
probe2gene <- c("1460302_at" = "THBS1 (S9 primary)",
                "1450377_at" = "THBS1 (2nd probe)",
                "1451910_a_at" = "CD6",
                "1419178_at"  = "CD3G",
                "1419669_at"  = "PRTN3")
sel <- names(probe2gene)[names(probe2gene) %in% rownames(X)]

# ------------------------- 1) pooled limma analysis --------------------------
grp    <- factor(meta$group, levels = c("ref", "test"))
design <- model.matrix(~ grp)
fit    <- eBayes(lmFit(X, design))
tt     <- topTable(fit, coef = "grptest", number = Inf, sort.by = "none")

pooled <- data.frame(gene  = unname(probe2gene[sel]),
                     probe = sel,
                     logFC = round(tt[sel, "logFC"], 3),
                     P     = signif(tt[sel, "P.Value"], 3),
                     adjP  = signif(tt[sel, "adj.P.Val"], 3),
                     row.names = NULL, stringsAsFactors = FALSE)

# ------------------------- 2) time-point logFC (no tests) --------------------
time_fc <- sapply(c(P14 = "P14", P29 = "P29"), function(tp) {
  i_test <- meta$time == tp & meta$group == "test"
  i_ref  <- meta$time == tp & meta$group == "ref"
  round(rowMeans(X[, i_test, drop = FALSE]) - rowMeans(X[, i_ref, drop = FALSE]), 3)
})
timecourse <- data.frame(gene = pooled$gene, probe = pooled$probe,
                         logFC_P14 = time_fc[sel, "P14"],
                         logFC_P29 = time_fc[sel, "P29"],
                         row.names = NULL)

# ------------------- manuscript S9 values for comparison ---------------------
s9 <- data.frame(
  gene  = c("THBS1 (S9 primary)", "CD6", "CD3G", "PRTN3"),
  logFC = c(0.972, -0.437, -1.219, -0.136),
  P     = c(0.012, 0.016, 0.034, 0.171),
  adjP  = c(0.037, 0.048, 0.084, 0.284),
  P14   = c(1.34, -0.19, -0.74, -0.10),
  P29   = c(0.60, -0.68, -1.70, -0.17), stringsAsFactors = FALSE)

# --------------------------------- report ------------------------------------
cat("\n================ POOLED (P14+P29): computed vs S9 ================\n")
comp <- merge(pooled, s9, by = "gene", suffixes = c(".computed", ".S9"), sort = FALSE)
comp$logFC_ok <- abs(comp$logFC.computed - comp$logFC.S9) <= 0.005
comp$P_ok     <- abs(comp$P.computed     - comp$P.S9)     <= 0.002
comp$adjP_ok  <- abs(comp$adjP.computed  - comp$adjP.S9)  <= 0.002
print(comp[, c("gene", "probe", "logFC.computed", "logFC.S9", "logFC_ok",
               "P.computed", "P.S9", "P_ok", "adjP.computed", "adjP.S9", "adjP_ok")])

cat("\n================ TIME-POINT logFC: computed vs S9 ================\n")
tcomp <- merge(timecourse, s9, by = "gene", suffixes = c(".computed", ".S9"), sort = FALSE)
tcomp$P14_ok <- abs(tcomp$logFC_P14 - tcomp$P14) <= 0.01
tcomp$P29_ok <- abs(tcomp$logFC_P29 - tcomp$P29) <= 0.01
print(tcomp[, c("gene", "probe", "logFC_P14", "P14", "P14_ok", "logFC_P29", "P29", "P29_ok")])

write.csv(pooled,     file.path(OUT_DIR, "GSE25286_pooled_limma_results.csv"),     row.names = FALSE)
write.csv(timecourse, file.path(OUT_DIR, "GSE25286_timecourse_logFC_results.csv"), row.names = FALSE)
cat("\nCSV results written to:\n  ", OUT_DIR, "\n", sep = "")

# Interpretation hint
cat("\nNOTE: logFC (mean difference on log2 scale) must match S9 exactly.\n")
cat("P / adjP should be very close; small deviations indicate a slightly different\n")
cat("preprocessing/annotation version of the matrix used to generate S9.\n")
