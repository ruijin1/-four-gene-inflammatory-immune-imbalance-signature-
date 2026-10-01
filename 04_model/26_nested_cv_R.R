# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 04_model/26_nested_cv_R.R
#  Status : ORIGINAL - re-run and validated (R 4.4.3)
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
# 26  Nested cross-validation of the four-gene signature (selection inside folds)
# -----------------------------------------------------------------------------
# Purpose: quantify the optimism introduced by the feature-selection step itself.
#          The original bootstrap optimism correction and LOOCV only refit the
#          final four-gene logistic model; the LASSO + SVM-RFE selection was not
#          repeated, so those internal estimates are optimistically biased.
#
# Three nested variants (outer loop = leave-one-out, 63 folds):
#   V0  refit the four-gene logistic model only (sanity check; expect ~0.873)
#   V1  redo LASSO + SVM-RFE with the feature pool fixed at the 57 DE-IRGs
#   V2  redo differential expression (limma eBayes) + LASSO + SVM-RFE
#
# Engines match the published 09 script:
#   LASSO  = glmnet::cv.glmnet(alpha = 1, nfolds = 10, s = "lambda.min")
#   SVM    = e1071::svm(kernel = "linear", cost = 1, scale = TRUE)
#   SVM-RFE = |w| recursive elimination ranking, then 10-fold CV over the
#             subset-size grid (this is what caret::rfe(method = "svmLinear",
#             rfeControl = rfeControl(method = "cv", number = 10), sizes = ...)
#             does; caret is not installed here, so it is implemented directly).
#
# NOTE: this script is pure ASCII on purpose. When R is launched from a shell
#       whose locale is not UTF-8, R cannot resolve non-ASCII paths, so inputs
#       are staged in an ASCII directory. Run with:
#   Rscript nested_cv_R.R
# =============================================================================
# STATUS (2026-09-26): RUN AND VALIDATED on R 4.4.3.
#   Implementation check: the four genes' training-set means match the Python run to
#   three decimals, the apparent four-gene AUC is 0.9284, and V0 reproduces the
#   reported leave-one-out AUC of 0.873. Results: V1 = 0.870, V2 = 0.794.
suppressPackageStartupMessages({
  library(readxl); library(limma); library(glmnet); library(e1071); library(writexl)
})

IN  <- P("data/derived/nested_cv_inputs")
OUT <- P("results/nested_cv")
SEED <- 20250101
SIZES <- c(1:20, 30, 40, 50)
DE_ADJP <- 0.05
DE_LFC  <- 0.58
TRUE4 <- c("CD3G", "PRTN3", "CD6", "THBS1")

## ---------- 1. data ----------
df <- as.data.frame(read_excel(file.path(IN, "expr.xlsx"), sheet = "S1"))
genes <- as.character(df[[1]])
m <- as.matrix(df[, -1]); storage.mode(m) <- "double"; rownames(m) <- genes
IRG  <- as.character(read_excel(file.path(IN, "irg.xlsx"))[[1]])
DE57 <- as.character(read_excel(file.path(IN, "de57.xlsx"))[[1]])
gp <- as.data.frame(read_excel(file.path(IN, "group.xlsx"), sheet = "norma and msbpd"))
smap <- setNames(as.character(gp[[2]]), as.character(gp[[1]]))
cols <- colnames(m)[colnames(m) %in% names(smap)]
m <- m[, cols]
y <- factor(ifelse(smap[cols] == "MS BPD", "MS-BPD", "control"),
            levels = c("control", "MS-BPD"))
cat(sprintf("matrix: %d genes x %d samples (%d cases / %d controls)\n",
            nrow(m), ncol(m), sum(y == "MS-BPD"), sum(y == "control")))

## ---------- 2. helpers ----------
auc_manual <- function(y, p) {
  yv <- as.integer(y == levels(y)[2]); o <- order(p)
  r <- rank(p[o], ties.method = "average"); n1 <- sum(yv); n0 <- length(yv) - n1
  (sum(r[yv[o] == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

strat_folds <- function(y, k, seed) {
  set.seed(seed); idx <- seq_along(y); f <- integer(length(y))
  for (lv in levels(y)) {
    ii <- sample(idx[y == lv]); f[ii] <- rep(seq_len(k), length.out = length(ii))
  }
  f
}

cv_acc <- function(x, y, k = 10, seed = SEED) {
  f <- strat_folds(y, k, seed); acc <- numeric(0)
  for (i in seq_len(k)) {
    tr <- f != i; te <- f == i
    if (length(unique(y[tr])) < 2 || sum(te) == 0) next
    fit <- svm(x[tr, , drop = FALSE], y[tr], kernel = "linear", cost = 1, scale = TRUE)
    acc <- c(acc, mean(predict(fit, x[te, , drop = FALSE]) == y[te]))
  }
  if (!length(acc)) return(0)
  mean(acc)
}

svm_rfe_rank <- function(x, y) {
  feats <- colnames(x); rk <- character(0)
  while (length(feats) > 1) {
    fit <- svm(x[, feats, drop = FALSE], y, kernel = "linear", cost = 1, scale = TRUE)
    w <- abs(as.numeric(t(fit$coefs) %*% fit$SV))
    drop <- feats[which.min(w)]
    rk <- c(drop, rk); feats <- setdiff(feats, drop)
  }
  c(feats, rk)
}

svm_rfe_select <- function(x, y) {
  if (ncol(x) <= 1) return(colnames(x))
  ord <- svm_rfe_rank(x, y)
  grid <- SIZES[SIZES <= length(ord)]; if (!length(grid)) grid <- length(ord)
  best <- list(size = grid[1], acc = -1)
  for (s in grid) {
    a <- cv_acc(x[, ord[1:s], drop = FALSE], y, 10)
    if (a > best$acc + 1e-12) best <- list(size = s, acc = a)
  }
  list(genes = ord[1:best$size], size = best$size)
}

lasso_select <- function(x, y) {
  set.seed(SEED)
  lam <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
  co <- coef(lam, s = "lambda.min")
  setdiff(rownames(co)[which(as.numeric(co) != 0)], "(Intercept)")
}

de_irgs <- function(X, yy) {
  design <- model.matrix(~ yy)
  fit <- eBayes(lmFit(X, design))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  rownames(tt)[tt$adj.P.Val < DE_ADJP & abs(tt$logFC) > DE_LFC]
}

## ---------- 2b. data alignment check (must match the Python run) ----------
cat("colnames[1:5]:", paste(head(colnames(m), 5), collapse = ", "), "\n")
cat("gene means:", paste(sprintf("%s=%.3f", TRUE4, rowMeans(m[TRUE4, , drop = FALSE])),
                        collapse = " "), "\n")
cat("(expect CD3G=10.207 PRTN3=8.798 CD6=7.998 THBS1=10.727)\n")
x_all <- t(m[TRUE4, , drop = FALSE])
fit_all <- glm(outcome ~ ., data = data.frame(x_all, outcome = y), family = binomial)
cat(sprintf("apparent AUC (all 63, four genes) = %.4f  (expect 0.9284)\n",
            auc_manual(y, predict(fit_all, type = "response"))))

## ---------- 3. leave-one-out ----------
run_nested <- function(mode) {
  n <- length(y); prob <- rep(NA_real_, n)
  cnt <- list(); de_n <- integer(0); hub_n <- integer(0)
  for (i in seq_len(n)) {
    tr <- setdiff(seq_len(n), i)
    if (mode == "V0") {
      pool <- TRUE4
    } else {
      pool <- if (mode == "V1") DE57 else intersect(de_irgs(m[, tr, drop = FALSE], y[tr]), IRG)
      pool <- intersect(pool, rownames(m))
      if (length(pool) > 200) pool <- pool[seq_len(200)]
      de_n <- c(de_n, length(pool))
    }
    x_tr <- t(m[pool, tr, drop = FALSE]); x_te <- t(m[pool, i, drop = FALSE])
    if (ncol(x_tr) < 2) { prob[i] <- 0.5; next }
    if (mode == "V0") {
      hub <- pool
    } else {
      gl <- lasso_select(x_tr, y[tr]); gs <- svm_rfe_select(x_tr, y[tr])$genes
      hub <- intersect(gl, gs); if (!length(hub)) hub <- gs[1]
    }
    hub_n <- c(hub_n, length(hub))
    for (g in hub) cnt[[g]] <- if (is.null(cnt[[g]])) 1 else cnt[[g]] + 1
    dd <- data.frame(x_tr[, hub, drop = FALSE])
    dd$outcome <- y[tr]                      # NOTE: do NOT write y[tr] inside the formula:
    fit <- glm(outcome ~ ., data = dd,        # with a column named y present in data, the
               family = binomial)            # term y[tr] would be re-subset and misalign.
    prob[i] <- predict(fit, newdata = as.data.frame(x_te[, hub, drop = FALSE]),
                       type = "response")
  }
  aucv <- auc_manual(y, prob)
  tbl <- sort(unlist(cnt), decreasing = TRUE)
  cat(sprintf("\n[%s] AUC = %.4f\n", mode, aucv))
  if (length(de_n)) cat(sprintf("  fold feature pool: median %g (%d-%d); median hub size %g\n",
                               median(de_n), min(de_n), max(de_n), median(hub_n)))
  cat("  true four selected in N folds:",
      paste(sprintf("%s=%d", TRUE4,
                    sapply(TRUE4, function(g) if (is.null(cnt[[g]])) 0 else cnt[[g]])),
            collapse = ", "), "\n")
  cat("  top genes:", paste(sprintf("%s x%d", head(names(tbl), 10), head(tbl, 10)),
                           collapse = ", "), "\n")
  list(auc = aucv, prob = prob, counts = tbl, de_n = de_n, hub_n = hub_n)
}

res <- lapply(c("V0", "V1", "V2"), run_nested); names(res) <- c("V0", "V1", "V2")

cmp <- data.frame(
  Analysis = c("Apparent (training, fixed 4 genes)",
               "Bootstrap optimism-corrected (fixed 4 genes)",
               "LOOCV, coefficients refitted (fixed 4 genes) [V0]",
               "Nested: LASSO + SVM-RFE redone on 57 DE-IRGs [V1]",
               "Nested: DE + LASSO + SVM-RFE redone [V2]",
               "External GSE220135 day 14 (frozen)",
               "External GSE220135 day 28 (frozen)"),
  AUC = round(c(0.9284, 0.902, res$V0$auc, res$V1$auc, res$V2$auc, 0.795, 0.769), 4))
print(cmp, row.names = FALSE)

if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
write_xlsx(list(`Performance comparison` = cmp,
                `Per-sample predictions` = data.frame(sample = cols, group = as.character(y),
                                                      V0 = res$V0$prob, V1 = res$V1$prob,
                                                      V2 = res$V2$prob)),
           file.path(OUT, "nested_cv_R_output.xlsx"))
cat("\nwritten:", file.path(OUT, "nested_cv_R_output.xlsx"), "\n")
