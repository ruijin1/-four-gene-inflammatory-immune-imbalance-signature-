# Try to reproduce the manuscript's SVM-RFE 8 genes on the authentic 28-gene module.
# Variant: 10-fold CV aggregated ranking (caret-style), linear svm (e1071).
library(e1071)
tmp <- "D:/wb_ascii_tmp"
X <- as.matrix(read.csv(file.path(tmp,"mod28_expr.csv"), row.names = 1, check.names = FALSE))
y <- as.factor(read.csv(file.path(tmp,"mod28_y.csv"))$y)
target8 <- c("CD3G","THBS1","CD6","BCL11B","CD3E","PRTN3","CD27","ALOX5")

rfe_order_fold <- function(Xtr, ytr) {
  feats <- colnames(Xtr); out <- c()
  while (length(feats) > 1) {
    sc <- scale(Xtr[, feats, drop = FALSE])
    m <- svm(sc, ytr, kernel = "linear", cost = 1)
    w <- drop(t(m$coefs) %*% m$SV)
    least <- feats[which.min(w^2)]
    out <- c(out, least); feats <- setdiff(feats, least)
  }
  rev(c(out, feats))          # most important first
}

set.seed(20250101)
folds <- split(seq_len(nrow(X)), cut(seq_len(nrow(X)), breaks = 10))
if (any(sapply(folds, function(f) length(unique(y[f])) < 2))) {
  # stratified folds
  folds <- list()
  for (cl in levels(y)) {
    idx <- which(y == cl); folds[[cl]] <- split(idx, cut(seq_along(idx), 5))
  }
  fl <- list()
  for (k in 1:5) fl[[k]] <- unlist(sapply(folds, function(f) f[[k]]))
  folds <- fl
}

ranks <- matrix(0, nrow = ncol(X), ncol = length(folds),
                dimnames = list(colnames(X), NULL))
cv_err <- matrix(NA, nrow = 28, ncol = length(folds))
for (f in seq_along(folds)) {
  te <- folds[[f]]; tr <- setdiff(seq_len(nrow(X)), te)
  ord <- rfe_order_fold(X[tr, ], y[tr])
  ranks[, f] <- match(colnames(X), ord)
  for (k in 1:28) {
    gs <- ord[1:k]
    sc_tr <- scale(X[tr, gs, drop = FALSE], center = TRUE, scale = TRUE)
    m <- svm(sc_tr, y[tr], kernel = "linear", cost = 1)
    mu <- attr(sc_tr, "scaled:center"); sd <- attr(sc_tr, "scaled:scale")
    sc_te <- scale(X[te, gs, drop = FALSE], center = mu, scale = sd)
    pr <- predict(m, sc_te)
    cv_err[k, f] <- mean(pr != y[te])
  }
}
kstar <- which.min(rowMeans(cv_err))
mean_rank <- rowMeans(ranks)
topk <- names(sort(mean_rank))[1:kstar]   # smaller mean rank = more important
cat("best size k* =", kstar, " (CV err", round(min(rowMeans(cv_err)), 3), ")\n")
cat("top-k by mean rank:", paste(topk, collapse = ", "), "\n")
cat("match target8:", setequal(topk, target8), "| overlap:", length(intersect(topk, target8)), "/8\n")
cat("\nmean-rank order (first 12):\n")
print(head(sort(mean_rank), 12))
cat("\nCV error by size (k = 1..15):", paste(round(rowMeans(cv_err)[1:15], 3), collapse = ", "), "\n")
