# Robustness of the 8-gene set under the identified configuration
# (20 resamples x 80% subsample, linear SVM cost 0.05, standardized, rank averaged)
tmp <- "D:/wb_ascii_tmp"
Xall <- as.matrix(read.csv(file.path(tmp, "svmrfe_orig_expr.csv"), row.names = 1, check.names = FALSE))
lab  <- read.csv(file.path(tmp, "svmrfe_orig_lab.csv"), stringsAsFactors = FALSE)
yy   <- factor(lab$grp, levels = c("Normal", "MS BPD"))
tgt  <- read.csv(file.path(tmp, "svmrfe_target_avgrank.csv"), stringsAsFactors = FALSE)
TARGET8 <- c("CD3G","THBS1","CD6","BCL11B","CD3E","PRTN3","CD27","ALOX5")
suppressMessages(library(e1071))

fold_rank <- function(X, y, cost, scale) {
  genes <- colnames(X); elim <- character(0)
  while (length(genes) > 1) {
    m <- svm(X[, genes, drop = FALSE], y, kernel = "linear", cost = cost, scale = scale)
    w <- as.numeric(t(m$coefs) %*% m$SV); k <- which.min(w^2)
    elim <- c(elim, genes[k]); genes <- genes[-k]
  }
  elim <- c(elim, genes); setNames(rev(seq_along(elim)), elim)[colnames(X)]
}
bag_rank <- function(frac, nsub, cost, scale, seed) {
  set.seed(seed); n <- nrow(Xall)
  R <- matrix(NA, nsub, ncol(Xall), dimnames = list(NULL, colnames(Xall)))
  for (i in seq_len(nsub)) {
    ix <- sample(n, floor(frac * n))
    R[i, ] <- fold_rank(Xall[ix, , drop = FALSE], yy[ix], cost, scale)
  }
  colMeans(R)
}
out <- list()
for (frac in c(0.7, 0.8, 0.9)) {
  for (seed in 1:10) {
    a <- bag_rank(frac, 20, 0.05, TRUE, seed)
    t8 <- names(sort(a))[1:8]
    out[[length(out)+1]] <- data.frame(frac = frac, seed = seed,
      hit8 = sum(t8 %in% TARGET8), exact = setequal(t8, TARGET8),
      spearman = round(cor(a[tgt$gene], tgt$AvgRank, method = "spearman"), 3),
      top8 = paste(t8, collapse = " "), stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, out)
cat("=== per-config ===\n"); print(R[, 1:5], row.names = FALSE)
cat("\n=== summary ===\n")
agg <- aggregate(cbind(hit8, exact, spearman) ~ frac, data = R,
                 FUN = function(v) round(mean(v), 2))
print(agg, row.names = FALSE)
cat("\nexact-set recovery rate overall:", mean(R$exact), " (", sum(R$exact), "/", nrow(R), ")\n")
cat("mean hits:", round(mean(R$hit8), 2), "/8\n")
write.csv(R, file.path(tmp, "svmrfe_robustness.csv"), row.names = FALSE)
