# Grid 2: stratified (class-balanced) folds, finer cost grid; also 80% bagging variant
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
  elim <- c(elim, genes)
  setNames(rev(seq_along(elim)), elim)[colnames(X)]
}

strat_folds <- function(y, k, seed) {          # caret::createFolds style
  set.seed(seed)
  lapply(split(seq_along(y), y), function(ix) {
    f <- (seq_along(ix) - 1) %% k + 1
    split(ix, sample(f))
  }) -> per_class
  out <- integer(length(y))
  for (i in seq_along(per_class)) {
    idx <- as.integer(names(per_class[[i]])); out[idx] <- 1:k
  }
  out
}

run <- function(k, cost, scale, seed, mode) {
  n <- nrow(Xall)
  if (mode == "cv") {
    f <- strat_folds(yy, k, seed); kk <- k
  } else {                                     # bagging: k random 80% subsets
    set.seed(seed)
    bs <- lapply(seq_len(k), function(i) sample(n, floor(0.8 * n)))
    R <- matrix(NA, k, ncol(Xall), dimnames = list(NULL, colnames(Xall)))
    for (i in seq_len(k)) R[i, ] <- fold_rank(Xall[bs[[i]], , drop = FALSE], yy[bs[[i]]], cost, scale)
    return(colMeans(R))
  }
  R <- matrix(NA, kk, ncol(Xall), dimnames = list(NULL, colnames(Xall)))
  for (i in seq_len(kk)) {
    te <- which(f == i)
    R[i, ] <- fold_rank(Xall[-te, , drop = FALSE], yy[-te], cost, scale)
  }
  colMeans(R)
}

res <- list()
for (mode in c("cv", "bag")) {
  for (k in if (mode == "cv") c(5, 10) else c(20)) {
    for (cost in c(0.05, 0.1, 0.25, 0.5, 1, 2)) {
      for (scale in c(TRUE, FALSE)) {
        for (seed in 1:4) {
          a <- if (mode == "cv") run(k, cost, scale, seed, "cv") else run(20, cost, scale, seed, "bag")
          top8 <- names(sort(a))[1:8]
          res[[length(res)+1]] <- data.frame(mode = mode, k = k, cost = cost, scale = scale, seed = seed,
            hit8 = sum(top8 %in% TARGET8),
            spearman = round(cor(a[tgt$gene], tgt$AvgRank, method = "spearman"), 3),
            top8 = paste(top8, collapse = " "), stringsAsFactors = FALSE)
        }
      }
    }
  }
}
R <- do.call(rbind, res); R <- R[order(-R$hit8, -R$spearman), ]
cat("=== TOP 15 ===\n"); print(head(R, 15), row.names = FALSE)
cat("\nmax hit8:", max(R$hit8), " max spearman:", max(R$spearman), "\n")
ex <- R[R$hit8 == 8, ]
if (nrow(ex)) {
  cat("\n*** EXACT 8/8 CONFIGS: ", nrow(ex), " ***\n"); print(head(ex, 10), row.names = FALSE)
  b <- ex[which.max(ex$spearman), ]
  a <- if (b$mode == "cv") run(b$k, b$cost, b$scale, b$seed, "cv") else run(20, b$cost, b$scale, b$seed, "bag")
  cmp <- data.frame(gene = names(a), reproduced = round(a, 1), author = tgt$AvgRank[match(names(a), tgt$gene)])
  cmp <- cmp[order(cmp$reproduced), ]
  cat("\nBEST EXACT CONFIG:", b$mode, "k=", b$k, "cost=", b$cost, "scale=", b$scale, "seed=", b$seed, "\n")
  print(cmp, row.names = FALSE)
  cat("\nmean |dev|:", round(mean(abs(cmp$reproduced - cmp$author)), 2), "\n")
  write.csv(cmp, file.path(tmp, "svmrfe_avgrank_exact.csv"), row.names = FALSE)
}
