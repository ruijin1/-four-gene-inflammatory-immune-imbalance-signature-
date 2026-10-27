# Targeted reproduction of the author's original SVM-RFE AvgRank on the authentic
# 28-gene x 63-sample matrix (bit-identical to the author's SVM.csv input).
tmp <- "D:/wb_ascii_tmp"
Xall <- as.matrix(read.csv(file.path(tmp, "svmrfe_orig_expr.csv"), row.names = 1, check.names = FALSE))
lab  <- read.csv(file.path(tmp, "svmrfe_orig_lab.csv"), stringsAsFactors = FALSE)
yy   <- factor(lab$grp, levels = c("Normal", "MS BPD"))
tgt  <- read.csv(file.path(tmp, "svmrfe_target_avgrank.csv"), stringsAsFactors = FALSE)
TARGET8 <- c("CD3G","THBS1","CD6","BCL11B","CD3E","PRTN3","CD27","ALOX5")
suppressMessages(library(e1071))
cat("X:", dim(Xall), " levels:", levels(yy), "\n\n")

# ranking within one training set: eliminate worst (smallest w^2) first
fold_rank <- function(X, y, cost, scale) {
  genes <- colnames(X); elim <- character(0)
  while (length(genes) > 1) {
    m  <- svm(X[, genes, drop = FALSE], y, kernel = "linear", cost = cost, scale = scale)
    w  <- as.numeric(t(m$coefs) %*% m$SV)
    k  <- which.min(w^2)
    elim <- c(elim, genes[k]); genes <- genes[-k]
  }
  elim <- c(elim, genes)
  r <- setNames(rev(seq_along(elim)), elim)   # last survivor = rank 1
  r[colnames(X)]
}

avg_rank <- function(X, y, k, cost, scale, seed) {
  n <- nrow(X)
  if (is.na(seed)) {
    f <- (seq_len(n) - 1) %% k + 1
  } else {
    set.seed(seed); f <- sample((seq_len(n) - 1) %% k + 1)
  }
  R <- matrix(NA, nrow = k, ncol = ncol(X), dimnames = list(NULL, colnames(X)))
  for (i in seq_len(k)) {
    te <- which(f == i)
    Xtr <- X[-te, , drop = FALSE]; ytr <- y[-te]
    R[i, ] <- fold_rank(Xtr, ytr, cost, scale)
  }
  colMeans(R)
}

res <- list()
for (k in c(5, 10)) {
  for (cost in c(1, 0.1, 10)) {
    for (scale in c(TRUE, FALSE)) {
      for (seed in c(NA, 1, 2, 3, 123)) {
        a <- avg_rank(Xall, yy, k, cost, scale, seed)
        top8 <- names(sort(a))[1:8]
        hit  <- sum(top8 %in% TARGET8)
        rho  <- cor(a[tgt$gene], tgt$AvgRank, method = "spearman")
        res[[length(res) + 1]] <- data.frame(k = k, cost = cost, scale = scale,
          seed = ifelse(is.na(seed), "seq", seed), hit8 = hit, spearman = round(rho, 3),
          top8 = paste(top8, collapse = " "), stringsAsFactors = FALSE)
      }
    }
  }
}
R <- do.call(rbind, res)
R <- R[order(-R$hit8, -R$spearman), ]
cat("=== ranked by top-8 hits then Spearman vs author AvgRank ===\n")
print(head(R, 12), row.names = FALSE)
cat("\nmax top8 hits:", max(R$hit8), " max spearman:", max(R$spearman), "\n\n")

best <- R[which.max(R$spearman), ]
k <- best$k; cost <- best$cost; scale <- best$scale
seed <- if (best$seed == "seq") NA else as.integer(best$seed)
a <- avg_rank(Xall, yy, k, cost, scale, seed)
cat("BEST CONFIG: k=", k, " cost=", cost, " scale=", scale, " seed=", best$seed, "\n", sep = "")
cmp <- data.frame(gene = names(a), reproduced = round(a, 1), author = tgt$AvgRank[match(names(a), tgt$gene)])
cmp <- cmp[order(cmp$reproduced), ]
print(cmp, row.names = FALSE)
cat("\nmean abs deviation from author AvgRank:",
    round(mean(abs(cmp$reproduced - cmp$author)), 2), "\n")
write.csv(cmp, file.path(tmp, "svmrfe_avgrank_best.csv"), row.names = FALSE)
