# Try to reproduce the confirmed original SVM-RFE top-8 on the authentic 28-gene module
# Target (author-confirmed): CD3G THBS1 CD6 BCL11B CD3E PRTN3 CD27 ALOX5
tmp <- "D:/wb_ascii_tmp"
d  <- read.csv(file.path(tmp, "mod28_expr.csv"), row.names = 1, check.names = FALSE)
y  <- as.factor(read.csv(file.path(tmp, "mod28_y.csv"))$y)
X  <- as.matrix(d)
TARGET <- c("CD3G","THBS1","CD6","BCL11B","CD3E","PRTN3","CD27","ALOX5")

suppressMessages(library(e1071))

svm_rfe <- function(X, y, cost, scale) {
  genes <- colnames(X); rank <- character(0)
  while (length(genes) > 1) {
    m <- svm(X[, genes, drop = FALSE], y, kernel = "linear", cost = cost, scale = scale)
    w <- t(m$coefs) %*% m$SV          # normal vector
    sc <- as.numeric(w^2)
    k <- which.min(sc)
    rank <- c(rank, genes[k])
    genes <- genes[-k]
  }
  rev(c(rank, genes))                  # best first
}

cat("Target top8 :", paste(TARGET, collapse = " "), "\n\n")
for (cost in c(1, 0.5, 2, 10)) {
  for (scale in c(TRUE, FALSE)) {
    r <- svm_rfe(X, y, cost, scale)
    top8 <- r[1:8]
    hit <- sum(top8 %in% TARGET)
    cat("cost=%-4s scale=%-5s -> hits %d/8 | %s\n",
        cost, scale, hit, paste(top8, collapse = " "))
    if (identical(sort(top8), sort(TARGET))) cat("   *** EXACT MATCH ***\n")
  }
}

# also: remove-half-per-iteration variant (common speed-up), then refine
svm_rfe_half <- function(X, y, cost) {
  genes <- colnames(X); rank <- character(0)
  while (length(genes) > 1) {
    m <- svm(X[, genes, drop = FALSE], y, kernel = "linear", cost = cost, scale = TRUE)
    w <- t(m$coefs) %*% m$SV
    sc <- as.numeric(w^2)
    ord <- order(sc)                 # worst first
    n <- max(1, length(genes) %/% 2)
    drop <- ord[1:n]
    rank <- c(rank, genes[rev(drop)]) # keep order: better of dropped are later
    genes <- genes[-drop]
  }
  rev(c(rank, genes))
}
for (cost in c(1, 2)) {
  r <- svm_rfe_half(X, y, cost)
  top8 <- r[1:8]
  cat("half-variant cost=%s -> hits %d/8 | %s\n", cost, sum(top8 %in% TARGET), paste(top8, collapse = " "))
  if (identical(sort(top8), sort(TARGET))) cat("   *** EXACT MATCH ***\n")
}
