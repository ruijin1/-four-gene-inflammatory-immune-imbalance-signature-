# Reproduce the published LASSO / SVM-RFE selection with the original R toolchain.
# Input: mod28_expr.csv (28 authentic Cytoscape module genes x 63 samples), mod28_y.csv
#        deirg57_expr.csv (57 DE-IRGs) as the alternative LASSO input
# Output: lasso_repro.csv (per-seed selections), svmrfe28_ranking.csv

library(glmnet)
library(e1071)

tmp <- "D:/wb_ascii_tmp"
readx <- function(f) {
  d <- read.csv(file.path(tmp, f), row.names = 1, check.names = FALSE)
  as.matrix(d)   # rows = samples, cols = genes
}
four <- c("CD3G","PRTN3","CD6","THBS1")

run_lasso <- function(x, y, seeds, label) {
  out <- data.frame()
  for (s in seeds) {
    set.seed(s)
    cv <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
    for (rule in c("lambda.min","lambda.1se")) {
      co <- as.numeric(coef(cv, s = rule))
      g  <- rownames(coef(cv))[which(co != 0)]
      g  <- setdiff(g, "(Intercept)")
      out <- rbind(out, data.frame(input = label, seed = s, rule = rule,
                                   n_sel = length(g),
                                   four_all = paste(sort(intersect(g, four)), collapse = ","),
                                   sel = paste(sort(g), collapse = ",")))
    }
  }
  out
}

x28 <- readx("mod28_expr.csv"); y28 <- as.factor(read.csv(file.path(tmp,"mod28_y.csv"))$y)
x57 <- readx("deirg57_expr.csv")
seeds <- 20250101:(20250101+99)
r28 <- run_lasso(x28, y28, seeds, "module28")
r57 <- run_lasso(x57, y28, seeds, "deirg57")
res <- rbind(r28, r57)
write.csv(res, file.path(tmp,"lasso_repro.csv"), row.names = FALSE)

cat("== LASSO (glmnet, 100 seeds) ==\n")
for (lb in c("module28","deirg57")) for (rule in c("lambda.min","lambda.1se")) {
  s <- res[res$input==lb & res$rule==rule,]
  exact4 <- sum(s$sel == paste(sort(four),collapse=","))
  cat(sprintf("%-9s %-10s : n_sel median=%g range=[%g,%g]  exact-4-gene sets=%d/100  four-genes-all-selected=%d/100\n",
      lb, rule, median(s$n_sel), min(s$n_sel), max(s$n_sel), exact4,
      sum(s$four_all == paste(sort(four),collapse=","))))
}
cat("\n-- module28 lambda.min: top 12 most frequent selections --\n")
tb <- sort(table(r28$rule=="lambda.min" & FALSE)) # placeholder
s <- r28[r28$rule=="lambda.min",]
fr <- sort(table(s$sel), decreasing = TRUE)
for (i in 1:min(12,length(fr))) cat(sprintf("  %3dx  %s\n", fr[i], names(fr)[i]))

# frequency of each of the four genes being retained (lambda.min, module28)
for (g in four) cat(sprintf("  %s retained in %d/100 (module28, lambda.min)\n", g, sum(grepl(g, s$sel))))

# ---- SVM-RFE (classic Guyon RFE, linear kernel e1071) on the 28 module genes ----
svmrfe <- function(X, y) {
  genes <- colnames(X); feats <- genes
  while (length(feats) > 1) {
    m <- svm(X[, feats, drop=FALSE], y, kernel="linear", cost=1, scale=TRUE)
    w <- colSums(t(m$coefs) * m$X.scale$`scaled:center` / m$X.scale$`scaled:scale`) # not needed; use weight via t(coef)
    wt <- as.numeric(t(m$coefs) %*% m$SV)          # hyperplane normal in scaled space? use standard:
    # simpler: ranks by |w|^2 using computed weights on scaled data
    sc <- scale(X[, feats, drop=FALSE])
    m2 <- svm(sc, y, kernel="linear", cost=1)
    w2 <- drop(t(m2$coefs) %*% m2$SV)
    rnk <- order(w2^2)                              # smallest |w| first = least important
    drop_g <- feats[rnk[1]]
    if (length(feats) > 20 && length(feats) - 1 >= 1) drop_n <- 1 else drop_n <- 1
    feats <- setdiff(feats, drop_g)
  }
  invisible(NULL)
}

# proper RFE returning full elimination order (drop one at a time; 28 genes is cheap)
svmrfe_order <- function(X, y) {
  feats <- colnames(X)
  order_out <- c()
  while (length(feats) > 1) {
    sc <- scale(X[, feats, drop = FALSE])
    m <- svm(sc, y, kernel = "linear", cost = 1)
    w <- drop(t(m$coefs) %*% m$SV)
    least <- names(which.min(w^2))
    order_out <- c(order_out, least)
    feats <- setdiff(feats, least)
  }
  c(rev(order_out), feats)   # most important first
}

ord <- svmrfe_order(x28, y28)
cat("\n== SVM-RFE (e1071 linear, one-at-a-time) top 10 ==\n")
cat(paste(head(ord, 10), collapse = ", "), "\n")
write.csv(data.frame(rank = seq_along(ord), gene = ord),
          file.path(tmp,"svmrfe28_ranking.csv"), row.names = FALSE)
