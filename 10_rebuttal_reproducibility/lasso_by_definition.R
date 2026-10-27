# Per-definition LASSO sensitivity (R glmnet, original toolchain).
# Inputs (D:/wb_ascii_tmp):
#   py_modules.csv            - definitions + module gene lists (Python CytoHubba reimpl)
#   mod28_expr.csv            - authentic Cytoscape module expression (Def-1 anchor)
#   mod28_y.csv               - class labels
#   modpy_<i>_expr.csv        - per-definition module expression (i = row order of py_modules.csv)
# Outputs:
#   lasso_by_definition.csv   - full-data LASSO: modal selection + frequencies (100 seeds)
#   lasso_bootstrap.csv       - 100 resamples x definitions: per-gene selection frequency
library(glmnet)

tmp <- "D:/wb_ascii_tmp"
four <- c("CD3G","PRTN3","CD6","THBS1")
N_SEED <- 100; N_BOOT <- 100
SEED0 <- 20250101

defs <- read.csv(file.path(tmp,"py_modules.csv"), fileEncoding = "UTF-8-BOM", check.names = FALSE)
defs <- rbind(data.frame(def_id = "Def-1a", n_module = 28, module = ""), defs)
y    <- as.factor(read.csv(file.path(tmp,"mod28_y.csv"))$y)

getx <- function(id) {
  if (id == "Def-1a")
    as.matrix(read.csv(file.path(tmp,"mod28_expr.csv"), row.names = 1, check.names = FALSE))
  else {
    i <- sub("^Def-", "", id)
    as.matrix(read.csv(file.path(tmp,sprintf("modpy_%s_expr.csv", i)), row.names = 1, check.names = FALSE))
  }
}

sel_at_min <- function(x, y, seed) {
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
  co <- as.numeric(coef(cv, s = "lambda.min"))
  setdiff(rownames(coef(cv))[which(co != 0)], "(Intercept)")
}

# ---------- full data: modal selection over seeds ----------
rowsA <- list()
for (i in seq_len(nrow(defs))) {
  x <- getx(defs$def_id[i])
  # loop seeds
  sel_list <- lapply(SEED0:(SEED0+N_SEED-1), function(s) sel_at_min(x, y, s))
  sets <- sapply(sel_list, function(v) paste(sort(v), collapse = ","))
  tb <- sort(table(sets), decreasing = TRUE)
  modal <- names(tb)[1]
  fr4 <- sapply(four, function(gg) round(100*mean(sapply(sel_list, function(v) gg %in% v)), 1))
  rowsA[[i]] <- data.frame(definition = defs$def_id[i], n_module = ncol(x),
                           modal_set = modal, modal_freq_pct = round(100*tb[1]/N_SEED, 1),
                           n_modal_genes = length(strsplit(modal, ",")[[1]]),
                           four_all_in_modal = all(four %in% strsplit(modal, ",")[[1]]),
                           CD3G_pct = fr4["CD3G"], PRTN3_pct = fr4["PRTN3"],
                           CD6_pct = fr4["CD6"], THBS1_pct = fr4["THBS1"])
}
A <- do.call(rbind, rowsA)
write.csv(A, file.path(tmp,"lasso_by_definition.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("== Full-data LASSO (glmnet lambda.min, 100 seeds) ==\n")
print(A[, c("definition","n_module","modal_set","modal_freq_pct","four_all_in_modal",
            "CD3G_pct","PRTN3_pct","CD6_pct","THBS1_pct")], row.names = FALSE)

# ---------- bootstrap: per-gene selection frequency ----------
i1 <- which(y == 1); i0 <- which(y == 0)
cnt <- matrix(0, nrow = nrow(defs), ncol = length(four), dimnames = list(defs$def_id, four))
valid <- rep(0, nrow(defs))
for (b in 1:N_BOOT) {
  smp <- c(sample(i1, length(i1), TRUE), sample(i0, length(i0), TRUE))
  yb  <- y[smp]
  for (i in seq_len(nrow(defs))) {
    x <- getx(defs$def_id[i])[smp, , drop = FALSE]
    if (length(unique(yb)) < 2) next
    sel <- try(sel_at_min(x, yb, SEED0 + 1000 + b), silent = TRUE)
    if (inherits(sel, "try-error")) next
    valid[i] <- valid[i] + 1
    for (k in seq_along(four)) if (four[k] %in% sel) cnt[i, k] <- cnt[i, k] + 1
  }
}
B <- data.frame(definition = rownames(cnt), valid = valid,
                sapply(four, function(g) round(100*cnt[, g]/pmax(valid,1), 1)),
                check.names = FALSE)
B$four_all_pct_lower_bound <- round(100 * apply(cnt, 1, min) / pmax(valid, 1), 1)
write.csv(B, file.path(tmp,"lasso_bootstrap.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\n== Bootstrap (100 resamples, glmnet lambda.min) per-gene selection % ==\n")
print(B, row.names = FALSE)
