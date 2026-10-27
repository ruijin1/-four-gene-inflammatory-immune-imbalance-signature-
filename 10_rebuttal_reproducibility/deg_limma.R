# limma moderated t-test for the day-14 discovery matrix (GSE32472)
# ASCII paths only (Windows R locale cannot parse CJK paths).
suppressMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
tmp <- args[1]

e <- as.matrix(read.csv(file.path(tmp, "deg_expr.csv"), row.names = 1, check.names = FALSE))
l <- read.csv(file.path(tmp, "deg_lab.csv"), row.names = 1, check.names = FALSE)
l <- l[colnames(e), , drop = FALSE]
g <- factor(l$grp, levels = c(0, 1))

design <- model.matrix(~ g)
colnames(design) <- c("Intercept", "MSBPD")

fit <- lmFit(e, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = "MSBPD", number = Inf, adjust.method = "BH", sort.by = "P")

out <- data.frame(Gene = rownames(tt), logFC = tt$logFC,
                  P = tt$P.Value, FDR = tt$adj.P.Val, row.names = NULL)
write.csv(out, file.path(tmp, "deg_limma.csv"), row.names = FALSE)

cat("n genes tested:", nrow(out), "\n")
cat("FDR<0.05 & |logFC|>0.5 :", sum(out$FDR < 0.05 & abs(out$logFC) > 0.5), "\n")
cat("FDR<0.05               :", sum(out$FDR < 0.05), "\n")
cat("|logFC|>0.5            :", sum(abs(out$logFC) > 0.5), "\n")
four <- c("CD3G", "PRTN3", "CD6", "THBS1")
cat("--- four signature genes ---\n")
print(out[out$Gene %in% four, ])
