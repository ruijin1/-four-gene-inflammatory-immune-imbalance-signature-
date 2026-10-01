# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 04_model/09_lasso_svmrfe_roc.R
#  Status : RECONSTRUCTED from the Methods section - NOT validated in the original environment
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
# 09  LASSO (glmnet 4.1-7) + SVM-RFE + ROC (pROC)
# Manuscript: Methods, "Identification of diagnostic biomarkers"; Fig 5A-C
# [RECONSTRUCTED] Rewritten from the written description in the Methods section, for code release (Zenodo/GitHub) and reviewer reproducibility.
# NOTE: not validated in the original authors' analysis environment. Verify locally against the existing result files before use.
# =============================================================================
library(glmnet); library(e1071); library(caret); library(pROC)

expr  <- t(as.matrix(read.csv("data/GSE32472_day14_expr.csv", row.names = 1)))
grp   <- read.csv("data/GSE32472_day14_group.csv")
y     <- factor(grp$group, levels = c("control", "MS-BPD"))
x     <- expr[, read.csv("results/DE_IRG_57.csv")$DE_IRG]     # the 57 DE-IRGs

## ---- LASSO: 10-fold cross-validation for lambda selection ----
set.seed(20250101)
lam <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
co  <- coef(lam, s = "lambda.min")
gene_lasso <- rownames(co)[which(as.numeric(co) != 0)][-1]
pdf("results/lasso_cv.pdf"); plot(lam); dev.off()

## ---- SVM-RFE ----
ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 10)
rfe1 <- rfe(x, y, sizes = c(1:20, 30, 40, 50), rfeControl = ctrl, method = "svmLinear")
gene_svm <- predictors(rfe1)

## ---- Intersection = the four consensus genes ----
hub <- intersect(gene_lasso, gene_svm)          # expected: CD3G / PRTN3 / CD6 / THBS1

## ---- ROC / AUC / sensitivity / specificity (Youden-optimal cut-off) ----
fit  <- glm(as.formula(paste("y ~", paste(hub, collapse = "+"))), data = data.frame(x, y), family = binomial)
prob <- predict(fit, type = "response")
roc1 <- roc(y, prob, levels = c("control", "MS-BPD"))
ci.auc(roc1); coords(roc1, "best", ret = c("sensitivity", "specificity"))
pdf("results/roc_four_gene.pdf"); plot(roc1, print.auc = TRUE); dev.off()
