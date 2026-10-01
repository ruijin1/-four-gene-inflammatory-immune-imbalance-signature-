# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 04_model/10_nomogram_calibration_dca.R
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
# 10  Nomogram (rms 6.4-0) + calibration curve + DCA (exploratory)
# Manuscript: Methods, "Construction of the nomogram"; Fig 5D-H
# NOTE: net benefit depends on outcome prevalence; this is a case-control design, so the DCA is reported as exploratory only
# [RECONSTRUCTED] Rewritten from the written description in the Methods section, for code release (Zenodo/GitHub) and reviewer reproducibility.
# NOTE: not validated in the original authors' analysis environment. Verify locally against the existing result files before use.
# =============================================================================
library(rms); library(rmda); library(ggplot2)

expr <- t(as.matrix(read.csv("data/GSE32472_day14_expr.csv", row.names = 1)))
grp  <- read.csv("data/GSE32472_day14_group.csv")
hub  <- c("CD3G", "PRTN3", "CD6", "THBS1")
d    <- data.frame(expr[, hub], group = as.integer(grp$group == "MS-BPD"))

dd <- datadist(d); options(datadist = "dd")
f  <- lrm(group ~ CD3G + PRTN3 + CD6 + THBS1, data = d, x = TRUE, y = TRUE)
print(f)                                        # coefficients: -3.375, PRTN3 +0.610, CD3G -0.820, CD6 -1.402, THBS1 +1.560

## Nomogram
pdf("results/nomogram.pdf")
plot(nomogram(f, fun = plogis, funlabel = "Predicted probability of MS-BPD"))
dev.off()

## Calibration curve
pdf("results/calibration.pdf"); plot(calibrate(f, B = 1000)); dev.off()

## DCA (exploratory)
dca <- decision_curve(group ~ CD3G + PRTN3 + CD6 + THBS1, data = d,
                      thresholds = seq(0, 0.8, by = 0.01), bootstraps = 1000)
pdf("results/dca.pdf"); plot_decision_curve(dca); dev.off()
