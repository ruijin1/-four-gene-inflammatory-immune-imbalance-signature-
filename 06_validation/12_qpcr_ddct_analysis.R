# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 06_validation/12_qpcr_ddct_analysis.R
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
# 12  A549 hyperoxia/normoxia qRT-PCR: 2^(-ddCt) relative quantification + two-group t test
# Manuscript: Methods, "RNA isolation and quantitative real-time PCR"; Fig 5I / relevant panels of Fig 8
# input format: Ct table (samples x genes), three biological replicates x three technical replicates
# [RECONSTRUCTED] Rewritten from the written description in the Methods section, for code release (Zenodo/GitHub) and reviewer reproducibility.
# NOTE: not validated in the original authors' analysis environment. Verify locally against the existing result files before use.
# =============================================================================
library(ggplot2); library(dplyr)

ct <- read.csv("data/qPCR_Ct.csv")             # columns: sample, group (normoxia/hyperoxia), gene, Ct

## 1. Detection: Ct > 35 is below the limit of detection
ct <- ct %>% mutate(detected = Ct <= 35)

## 2. ddCt: beta-actin as reference gene, normoxia group as calibrator
dd <- ct %>%
  group_by(sample) %>%
  mutate(dCt = Ct[gene == "ACTB"] - Ct) %>%
  ungroup() %>%
  filter(gene != "ACTB") %>%
  group_by(gene, group) %>%
  mutate(mean_dCt_ctrl = mean(dCt[group == "normoxia"], na.rm = TRUE),
         ddCt = dCt - mean_dCt_ctrl,
         rel  = 2^(-ddCt)) %>%
  ungroup()

## 3. Between-group comparison: unpaired t test, p < 0.05
summ <- dd %>% group_by(gene, group) %>%
  summarise(mean = mean(rel), sd = sd(rel), n = n(), .groups = "drop")
p_t  <- sapply(unique(dd$gene), function(g)
          t.test(rel ~ group, data = dd[dd$gene == g, ])$p.value)
print(p_t)

p <- ggplot(summ, aes(group, mean, fill = group)) +
     geom_col() + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .2) +
     facet_wrap(~ gene, scales = "free_y") + theme_bw()
ggsave("results/qpcr_barplot.pdf", p, width = 6, height = 4)
# manuscript result: THBS1 significantly up-regulated under hyperoxia (P < 0.01); CD3G/CD6/PRTN3 not detected in A549 (Ct > 35)
