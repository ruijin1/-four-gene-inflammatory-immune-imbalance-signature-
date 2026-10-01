# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 01_data_prep/01_geo_download_preprocess_GSE32472.R
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
# 01  GSE32472 data acquisition and preprocessing (day-14 = training set / day-28 = temporal validation set)
# Manuscript: Methods, "Data collection and preprocessing"; Additional files 1 and 2
# [RECONSTRUCTED] Rewritten from the written description in the Methods section, for code release (Zenodo/GitHub) and reviewer reproducibility.
# NOTE: not validated in the original authors' analysis environment. Verify locally against the existing result files before use.
# =============================================================================
library(GEOquery)
library(limma)

set.seed(20250101)
dir.create("data", showWarnings = FALSE)

## ---- 1. Download and read the GPL6244 platform data (skip if already downloaded) ----
gse <- getGEO("GSE32472", GSEMatrix = TRUE, getGPL = TRUE)
eset <- gse[[1]]
expr  <- exprs(eset)              # probes x samples
pheno <- pData(eset)

## ---- 2. Annotate probes to gene symbols; for duplicated genes keep the probe with the highest mean expression ----
gpl <- getGEO("GPL6244", destdir = "data", AnnotGPL = TRUE)
annot <- Table(gpl)[, c("ID", "Gene symbol")]
annot <- annot[annot$`Gene symbol` != "" & !is.na(annot$`Gene symbol`), ]
annot$sym <- sub(" ///.*", "", annot$`Gene symbol`)

expr <- expr[rownames(expr) %in% annot$ID, ]
expr <- avereps(expr, ID = annot$sym[match(rownames(expr), annot$ID)])
expr <- log2(expr + 1)           # if the raw values are MAS5 signals, quantile-normalise after log2

## ---- 3. limma background correction + quantile normalisation (manuscript: limma 3.54.1 / R 4.2.1) ----
expr_norm <- normalizeBetweenArrays(expr, method = "quantile")

## ---- 4. Split by day of sampling into the day-14 training set and the day-28 temporal validation set ----
meta <- data.frame(
  sample = colnames(expr_norm),
  day    = pheno$`postnatal day:ch1`,
  group  = pheno$`disease state:ch1`
)
day14 <- meta$sample[meta$day == "14"]
day28 <- meta$sample[meta$day == "28"]

train_expr <- expr_norm[, day14];  train_grp <- meta$group[match(day14, meta$sample)]
val28_expr <- expr_norm[, day28];  val28_grp <- meta$group[match(day28, meta$sample)]

write.csv(train_expr,  "data/GSE32472_day14_expr.csv")
write.csv(val28_expr,  "data/GSE32472_day28_expr.csv")
write.csv(data.frame(sample = day14, group = train_grp), "data/GSE32472_day14_group.csv", row.names = FALSE)
write.csv(data.frame(sample = day28, group = val28_grp), "data/GSE32472_day28_group.csv", row.names = FALSE)
# expected: day-14 = 24 MS-BPD + 39 controls; day-28 = 24 MS-BPD + 38 controls
print(table(train_grp)); print(table(val28_grp))

## ---- 5. IRG list: GeneCards v5.22 + MSigDB v2024.1.Hs, search term "inflammatory response" ----
## merged and de-duplicated from the manual/web searches (Additional file 3; 3,336 IRGs)
# irg <- read.csv("data/IRG_union_3336.csv")$symbol
# intersection of the IRGs with the genes detected in the training set -> 1,498 (Additional file 5)
