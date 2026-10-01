# Environment

## 1. Software versions

| | Version | Notes |
|---|---|---|
| R | 4.2.1 | original analysis environment: preprocessing, enrichment, model building, CIBERSORT |
| R | 4.4.3 | used to run `04_model/26_nested_cv_R.R` (nested cross-validation) |
| Python | 3.13 | all cohort-level validation, sensitivity, single-cell and HPA analyses |
| Cytoscape | 3.10.3 + CytoHubba | network figures (steps 7 and 8) — GUI, not scriptable |

## 2. R packages

| Package | Version (original) | Used in |
|---|---|---|
| limma | 3.54.1 | preprocessing, DEGs, single-gene GSEA, GSVA comparison, mouse re-analysis |
| clusterProfiler | 4.4.4 | GO / KEGG, ORA, GSEA |
| msigdbr | — | MSigDB C2.CP and Hallmark gene sets |
| GSVA | 1.44.3 | Hallmark GSVA |
| glmnet | 4.1-7 | LASSO, and the SVM-RFE reimplementation in the nested CV script |
| e1071 | — | SVM-RFE (linear kernel) |
| rms | 6.4-0 | nomogram, calibration |
| rmda | — | decision curve analysis |
| pROC | — | ROC / AUC |
| pheatmap | 1.0.12 | heatmaps |
| ggplot2 | 3.4.2 | all enrichment and model figures |
| dplyr | — | data handling |
| GEOquery | — | downloading GSE32472 |
| readxl / writexl | — | reading and writing `.xlsx` |
| VennDiagram | 1.8.2 | DEG ∩ IRG Venn diagram |
| CIBERSORT | 1.06 | immune deconvolution (needs the LM22 matrix, see DATA.md) |
| caret | — | **not installed** in the release environment. `09_lasso_svmrfe_roc.R` lists it, but `26_nested_cv_R.R` implements the equivalent SVM-RFE directly with `glmnet` + `e1071`, so the step can be run without it. |

## 3. Python packages

See [`../requirements.txt`](../requirements.txt). Core: numpy, pandas, scipy, scikit-learn,
matplotlib, statsmodels, openpyxl. The single-cell steps additionally need scanpy, anndata and
pydeseq2; figure assembly needs Pillow.

## 4. Pitfalls encountered — worth knowing before you run anything

**(a) `Rscript` on Windows cannot open source files whose path contains non-ASCII characters.**
It fails with `cannot open the connection`, and `dir.exists()` silently returns `FALSE` for such
paths. This is why all R source files here are ASCII-only (comments included), and why
`26_nested_cv_R.R` reads and writes to an ASCII working directory. RStudio is not affected. If your
clone lives under a non-ASCII path such as `D:/my workspace/...`, either run R from RStudio or move the
clone to an ASCII path (e.g. `C:/bpd-signature`).

**(b) `scanpy.pp.scrublet` requires `scikit-image`.** Without it (and with `threshold=None`) scanpy
raises `ValueError` rather than degrading gracefully. Install `scikit-image` if you want to run
`08_singlecell_hpa/42_scrublet_doublet_check.py` or `44_doublet_impact.py`.

**(c) In the nested cross-validation, the fold's response variable must be placed in the data frame
before fitting.** Writing

```r
glm(y[tr] ~ ., data = cbind(dd, y = y[tr]))     # WRONG
```

makes the formula look `y[tr]` up **in `data`** and subset it a second time, so the response and the
predictors become misaligned **without any error being raised** (the symptom was V0 = 0.8344 instead
of the correct 0.873). The correct form is

```r
dd$outcome <- y[tr]
glm(outcome ~ ., data = dd)                     # CORRECT
```

**(d) CIBERSORT requires the LM22 signature matrix** and a licence-gated download. It is not
redistributed here.

**(e) `statsmodels` has no `conf_int(method="profile")`.** Profile-likelihood intervals for the
logistic coefficients have to be obtained by constrained maximisation (this is how the confidence
intervals reported for the model coefficients were produced).
