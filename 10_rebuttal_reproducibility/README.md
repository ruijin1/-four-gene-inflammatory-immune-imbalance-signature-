# 10 · Reproducibility scripts for the reviewer response (Oct 2026)

This folder was added for the first revision. It contains the minimal scripts and
inputs needed to reproduce the numbers quoted in the response letter and in
**Supplementary Data S17 and S18**. Nothing here changes the analysis reported in
the original submission; every script re-runs an already-reported result.

The folder was first assembled for the gene-set sensitivity analysis (Data S17)
and later extended with the two analyses added in response to Reviewer 2
(GSE151974 hyperoxia-lung single-cell mapping) and to Reviewer 1 (whole-blood
marker-set check, Data S18).

## Scripts

| Script | Reproduces | Notes |
|---|---|---|
| `deg_limma.R` | 139 DEGs (84 up / 55 down) from GSE32472 day-14 | limma moderated t, Benjamini–Hochberg adjusted *P* < 0.05 and \|log2FC\| > 0.58; 57 DE-IRGs after intersecting the 3,336 inflammatory-response genes |
| `lasso_by_definition.R` | LASSO selection under the manuscript definition and four alternatives (five in total) | glmnet, `lambda.min`, 10-fold CV; modal solution across 100 seeds |
| `lasso_svm_repro.R` | LASSO + SVM-RFE on the 28-gene CytoHubba module | reproduces the four-gene signature CD3G / CD6 / PRTN3 / THBS1 |
| `svmrfe_avgrank_repro.R` | SVM-RFE average-rank aggregation | 20 resamples at 80%, linear kernel, `cost = 0.05`, standardised predictors |
| `svmrfe_avgrank_grid2.R` | grid check of the resample/cost settings used above | sensitivity of the retained gene set |
| `svmrfe_cv_repro.R` | SVM-RFE under cross-validation variants | comparison only; the resampling version is the reported one |
| `svmrfe_robustness.R` | stability of the top-8 gene set across 30 runs | reports per-gene selection frequency |
| `svmrfe_targeted_repro.R` | targeted check of individual genes (e.g. PRTN3) | |
| `make_supplementary_S17.py` | builds `Data S17.xlsx` from the CSV outputs | requires `openpyxl`, `pandas` |
| `scBPD_4gene_lung_mapping.py` | cellular localization of *Cd3g* / *Cd6* / *Thbs1* / *Prtn3* in hyperoxia-exposed newborn mouse lung (Reviewer 2, major 1) | GSE151974; 61,839 cells at P3 / P7 / P14; the 34 author-supplied cell-type annotations are used as supplied, without re-clustering, so the analysis is descriptive |
| `decompose_pseudobulk.py` | decomposition of the bulk *Cd3g* decline into cell-composition versus per-cell expression components (Reviewer 2, major 3) | shows the bulk change is attributable to a smaller T-cell fraction rather than to lower per-cell expression |
| `blood_marker_check.py` | deconvolution-free marker-set test for lung-derived contamination of the whole-blood signal (Reviewer 1, major 4) | run on the GSE32472 day-14 matrix; T-cell, neutrophil, platelet, immature-granulocyte, airway/alveolar epithelial and endothelial sets, reported as standardised differences with BH-adjusted FDR |
| `make_S18.py` | builds `Data S18.xlsx` from the CSV outputs | marker-set analysis plus the four GSE151974 summary tables; requires `pandas`, `openpyxl` |

## Input data

- `cytohubba_module28_from_cytoscape.csv` — the 28 genes exported from CytoHubba
  (intersection of MCC / EPC / Degree top-30 lists).
- `author_confirmed_original_selections.csv` — the gene lists obtained with LASSO
  and SVM-RFE in the original analysis, as provided by the authors.
- `SVM.csv`, `SVM_feature_scoring_table.xlsx` — the original SVM input matrix and
  feature-scoring table from the first submission.
- `irg_sensitivity_v4_*.csv` — the outputs reproduced by the scripts above and
  tabulated in Supplementary Data S17.
- `marker_check_results.csv` — marker-set standardised differences and adjusted
  *P* values, tabulated as the first sheet of Supplementary Data S18.
- `composition_shift.csv`, `de_within_celltypes.csv`,
  `decomposition_contribution.csv`, `pseudobulk_log2FC.csv` — the GSE151974
  outputs tabulated as sheets 2–5 of Supplementary Data S18.
- `4gene_localisation.csv` — per-cell-type expression and detection proportions
  for the four genes in the hyperoxia-exposed newborn mouse lung.

## Environment

R >= 4.4 with `limma`, `glmnet`, `e1071`; Python >= 3.10 with `pandas`, `openpyxl`,
`scanpy`/`anndata` for `scBPD_4gene_lung_mapping.py`. Random seeds are set
explicitly inside each script; the LASSO modal solution and the SVM-RFE gene *set*
are stable across seeds, while individual average-rank values vary slightly with
the seed.

Public input datasets are GSE32472, GSE220135, GSE25286, GSE212098 and GSE151974;
none is redistributed here. Paths inside the scripts point at the authors' local
working directories and must be adjusted before re-running.
