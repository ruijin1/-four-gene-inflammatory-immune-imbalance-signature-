# Data sources and directory layout

The input data are **not** redistributed with this repository. Every dataset is either a public GEO
deposit or a third-party reference resource. Download them and arrange them as described below.

## 1. Public datasets

| Dataset | Accession | Type | Used for |
|---|---|---|---|
| Discovery (training) cohort | **GSE32472** | Affymetrix GPL6244 microarray, whole blood | day-14 training set (24 MS-BPD / 39 controls); day-28 temporal window (24 / 38); case-definition sensitivity analyses |
| Independent external cohort | **GSE220135** | RNA-seq, whole blood | external validation at postnatal days 14 and 28 |
| Mouse hyperoxia model | **GSE25286** | lung tissue transcriptome | cross-species validation (Table S5) |
| Rat hyperoxia model | **GSE212098** | lung tissue, postnatal days 1-14 | independent cross-species validation (Table S8) |
| Additional neonatal cohort (optional) | E-MTAB-14507 | — | exploratory comparison only; not required to reproduce any reported result |

### ⚠️ Case definitions — important when interpreting the analyses

* **GSE32472** was graded by the original investigators with the consensus definition in force at
  recruitment (2001 NICHD/NHLBI): moderate and severe BPD defined by the level of oxygen requirement
  and the need for positive-pressure support at 36 weeks' PMA. The manuscript reports this, and adds a
  re-grading analysis based on the oxygen status at 36 weeks' PMA recorded in the public series matrix
  (2018 NICHD oxygen criterion; Table S9).
* **GSE220135** assigns BPD status **binarily, without grading at 36 weeks' PMA**. Mild cases
  therefore enter the case group in that cohort. This is the main reason the external discrimination
  is lower than the training-set estimate, and it is the subject of the case-definition sensitivity
  analyses in `07_cohort_validation/auc_broad_grouping.py` and `external_delta_auc.py`.

## 2. Third-party reference resources

| Resource | Version | Where to get it | Used for |
|---|---|---|---|
| Inflammation-related gene (IRG) list, 3,336 genes | GeneCards v5.22 + MSigDB v2024.1.Hs | genecards.org, msigdb.org — query `"inflammatory response"` | defining the IRG universe (step 1) |
| CIBERSORT LM22 signature matrix | CIBERSORT 1.06 | cibersortx.stanford.edu | immune cell deconvolution (Fig 6) |
| PBMC single-cell reference | 10x Genomics v3, GRCh38-3.0.0 | 10xgenomics.com | cellular localisation; 10,667 cells retained after QC (Fig 10) |
| Human Protein Atlas immune-cell RNA-seq | — | proteinatlas.org, file `rna_immune_cell.tsv` | orthogonal validation of cellular origin (Fig S4) |

## 3. Directory layout expected by the scripts

Paths inside the scripts are written as `P("<relative path>")`, resolved against the repository root
Place the data as follows:

```
data/
├── raw/
│   ├── GSE32472/
│   │   ├── GSE32472_series_matrix.txt            # series matrix (download from GEO)
│   │   ├── 正常组与中重度BPD的比较/                 # MS-BPD vs control comparisons
│   │   │   ├── 12-17d/GSE32472_12-17d.csv        # day-14 expression matrix (63 samples)
│   │   │   ├── 22-27d/GSE32472_22-27d.xlsx       # day-28 window
│   │   │   └── 12-17d/炎症并集/Nomogram/预测值-系数.xlsx
│   │   ├── 正常组与BPD比较/
│   │   │   └── 12-17d/GSE32472_12-17d.csv        # day-14, all 97 samples
│   │   └── GSE25286_mouse_validation/complete/
│   │       ├── expression_matrix.csv              # original filename: 表达谱.csv
│   │       └── sample_metadata.xlsx               # original filename: 样本信息.xlsx
│   ├── GSE220135/
│   │   ├── GSE220135_14d.csv
│   │   └── GSE220135_28d.csv
│   ├── E-MTAB-14507/                              # optional
│   ├── GSE32472.xlsx                              # severity grades at 28 days
│   └── supporting_information/S9 Data.xlsx        # external cohort clinical annotations
└── derived/
    ├── training_set_clinical_genotype.csv         # 63 x 9, 24 MS-BPD / 39 control
    ├── gse220135_samples.csv
    ├── Clinical_subset_AUC_train.csv
    ├── Clinical_subset_AUC_external.csv
    ├── NRI_IDI_summary.json
    └── NRI_IDI_bootstrap_detail.csv
```

> **Note on Chinese directory names.** The intermediate directories inside `data/raw/GSE32472/` keep
> the names used in the original analysis (`正常组与中重度BPD的比较`, `炎症并集`, …) so that the scripts
> match the archived layout exactly. If you prefer ASCII names, rename them consistently and update
> the corresponding `P("...")` strings.

## 4. Input files that are derived rather than downloaded

A few inputs are intermediate tables produced by earlier steps of the workflow rather than public raw
data — for example `training_set_clinical_genotype.csv` (the 63-infant training set with gestational
age, birth weight, sex and outcome) and `gse220135_samples.csv`. These are small (≲ 20 KB), contain no
identifiable information, and are listed in section 3. They can be rebuilt from the public series
matrices with `01_data_prep/01_geo_download_preprocess_GSE32472.R` and
`07_cohort_validation/rnaseq_deseq2.py`.

## 5. What is deliberately not included

* **Identifying or clinical information.** Only the aggregate variables already published in the
  manuscript (gestational age, birth weight, sex, BPD grade) appear in any input table.
* **The manuscript files.** The `.docx` manuscripts, reviewer responses and submission packages are
  not part of this repository.
* **Figure files.** All figures written by the scripts go to `results/` and `figures/`, which are
  git-ignored. The published figures are available with the article.

## 6. Ethics

Both human datasets are fully de-identified public data, and the mouse and rat datasets are
re-analyses of publicly deposited data; no new animal procedures were performed. The re-analysis was
granted an exemption from ethical review by the Medical Ethics Committee of the authors'
institution (see the manuscript's Ethics statement).
