# Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for moderate-to-severe bronchopulmonary dysplasia

Analysis code accompanying the manuscript:

> **A Four-Gene Inflammatory Immune-Imbalance Signature (CD3G/PRTN3/CD6/THBS1) in Moderate-to-Severe Bronchopulmonary Dysplasia: A Machine-Learning and Cross-Platform Validation Study**

---

## 1. Contents

```
├── 01_data_prep/              GSE32472 acquisition, normalisation, limma DEGs
├── 02_enrichment/             GO/KEGG, ORA (MSigDB C2.CP), single-gene GSEA, GSVA (Hallmark)
├── 03_network/                PPI (STRING + CytoHubba) and miRNA/TF/drug networks — workflow documents
├── 04_model/                  LASSO + SVM-RFE, nomogram, calibration, DCA, nested cross-validation
├── 05_immune_infiltration/    CIBERSORT (LM22) immune landscape
├── 06_validation/             A549 qRT-PCR (2^-ΔΔCt), mouse GSE25286 limma re-analysis
├── 07_cohort_validation/      external cohort (GSE220135), discrimination, matching, IPTW,
│                              NRI/IDI, case-definition and birth-weight sensitivity analyses
├── 08_singlecell_hpa/         PBMC single-cell localisation (Scanpy), doublet checks, HPA validation
├── 09_plotting/               supporting plots
├── 10_rebuttal_reproducibility/  scripts and inputs behind the numbers quoted in the
│                              response letter (Data S17 and S18) — see its own README
├── lib/                       shared statistics module (DeLong, NRI, IDI, VIF, all-subset AUC)
├── docs/                      data sources, environment and workflow
├── data/                      input data — not included, see docs/DATA.md
├── results/                   script outputs (created on first run)
└── figures/                   figure files (created on first run)
```

**53 scripts** in total: 20 R, 31 Python, 2 workflow documents, plus the 14 input/output
tables under `10_rebuttal_reproducibility/data/`. The 40 scripts in `01`–`09` reproduce the
manuscript; the 13 scripts in `10_rebuttal_reproducibility/` reproduce the numbers quoted in
the response letter and tabulated as Supplementary Data S17 and S18.

**Not included, by design**

* The 20 manuscript-editing scripts used during revision. They operate on the `.docx` manuscript
  (redlining, reference renumbering, font unification, package auditing) and are **not needed to
  reproduce any result**.
* The regulatory-network figures (Figure 9A-C). They were built interactively in Cytoscape, so no
  script exists; the exact parameters, database queries and export settings are documented in
  `03_network/`.

## 2. Quick start

```bash
git clone <REPOSITORY-URL>
cd <repository>

# Point the scripts at the repository root
export BPD_REPO="$(pwd)"          # Linux / macOS
set BPD_REPO=%CD%                 # Windows cmd
$env:BPD_REPO = (Get-Location)    # Windows PowerShell
```

Place the input files as described in **[docs/DATA.md](docs/DATA.md)**, then run a script, e.g.

```bash
Rscript 04_model/26_nested_cv_R.R
python 07_cohort_validation/roc_diagnosis_4v3.py
```

If `BPD_REPO` is not set, each script falls back to its own parent directory, so running from the
repository root also works.

## 3. Script status — please read

Every script carries a `Status:` line in its header. Three categories:

| Status | Count | Meaning |
|---|---|---|
| **ORIGINAL** | 28 | The script as used in the analysis, with paths made portable. |
| **RECONSTRUCTED** | 10 | Rewritten from the *written description* in the Methods section, for code release and reviewer reproducibility. **These were not validated in the original authors' analysis environment.** Each one still needs to be run locally against the archived result files before being relied upon — in particular the AUC, coefficient and cell-fraction values. |
| **GUI WORKFLOW** | 2 | Cytoscape-based network construction; no script is possible. Parameters, queries and export settings are documented instead. |

> ⚠️ The reconstructed scripts cover steps 1-6 and 9-12 of the map below (data pre-processing
> through CIBERSORT and qRT-PCR). All cohort-level validation, sensitivity, single-cell and HPA
> analyses (steps 13-26) are original, run-verified code.

## 4. Code-to-figure map

| # | Analysis step | Script | Output | Status |
|---|---|---|---|---|
| 1 | GSE32472 download and preprocessing | `01_data_prep/01_geo_download_preprocess_GSE32472.R` | Data S1, S2 | reconstructed |
| 2 | limma DEGs + volcano / heatmap / Venn | `01_data_prep/02_deg_limma_volcano_heatmap_venn.R` | Fig 2A-C; Data S4 | reconstructed |
| 3 | GO / KEGG + logFC circle plot | `02_enrichment/03_go_kegg_circle_plot.R` | Fig 3A-B | reconstructed |
| 4 | ORA (MSigDB C2.CP) | `02_enrichment/04_ora_msigdb_C2CP.R` | Fig 3C | reconstructed |
| 5 | Single-gene GSEA (KEGG, pre-ranked) | `02_enrichment/05_single_gene_GSEA_KEGG.R` | Fig 7 | reconstructed |
| 6 | GSVA (Hallmark) | `02_enrichment/06_GSVA_hallmark.R` | Fig 8 | reconstructed |
| 7 | PPI (STRING 12.0) + CytoHubba intersection | `03_network/07_ppi_string_cytohubba_workflow.md` | Fig 4A-E | GUI workflow |
| 8 | miRNA / TF / drug regulatory networks | `03_network/08_mirna_tf_drug_network_workflow.md` | Fig 9A-C | GUI workflow |
| 9 | LASSO + SVM-RFE + ROC | `04_model/09_lasso_svmrfe_roc.R` | Fig 4F-H; Fig 5D-E | reconstructed |
| 10 | Nomogram + calibration + DCA | `04_model/10_nomogram_calibration_dca.R` | Fig 5A-C | reconstructed |
| 11 | CIBERSORT immune infiltration | `05_immune_infiltration/11_cibersort_immune_landscape.R` | Fig 6A-D | reconstructed |
| 12 | A549 qRT-PCR (2^-ΔΔCt) | `06_validation/12_qpcr_ddct_analysis.R` | Fig 5J | reconstructed |
| 13 | Mouse GSE25286 limma re-analysis | `06_validation/gse25286_limma_recalc.R` | Table S5 | original |
| 14 | Day-28 temporal validation (coefficients locked) | `07_cohort_validation/external_validation_gse32472.py` | Fig 5F-G | original |
| 15 | External cohort GSE220135 (DESeq2 vst) | `07_cohort_validation/rnaseq_deseq2.py` | Table S2 | original |
| 16 | ROC / Youden threshold | `07_cohort_validation/roc_diagnosis_4v3.py` | Fig 5D-E; Table S3 | original |
| 17 | External ΔAUC and model comparison | `07_cohort_validation/external_delta_auc.py`, `external_clinical_comparison.py` | Table S2 | original |
| 18 | 1000× bootstrap internal validation | `07_cohort_validation/bootstrap_auc_internal_validation.py`, `plot_bootstrap*.py` | Fig S1, S5; Data S6 | original |
| 19 | Gestational age / birth weight matching and stratification | `07_cohort_validation/ga_match_strat.py`, `auc_broad_grouping.py` | Table S3 | original |
| 20 | Clinical adjustment and incremental value | `07_cohort_validation/baseline_adjust.py`, `gene_vs_clinical_full.py` | Table S3, S4 | original |
| 21 | NRI / IDI (2000× bootstrap) | `07_cohort_validation/gen_nri_idi_data.py` (`lib/delong_nri_idi.py`) | Table S4 | original |
| 22 | Day-14 / day-28 increment re-computation | `07_cohort_validation/day28_increment_check.py`, `day28_increment_repro.py` | Table S4 | original |
| 23 | Single-cell localisation (Scanpy) | `08_singlecell_hpa/01_load_qc_integrate.py`, `02_annotate.py`, `02b_marker_table.py`, `10_ref_pbmc_localize.py`, `11_localize_stats_and_figs.py` | Fig 10 | original |
| 24 | HPA orthogonal validation | `08_singlecell_hpa/13_hpa_crossvalidation_figure.py` | Fig S4, Data S8 | original |
| 25 | Figure 10 panel assembly | `08_singlecell_hpa/12_composite_figure.py` | Fig 10 | original |
| 26 | Nested cross-validation (selection inside the fold) | `04_model/26_nested_cv_R.R` | Data S15 | original, run-verified |

A longer description of each step, including the thresholds and design decisions that are not
evident from the script names, is in **[docs/WORKFLOW.md](docs/WORKFLOW.md)**.

Four further scripts belong to the workflow but do not map to a single figure:

| Script | Purpose |
|---|---|
| `07_cohort_validation/merge_14d_28d_deg.py` | Merges the day-14 and day-28 GSE220135 differential-expression tables into one workbook. |
| `08_singlecell_hpa/42_scrublet_doublet_check.py` | Scrublet doublet detection at the expected rates used for the single-cell atlas (QC for Fig 10). |
| `08_singlecell_hpa/44_doublet_impact.py` | Effect of removing predicted doublets on the cell-type counts (QC for Fig 10). |
| `09_plotting/13_bootstrap_auc_curve.py` | Bootstrap AUC distribution curve (companion to step 18). |

The analyses added during revision — gene-set-definition sensitivity, SVM-RFE stability, the
GSE151974 single-cell localisation, the composition/expression decomposition and the
deconvolution-free marker-set check — sit in **`10_rebuttal_reproducibility/`**; its README maps
each script to the reviewer comment it answers and to the sheet of Data S17 or S18 it fills.

## 5. Reproducibility notes

* **Random seeds.** The original analysis fixed a single seed of `20250101`, but individual scripts as
  archived use their own seeds where a resampling step was added later (for example
  `np.random.default_rng(20260911)` in the NRI/IDI script, and the bootstrap replicates in
  `07_cohort_validation/`). Where a reported number depends on a seed, the seed is set explicitly at
  the top of the relevant script.
* **Fixed coefficients.** Throughout the external-validation and sensitivity analyses the four-gene
  coefficients are taken from the day-14 training set and are **never re-estimated**:
  `logit(P) = −3.375 + 0.610·PRTN3 − 0.820·CD3G − 1.402·CD6 + 1.560·THBS1`.
* **The animal datasets are re-analyses of publicly deposited data** (GSE25286, GSE212098); no new
  animal procedures were performed.
* **Windows command line and non-ASCII paths.** `Rscript` on Windows cannot open source files whose
  path contains non-ASCII characters (it fails with `cannot open the connection`). All R sources here
  are ASCII-only, comments included, as are all Python sources except for the Chinese directory, sheet
  and column names of the archived input files, which are kept verbatim in seven scripts under
  `07_cohort_validation/` (see [docs/DATA.md](docs/DATA.md)) — do not rename those inputs. If your
  clone path contains non-ASCII characters and you run R from the command line, either use RStudio or
  move the clone to an ASCII path (e.g. `C:/bpd-signature`).
* **Cytoscape.** The PPI and regulatory-network steps cannot be reproduced by script. Import the
  exported edge tables listed in `03_network/*.md` into Cytoscape 3.10.
* **Data are not included.** See [docs/DATA.md](docs/DATA.md).

## 6. Environment

Python 3.13 and R 4.4.3 (the original analysis environment used R 4.2.1 for the preprocessing steps).
Full package list and known pitfalls in **[docs/ENVIRONMENT.md](docs/ENVIRONMENT.md)**; Python
dependencies in [`requirements.txt`](requirements.txt).

## 7. Licence

MIT — see [LICENSE](LICENSE).

## 8. Citation

If you use this code, please cite the manuscript. The GEO datasets should be cited through
their original publications, listed in `docs/DATA.md`.

## 9. Contact

Correspondence regarding the code: please open an issue, or contact the corresponding author listed
in the manuscript.
