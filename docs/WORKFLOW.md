# Analysis workflow, decision rules and outputs

Companion to the code-to-figure map in the [README](../README.md). This document records the
**thresholds and design decisions** that are not evident from the script names alone, so that a
reader can tell what "reproducing a result" means for each step.

---

## Stage 1 — Discovery: preprocessing to the four-gene signature

| Step | Script | Decision rule / parameter |
|---|---|---|
| Preprocessing | `01_data_prep/01_geo_download_preprocess_GSE32472.R` | GPL6244 probes annotated to gene symbols; duplicate genes collapsed by highest mean expression; `log2` transform; limma background correction + quantile normalisation. Split by day of sampling into **day-14 = training set (24 MS-BPD / 39 controls)** and the day-28 temporal window (24 / 38). |
| IRG universe | same script | GeneCards v5.22 ∪ MSigDB v2024.1.Hs, query `"inflammatory response"` → **3,336 IRGs**; intersected with genes detected in the training set → **1,498**. |
| Differential expression | `01_data_prep/02_deg_limma_volcano_heatmap_venn.R` | limma, **adj.P < 0.05 and \|log2FC\| > 0.58** (FC > 1.5) → 139 DEGs (84 up, 55 down). DEG ∩ IRG → **57 DE-IRGs**. |
| GO / KEGG | `02_enrichment/03_go_kegg_circle_plot.R` | clusterProfiler ORA against the 1,498 detected IRGs as background; BH adjustment; the circle plot colours each DE-IRG by its logFC. |
| ORA (C2.CP) | `02_enrichment/04_ora_msigdb_C2CP.R` | `clusterProfiler::enricher`, MSigDB C2.CP v2025.1, Fisher's exact test + BH. |
| Single-gene GSEA | `02_enrichment/05_single_gene_GSEA_KEGG.R` | Median split of each hub gene → limma → rank genes by log2FC → pre-ranked GSEA (KEGG). **Reported if \|NES\| > 1.5 and FDR < 0.25** (exploratory threshold, stated in the manuscript). |
| GSVA | `02_enrichment/06_GSVA_hallmark.R` | Hallmark gene sets; limma comparison of high- vs low-expression groups; **\|t\| > 2 and BH-adjusted p < 0.05**. |
| PPI network | `03_network/07_ppi_string_cytohubba_workflow.md` | STRING v12.0, confidence ≥ 0.4, 48 of the 57 DE-IRGs with PPI evidence; CytoHubba **MCC**, **EPC** and **Degree** rankings (Top 30 each), intersection of the three lists. |
| Regulatory networks | `03_network/08_mirna_tf_drug_network_workflow.md` | miRTarBase (miRNA), JASPAR (TF), DGIdb 5.0 (drug); networks drawn and exported in Cytoscape 3.10.3. |
| Feature selection | `04_model/09_lasso_svmrfe_roc.R` | LASSO (10-fold CV for lambda) ∩ SVM-RFE → **4 consensus genes: CD3G, PRTN3, CD6, THBS1**. |
| Nomogram | `04_model/10_nomogram_calibration_dca.R` | `rms` nomogram, calibration curve, DCA. **DCA is exploratory only** — net benefit depends on outcome prevalence and this is a case-control design. |
| Nested CV | `04_model/26_nested_cv_R.R` | The whole selection pipeline (differential expression → LASSO → SVM-RFE) is repeated **inside each training fold**. V0 = 0.873 (reproduces the reported LOOCV AUC; correctness check), V1 = 0.870 (feature pool fixed at the 57 DE-IRGs), **V2 = 0.794** (everything re-done inside the fold). See [ENVIRONMENT.md](ENVIRONMENT.md) for the formula bug this script had to avoid. |
| Immune infiltration | `05_immune_infiltration/11_cibersort_immune_landscape.R` | CIBERSORT 1.06 with LM22; samples with deconvolution `P < 0.05` retained; group comparison by Mann-Whitney U with BH adjustment; Spearman correlations for cell-cell and cell-gene matrices. |
| In vitro | `06_validation/12_qpcr_ddct_analysis.R` | A549 hyperoxia vs normoxia; **Ct > 35 = below detection**; 2^-ΔΔCt with β-actin as reference and the normoxia group as calibrator; unpaired t test. |

## Stage 2 — Validation and sensitivity analyses

### The coefficients are locked

A single set of coefficients is estimated **once**, in the day-14 training set, and then applied
without re-estimation everywhere:

```
logit(P) = −3.375 + 0.610·PRTN3 − 0.820·CD3G − 1.402·CD6 + 1.560·THBS1
```

Every transport, matching, restriction and case-definition analysis holds these fixed. This is what
makes the comparisons across cohorts interpretable — and it is why the day-28 and external estimates
are out-of-sample even though they re-use the same infants.

| Analysis | Script | Design decision |
|---|---|---|
| Day-28 temporal window | `07_cohort_validation/external_validation_gse32472.py` | Same infants as the training set, later time point → **non-independent**, but out-of-sample for the coefficients. |
| External cohort | `07_cohort_validation/rnaseq_deseq2.py` | GSE220135 processed with DESeq2 (vst); coefficients applied unchanged. |
| ROC / thresholds | `07_cohort_validation/roc_diagnosis_4v3.py` | **Youden-optimal** cut-off. Sensitivity and specificity therefore differ between case definitions — the threshold is re-derived for each and reported per row in Table S3. |
| Model comparison | `07_cohort_validation/external_clinical_comparison.py` | All models (signature / perinatal-only / combined) are fitted **in the training set** and applied to the external cohort, so the comparison reflects transportability rather than refitting. An internal-refit "apparent" AUC is reported alongside as a reference. |
| Bootstrap internal validation | `07_cohort_validation/bootstrap_auc_internal_validation.py` | 1000 bootstrap replicates, **Harrell optimism correction**, plus a 0.632+ estimate as a supplementary check. |
| GA / BW matching | `07_cohort_validation/ga_match_strat.py` | Nearest-neighbour matching with calipers ±1 week, ±2 weeks and ±150 g; additionally restricted to GA ≥ 26 / ≥ 27 / ≥ 28 weeks. |
| Case-definition sensitivity | `07_cohort_validation/auc_broad_grouping.py`, `external_delta_auc.py` | Case group broadened to **all grades of BPD**, with the control group held constant, to quantify the attenuation seen in GSE220135 (which does not grade severity). |
| 2018 re-grading | manuscript Table S9 | Uses the `oxygen 36 weeks postmenstrual age` field of the GSE32472 series matrix. Because the earlier criteria already required oxygen at 36 weeks for the moderate/severe categories, the revision reclassifies the mild group only — the case group is unchanged. |
| Incremental value | `07_cohort_validation/baseline_adjust.py`, `gene_vs_clinical_full.py`, `gen_nri_idi_data.py` | Multivariable adjustment (GA, BW, sex), ΔAUC with DeLong test, and **continuous NRI / IDI with 2000 bootstrap replicates** (`lib/delong_nri_idi.py`). |
| Day-14 / day-28 increment | `07_cohort_validation/day28_increment_check.py`, `day28_increment_repro.py` | Recomputes the increment in both windows with fixed coefficients. |

### Interpreting the increment results

The manuscript reports that the signature does **not** outperform a perinatal-variables-only model on
AUC, that the AUC increment is marginal, and that the reclassification-based increment (NRI, IDI) is
significant in the training set but **not evident after birth-weight matching or inverse probability
weighting**. All of those numbers come from the scripts in this stage. The matching and IPTW
estimates in particular rest on small subsets (11 matched pairs; 63 infants for IPTW) and are
reported as imprecise rather than as evidence of absence.

## Stage 3 — Cross-species, single-cell and orthogonal validation

| Analysis | Script | Design decision |
|---|---|---|
| Mouse model | `06_validation/gse25286_limma_recalc.R` | limma re-analysis of GSE25286; BH-adjusted P values. |
| Rat model | (GSE212098 inputs) | Independent hyperoxia model, postnatal days 1-14; directional concordance across time points. |
| Single-cell reference | `08_singlecell_hpa/01_load_qc_integrate.py` | 10x PBMC v3 reference; QC **300-6,000 genes per cell and mitochondrial fraction ≤ 15%** → 10,667 cells retained; harmonypy integration, Leiden clustering. |
| Annotation | `02_annotate.py`, `02b_marker_table.py` | Ten canonical PBMC populations verified against eleven marker gene sets (Fig S3). |
| Localisation | `10_ref_pbmc_localize.py`, `11_localize_stats_and_figs.py` | Per-cell-type detection proportions and mean log-normalised expression for the four hub genes. ⚠️ **The atlas is a single healthy adult donor**, so cells were treated as independent observations; the enrichment confidence intervals are therefore narrower than donor-level replication would give. This caveat is stated in the manuscript. |
| HPA orthogonal validation | `13_hpa_crossvalidation_figure.py` | Human Protein Atlas immune-cell bulk RNA-seq (nTPM), 19 purified populations. |
| Doublet checks | `42_scrublet_doublet_check.py`, `44_doublet_impact.py` | Scrublet at expected rates 0.06 and 0.075 (both 527 cells, 4.94%); impact of removing predicted doublets on cell-type counts. |
| Panel assembly | `12_composite_figure.py` | Assembles Fig 10 from the per-panel outputs written by the localisation and HPA scripts. |

## Reproducing a specific reported number

1. Find the step in the README's code-to-figure map.
2. Check its **status** — if `RECONSTRUCTED`, the script has not been validated in the original
   environment and should be checked against the archived outputs first.
3. Read the threshold/parameter row above; most disagreements come from a different cut-off or a
   different case definition rather than from the model itself.
4. Seeds: see section 5 of the README.
