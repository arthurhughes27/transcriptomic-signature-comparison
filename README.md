# transcriptomic-signature-comparison

## Overview

This repository implements a reproducible R pipeline for a **multi-study, multi-vaccine analysis of the HIPC ImmuneSignatures2 (IS2) dataset**. It has two related goals:

1. **Comparative transcriptomic biomarker discovery** — identify and compare gene-set-level ("modular") transcriptomic signatures of vaccination across 13 vaccines (covering pathogens such as Influenza, Yellow Fever, Ebola, Tuberculosis, Meningococcus, Pneumococcus, Malaria, HIV, Hepatitis A/B, Smallpox, and Varicella Zoster) and ~30 underlying studies, at multiple timepoints post-vaccination (Days 1, 3, 7, etc.).
2. **Robustness analysis** — evaluate how sensitive these differential gene-set signatures are to methodological choices (e.g. differential gene-set analysis (DGSA) method, multiple-testing correction approach/scope, significance thresholds), by running the same comparisons with two independent DGSA methods — `dearseq` and `QuSAGE` — and comparing results.

The analysis works at the level of curated **gene sets** (Blood Transcriptional Modules, "BTM", and BloodGen3 Modules, "BG3M") rather than individual genes, and relates transcriptomic activity to post-vaccination **antibody/immune response** (maximum fold-change, MFC) as well as clinical/demographic covariates.

**Note on data:** raw and processed data (`data-raw/`, `data/`) and generated results (`output/`) are gitignored and are not stored in this repository for privacy reasons. The repository only contains the analysis code; all scripts assume these folders exist locally, populated by the user.

## Repository structure

```
R/                                          # Modular, reusable analysis functions (sourced, not a package)
├── load_all.R                              # source(fs::path("R", "load_all.R")) loads every function below
├── dgsa_common.R                           # Shared DGSA helpers: comparison listing, sample pairing, covariate matrices, gene-set correlation
├── dgsa_dearseq.R                          # Modular dearseq DGSA: run_dearseq_comparison() for one vaccine x timepoint comparison
├── dgsa_qusage.R                           # Modular QuSAGE DGSA: run_qusage_comparison(), same contract as run_dearseq_comparison()
├── specifications.R                        # Specification grid generators (Table 2.1): raw (66) and post-hoc (540) tiers
├── postprocessing.R                        # Tidies a raw DGSA results list and applies p-value adjustment (3 scopes x 6 methods)
├── robustness_metrics.R                    # Computes the robustness metric pi_{g,v,j} by accumulating counts across raw runs
└── synthetic_data.R                        # make_synthetic_is2_data(): small synthetic dataset for testing without real data

tests/                                      # testthat suite for R/ (run via `Rscript tests/run_tests.R`)
└── testthat/

analysis/
├── analysis_master.R                     # Top-level entry point: runs preprocessing (then reanalysis/descriptive)
├── preprocessing/
│   ├── preprocessing_master.R            # Runs all preprocessing steps in order
│   ├── preprocessing_clinical.R          # Clinical/demographic data cleaning + vaccine/study metadata
│   ├── preprocessing_expression.R        # Gene expression matrices (normalised/non-normalised)
│   ├── preprocessing_immuneresponse.R    # Derives maximum fold-change (MFC) immune response outcomes
│   ├── preprocessing_BTM.R               # Blood Transcriptional Modules gene sets
│   ├── preprocessing_BG3M.R              # BloodGen3 Modules gene sets
│   └── preprocessing_merging.R           # Merges clinical + immune response + expression into analysis-ready tables
├── descriptive/
│   ├── descriptive_master.R              # Runs descriptive analyses
│   ├── is2_bubble_plot.R                 # Main-text figure: study x timepoint sample bubble plot
│   └── is2_appendix_descriptives.R       # Appendix: covariate distributions + study-level sample size table
└── reanalysis/
    ├── reanalysis_master.R               # Runs the dearseq DGSA pipeline end-to-end
    ├── dearseq_dgsa.R                    # Thin driver: baseline dearseq specification, via R/dgsa_common.R + R/dgsa_dearseq.R
    ├── process_dearseq_dgsa_results.R    # Tidies dearseq results, assigns colours, applies 6 p-value correction methods x 3 scopes
    ├── qusage_dgsa.R                     # Parallel DGSA/meta-analysis pipeline using QuSAGE (based on Hagan et al. 2022)
    ├── dgsa_comparison_example.R         # Example comparison of dearseq vs QuSAGE under fixed hyperparameters
    └── plot_circos.R                     # Circos plots comparing DGSA results across methods/timepoints/directions

manuscript/figures/                       # Pipeline overview and robustness diagrams for the write-up
transcriptomic-signature-comparison.Rproj
```

Analysis outputs are written to `data/` (processed intermediate tables) and `output/results/` and `output/figures/` (final results and plots), organised by pipeline stage (`descriptive/`, `reanalysis/`, `qusage/`).

## Data pipeline

### 1. Preprocessing (`analysis/preprocessing/`)

Raw inputs live in `data-raw/` (not tracked in git) and include:
- `all_noNorm_eset.rds`, `all_norm_eset.rds`, `young_noNorm_eset.rds`, `young_norm_eset.rds`, `old_noNorm_eset.rds`, `young/old_noNorm_withResponse_eset.rds` — `Biobase::ExpressionSet` objects with expression matrices + phenotype/clinical data, split by age cohort and normalisation status.
- `elisa_*.xlsx`, `elispot_*.xlsx`, `neut_ab_titer_*.xlsx`, `hai_*.xlsx` — raw immune-response assay exports.
- `BTM_for_GSEA_20131008.gmt`, `BTM_functional_groups.txt` — Blood Transcriptional Module gene sets and their functional aggregate labels.
- `Suppl_File_1_BIOINF.xls` — BloodGen3 Module (BG3M) gene set definitions.

Processed outputs (in `data/`), one row per sample/participant unless noted:

| File | Grain | Key columns |
|---|---|---|
| `hipc_clinical.rds` | 1 row/sample | `participant_id`, `study_accession`, `study_accession_unique` (disambiguates studies spanning >1 vaccine, e.g. `SDY1260a/b`), `gender`, `race`, `ethnicity`, `age_imputed`, `pathogen`, `vaccine_type` (abbreviated: CJ, IN, IN/RP, LV, PS, RVV, RP), `vaccine_name`, `vaccine_name_short`, `vaccine_colour`, `study_colour`, `study_time_collected`, `time_post_last_vax` |
| `all_norm_expr.rds`, `young_noNorm_expr.rds`, `young_norm_expr.rds` | 1 row/sample | `participant_id`, `study_time_collected`, then one column per gene (lowercase HGNC symbol, e.g. `a1cf` … `zzz3`) |
| `hipc_immResp.rds` | 1 row/participant | `participant_id`, then per-assay blocks prefixed `immResp_MFC_{anyAssay,nAb,elisa,elispot,hai}_*` with `_response_strain_analyte`, `_pre_time`, `_post_time`, `_pre_value`, `_post_value`, `_MFC`, `_log2_MFC` |
| `BTM_processed.rds`, `BG3M_processed.rds` | gene-set object | `genesets` (list of lowercase gene symbol vectors), `geneset.names`, `geneset.descriptions`, `geneset.names.descriptions`, `geneset.aggregates` (functional category, factor for BTM) |
| `hipc_merged_all_norm.rds`, `hipc_merged_young_noNorm.rds`, `hipc_merged_young_norm.rds` | 1 row/sample | full outer/right join of clinical + immune response + expression, keyed by `participant_id` + `study_time_collected` |

Key preprocessing decisions: "Unknown"/"Not Specified" gender & race are collapsed; MFC is computed as `post/pre` fold-change (log2-transformed) using the most recent pre-vaccination sample as baseline, restricted to a Day 21–35 (and ≤0) response window; duplicate baseline samples per participant are deduplicated and their timepoint recoded to 0.

### 2. Descriptive analysis (`analysis/descriptive/`)

Generates the dataset-description material for the thesis chapter, from `hipc_merged_all_norm.rds`, coloured consistently using the `vaccine_colour`/`study_colour` palettes defined during preprocessing:

- **`is2_bubble_plot.R`** — the main-text figure: a bubble plot with studies on the y-axis (coloured by vaccine), days post-vaccination on the x-axis, and bubble size proportional to the number of transcriptomic samples available per study x timepoint.
- **`is2_appendix_descriptives.R`** — Appendix A material: per-study covariate distributions (age, gender, race, immune-response assay availability and pre/post-vaccination distributions for nAb/HAI/ELISA) and a study-level sample-size summary table (participants, samples, timepoints sampled), saved to `output/tables/descriptive/` as both CSV and a `\input{}`-ready LaTeX table.

### 3. Reanalysis / DGSA (`analysis/reanalysis/`, `R/`)

Core biomarker-discovery and robustness analysis. Two independent differential gene-set analysis (DGSA) methods are run over the **same** vaccine x timepoint comparisons (post-vaccination vs. pre-vaccination self-baseline, adjusting for age, sex, and study).

To keep this tractable across many analytical specifications (Chapter 2's robustness/specification analysis), the actual DGSA logic lives in modular, reusable functions under `R/` rather than in the `analysis/` driver scripts:

- `R/dgsa_common.R` — shared building blocks used by every DGSA method: `list_valid_comparisons()` (which vaccine x timepoint pairs have usable data), `filter_paired_samples()` (pre-/post-vaccination sample pairing per participant), `build_covariate_matrix()` (design matrix for an arbitrary covariate subset, with automatic collinearity removal), and `calculate_gs_correlation()` (gene-set-level correlation with immune response).
- `R/dgsa_dearseq.R` — `run_dearseq_comparison(vax, day, hipc, BTM, gene_names, covariates, which_weights, gene_based_weights)`: runs one comparison with `dearseq::dgsa_seq()`, parameterised by the dearseq-specific hyperparameters from Table 2.1 (covariate set; mean-variance weighting method/level), plus gene-set scoring (`calculate_scores()`).
- `R/dgsa_qusage.R` — `run_qusage_comparison(vax, day, hipc, BTM, gene_names, equal_variance, sample_scope)`: the QuSAGE equivalent, with the same contract (same inputs/outputs) as `run_dearseq_comparison()`. Runs `qusage::qusage()` independently per contributing study and meta-analyses across studies via `qusage::combinePDFs()` (Meng et al. 2019), parameterised by the QuSAGE-specific hyperparameter from Table 2.1 (equal-variance assumption). Replaces the original, monolithic `qusage_dgsa.R` (adapted from Hagan et al. 2022) — see the design notes at the top of the file for the (deliberate, documented) differences from that script, including why `sample_scope = "paired"` (the default, used by the main pipeline) gives systematically different — not wrong, but genuinely different — p-values/scores than the original script, and how `sample_scope = "study"` reproduces the original script's sample selection exactly for validation purposes.
- `R/specifications.R` — `build_raw_specification_grid()` / `build_posthoc_specification_grid()` / `build_full_specification_grid()`: builds the Table 2.1 specification space as two independent tiers (66 "raw" specifications requiring an actual DGSA run; 540 cheap "post-hoc" p-value-adjustment/threshold combinations applied afterwards), so the full 35,640-specification grid can be evaluated without running DGSA 35,640 times.
- `R/postprocessing.R` — `build_tidy_dgsa_results()`: tidies a raw results list (from either `run_dearseq_comparison()` or `run_qusage_comparison()`) into one long-format table and applies p-value adjustment across all 3 scopes x 6 methods.
- `R/robustness_metrics.R` — `count_significant_specifications()` / `accumulate_robustness_counts()` / `compute_robustness_metric()`: computes the robustness metric pi_{g,v,j} (Section 2.2.4) by streaming-accumulating significance counts across raw runs, without materialising the full specification x comparison table.
- `R/synthetic_data.R` — `make_synthetic_is2_data()`, a small synthetic dataset generator matching the `hipc_merged_*.rds`/`BTM_processed.rds` schema, used by the `tests/` suite (and available for local smoke-testing) without needing access to the real, private IS2 data.

`analysis/reanalysis/dearseq_dgsa.R` is a thin driver: it loads the data, sources `R/load_all.R`, and loops `run_dearseq_comparison()` over every valid comparison at the baseline specification (bolded options in Table 2.1), checkpointing results incrementally. The Stage 1 (baseline method comparison) driver, and the Stage 2 specification-analysis driver scripts, are planned as later additions.

- **`dearseq_dgsa.R`** — thin driver for the baseline dearseq specification; see `R/dgsa_dearseq.R` above.
- **`process_dearseq_dgsa_results.R`** — reshapes the raw list into one tidy dataframe (`output/results/reanalysis/dearseq_dgsa_results_processed.rds`), one row per gene set x comparison, and appends p-value corrections for **6 methods** (`holm`, `hochberg`, `hommel`, `bonferroni`, `BH`, `BY`) applied at **3 scopes** (`global`, `withinTime`, `withinComparison`) — this multiplicity is itself part of the robustness assessment. Not yet rewired onto `R/postprocessing.R` (planned alongside the Stage 1 driver).
- **`qusage_dgsa.R`** — the original, monolithic pipeline adapted from Hagan et al. (2022); superseded by `R/dgsa_qusage.R` and pending deletion once its replacement is validated (see `validate_qusage_refactor.R`, a temporary script comparing the two implementations' output for one comparison — also pending deletion).
- **`dgsa_comparison_example.R`** — worked example comparing dearseq and QuSAGE outputs under one fixed set of hyperparameters.
- **`plot_circos.R`** — `circlize`-based circos plots comparing the two methods' significant gene sets (by aggregate/functional category, direction of regulation, and correlation with immune response) across vaccines and timepoints, with configurable p-value correction method/scope, significance threshold, and score-filtering options — used to visualise where the two DGSA methods agree or diverge (the robustness question).

### Combined results dataframe schema

Both `dearseq_dgsa_results_processed.rds` and the QuSAGE equivalent share a common schema so they can be row-bound for comparison: `comparison`, `condition` (vaccine), `time` (day post-vax), `condition.colour`, `gs.name`, `gs.description`, `gs.name.description`, `gs.aggregate`, `gs.colour`, `activation.score`, `fc.score`, `mean.corr`, `corr.mean`, `rawPval`, `method` (`"dearseq"` / `"qusage"`), plus `{scope}.adjPval_{correction}` columns for every scope x correction combination.

## Running the pipeline

```r
# from the project root, with data-raw/ populated
source("analysis/analysis_master.R")
```

This currently runs preprocessing end-to-end (`preprocessing_master.R`); the reanalysis (`reanalysis_master.R`, which runs `dearseq_dgsa.R` → `process_dearseq_dgsa_results.R` → `qusage_dgsa.R`) and descriptive (`descriptive_master.R`) scripts can be sourced separately once preprocessing has produced the files in `data/`. Some reanalysis scripts (e.g. dearseq DGSA, QuSAGE) are long-running and checkpoint their results to `output/results/` so interrupted runs can resume rather than restart.
