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
├── comparison_metrics.R                    # Stage 1: dearseq vs QuSAGE concordance metrics (correlation, agreement, kappa)
├── plot_helpers.R                          # Shared colour assignment + default condition/aggregate configuration
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
    ├── reanalysis_master.R               # Runs the full Stage 1 pipeline end-to-end
    ├── dearseq_dgsa.R                    # Thin driver: baseline dearseq specification, via R/dgsa_common.R + R/dgsa_dearseq.R
    ├── process_dearseq_dgsa_results.R    # Thin driver: tidy + adjust (R/postprocessing.R) + colour (R/plot_helpers.R)
    ├── qusage_dgsa.R                     # Thin driver: baseline QuSAGE specification, via R/dgsa_common.R + R/dgsa_qusage.R
    ├── process_qusage_dgsa_results.R     # Thin driver: tidy + adjust (R/postprocessing.R) + colour (R/plot_helpers.R)
    ├── dgsa_comparison_example.R         # Stage 1: dearseq vs QuSAGE baseline comparison (R/comparison_metrics.R)
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
- `R/dgsa_qusage.R` — `run_qusage_comparison(vax, day, hipc, BTM, gene_names, equal_variance, sample_scope)`: the QuSAGE equivalent, with the same contract (same inputs/outputs) as `run_dearseq_comparison()`. Runs `qusage::qusage()` independently per contributing study and meta-analyses across studies via `qusage::combinePDFs()` (Meng et al. 2019), parameterised by the QuSAGE-specific hyperparameter from Table 2.1 (equal-variance assumption). Replaces the original, monolithic `qusage_dgsa.R` (adapted from Hagan et al. 2022, since deleted): the goal here is a faithful, modular reproduction of Hagan et al.'s approach, not a methodologically "corrected" one, so `sample_scope` defaults to `"study"`, which reproduces the original script's (permissive) sample selection exactly — see the design notes at the top of the file. A stricter alternative (`sample_scope = "paired"`, matching dearseq's sample selection) is available but not used by default.
- `R/specifications.R` — `build_raw_specification_grid()` / `build_posthoc_specification_grid()` / `build_full_specification_grid()`: builds the Table 2.1 specification space as two independent tiers (66 "raw" specifications requiring an actual DGSA run; 540 cheap "post-hoc" p-value-adjustment/threshold combinations applied afterwards), so the full 35,640-specification grid can be evaluated without running DGSA 35,640 times.
- `R/postprocessing.R` — `build_tidy_dgsa_results()`: tidies a raw results list (from either `run_dearseq_comparison()` or `run_qusage_comparison()`) into one long-format table and applies p-value adjustment across all 3 scopes x 6 methods.
- `R/robustness_metrics.R` — `count_significant_specifications()` / `accumulate_robustness_counts()` / `compute_robustness_metric()`: computes the robustness metric pi_{g,v,j} (Section 2.2.4) by streaming-accumulating significance counts across raw runs, without materialising the full specification x comparison table.
- `R/comparison_metrics.R` — `compute_concordance_metrics()`: the Stage 1 (baseline method comparison) metrics — Spearman correlation of raw p-values and fold-change scores between dearseq and QuSAGE, plus agreement on significance calls (percent agreement, Cohen's kappa, both/either/neither-significant counts), both overall and per comparison.
- `R/plot_helpers.R` — `assign_dgsa_colours()` and the default vaccine-condition/gene-set-aggregate colour palettes, shared by the dearseq and QuSAGE result-processing scripts so their outputs use consistent colours and factor levels.
- `R/synthetic_data.R` — `make_synthetic_is2_data()`, a small synthetic dataset generator matching the `hipc_merged_*.rds`/`BTM_processed.rds` schema, used by the `tests/` suite (and available for local smoke-testing) without needing access to the real, private IS2 data.

The `analysis/reanalysis/` driver scripts are all thin: they load data, source `R/load_all.R`, and call into the modules above.

- **`dearseq_dgsa.R`** / **`qusage_dgsa.R`** — loop `run_dearseq_comparison()` / `run_qusage_comparison()` over every valid comparison at each method's baseline specification (bolded options in Table 2.1), checkpointing results incrementally to `output/results/reanalysis/{method}_dgsa_results_list.rds`.
- **`process_dearseq_dgsa_results.R`** / **`process_qusage_dgsa_results.R`** — `build_tidy_dgsa_results()` (tidy + p-value adjustment for **6 methods** — `holm`, `hochberg`, `hommel`, `bonferroni`, `BH`, `BY` — at **3 scopes** — `global`, `withinTime`, `withinComparison`, this multiplicity itself being part of the robustness assessment) followed by `assign_dgsa_colours()`, saved to `output/results/reanalysis/{method}_dgsa_results_processed.rds`.
- **`dgsa_comparison_example.R`** — Stage 1: combines both methods' processed results and runs `compute_concordance_metrics()`, saving a summary (overall + per-comparison concordance) to `output/results/reanalysis/dgsa_comparison_summary.rds`.
- **`plot_circos.R`** — `circlize`-based circos plots comparing the two methods' significant gene sets (by aggregate/functional category, direction of regulation, and correlation with immune response) across vaccines and timepoints, with configurable p-value correction method/scope, significance threshold, and score-filtering options — used to visualise where the two DGSA methods agree or diverge (the robustness question).

### Combined results dataframe schema

Both `dearseq_dgsa_results_processed.rds` and the QuSAGE equivalent share a common schema so they can be row-bound for comparison: `comparison`, `condition` (vaccine), `time` (day post-vax), `condition.colour`, `gs.name`, `gs.description`, `gs.name.description`, `gs.aggregate`, `gs.colour`, `activation.score`, `fc.score`, `mean.corr`, `corr.mean`, `rawPval`, `method` (`"dearseq"` / `"qusage"`), plus `{scope}.adjPval_{correction}` columns for every scope x correction combination.

### 4. Specification analysis (`analysis/specification_analysis/`)

Chapter 2, Section 2.2.4's robustness/specification analysis: evaluates the 35,640-specification grid from Table 2.1 (see `R/specifications.R` above) and computes the robustness metric pi_{g,v,j} for every gene set x vaccine x timepoint comparison, without ever running DGSA 35,640 times or materialising a table of that size (see `R/robustness_metrics.R` above for how). Three numbered driver scripts, run in order:

- **`01_build_specification_grid.R`** — builds and saves the raw (66) and post-hoc (540) specification grids to `output/results/specification_analysis/`.
- **`02_run_raw_specifications.R`** — **the expensive step.** Runs each of the 66 raw specifications across every valid comparison (dearseq's permutation test in particular is not cheap), checkpointing each specification's results to its own file (`output/results/specification_analysis/raw/{spec_label}.rds`); an interrupted run resumes rather than restarts, both across specifications and within one. Has a `SMOKE_TEST` switch at the top to restrict to a handful of comparisons first, to confirm the whole pipeline runs end-to-end before committing a laptop to the full run.
- **`03_apply_posthoc_and_robustness.R`** — tidies and p-value-adjusts each raw specification's results (`R/postprocessing.R`) and folds them into the robustness-metric accumulator (`R/robustness_metrics.R`), checkpointed by which raw specifications have been accumulated so far (`output/results/specification_analysis/robustness_accumulator_state.rds`). Safe to re-run at any point, including while `02` is still producing more results — saves `output/results/specification_analysis/robustness_metrics.rds`, warning if the result is still partial.

A fourth script producing the **specification heatmap** (the visualisation summary tool introduced alongside the robustness metric) is intentionally not yet built — its design hasn't been finalised.

## Running the pipeline

```r
# from the project root, with data-raw/ populated
source("analysis/analysis_master.R")
```

This currently runs preprocessing end-to-end (`preprocessing_master.R`); the reanalysis (`reanalysis_master.R`, which runs `dearseq_dgsa.R` → `process_dearseq_dgsa_results.R` → `qusage_dgsa.R` → `process_qusage_dgsa_results.R` → `dgsa_comparison_example.R` → `plot_circos.R`) and descriptive (`descriptive_master.R`) scripts can be sourced separately once preprocessing has produced the files in `data/`. The dearseq and QuSAGE driver scripts are long-running and checkpoint their results to `output/results/` so interrupted runs can resume rather than restart.
