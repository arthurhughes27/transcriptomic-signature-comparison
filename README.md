# transcriptomic-signature-comparison

A reproducible R pipeline for a multi-study, multi-vaccine analysis of the HIPC ImmuneSignatures2 (IS2) dataset (Chapter 2 of a PhD thesis). It has two goals:

1. **Comparative transcriptomic biomarker discovery** — identify gene-set-level ("modular") transcriptomic signatures of vaccination across 13 vaccines and ~30 studies, at multiple post-vaccination timepoints, and relate them to antibody response.
2. **Robustness analysis** — evaluate how sensitive those signatures are to methodological choices (differential gene-set analysis method, multiple-testing correction, significance thresholds) by re-running the analysis across a large grid of analytical specifications and two independent DGSA methods, `dearseq` and `QuSAGE`.

Analysis works at the level of curated gene sets (Blood Transcriptional Modules / BloodGen3 Modules) rather than individual genes.

**Note on data:** raw/processed data (`data-raw/`, `data/`) and generated results (`output/`) are gitignored, not stored in this repository. All scripts assume these folders exist locally, populated by the user.

## Repository structure

- **`R/`** — modular, reusable analysis functions (sourced via `R/load_all.R`, not a package): DGSA logic for each method, the specification grid, p-value adjustment, the robustness metric and its plots, and shared plotting helpers.
- **`tests/`** — testthat suite for `R/` (`Rscript tests/run_tests.R`), runnable against a synthetic dataset without real data.
- **`analysis/`** — driver scripts, organised by pipeline stage:
  - `preprocessing/` — raw data (`data-raw/`) → analysis-ready tables (`data/`).
  - `descriptive/` — dataset-description figures/tables for the thesis chapter.
  - `reanalysis/` — Stage 1: baseline dearseq vs QuSAGE comparison (circos plots, heatmaps, concordance metrics).
  - `specification_analysis/` — Stage 2: the full specification grid and the robustness metric it produces (heatmaps, distribution plots, summary tables).
- **`manuscript/figures/`** — diagrams for the write-up.

## Running the pipeline

```r
# from the project root, with data-raw/ populated
source("analysis/analysis_master.R")
```

This runs preprocessing end-to-end. The `reanalysis_master.R`, `descriptive_master.R`, and `specification_analysis/` scripts (numbered, run in order) can be sourced separately once `data/` has been produced. Long-running scripts checkpoint their results to `output/results/` so interrupted runs resume rather than restart.
