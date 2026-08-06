# =============================================================================
# Specification analysis — 01: build specification grid
# =============================================================================
# Builds and saves the raw (66) and post-hoc (540) specification grids
# (R/specifications.R), so downstream scripts (02, 03) and manual
# inspection share a single, versioned copy rather than rebuilding it ad
# hoc. See Chapter 2, Table 2.1.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

out_dir <- fs::path("output", "results", "specification_analysis")
fs::dir_create(out_dir)

# =============================================================================
# BUILD AND SAVE
# =============================================================================

raw_grid     <- build_raw_specification_grid()
posthoc_grid <- build_posthoc_specification_grid()

message(sprintf(
  "Raw specifications: %d (%d dearseq, %d qusage)",
  nrow(raw_grid), sum(raw_grid$method == "dearseq"), sum(raw_grid$method == "qusage")
))
message(sprintf("Post-hoc specifications: %d", nrow(posthoc_grid)))
message(sprintf("Full specification grid size: %d", nrow(raw_grid) * nrow(posthoc_grid)))

saveRDS(raw_grid, file = fs::path(out_dir, "raw_specification_grid.rds"))
saveRDS(posthoc_grid, file = fs::path(out_dir, "posthoc_specification_grid.rds"))

message("Specification grids saved to: ", out_dir)

rm(list = ls())
