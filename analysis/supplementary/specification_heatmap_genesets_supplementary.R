# =============================================================================
# Specification analysis — gene-set-level robustness heatmap (supplementary)
# =============================================================================
# Gene-set-level companion to the main-text aggregate heatmap
# (analysis/specification_analysis/04_specification_heatmap_aggregate.R):
# one row per gene set instead of one per aggregate, using the robustness
# metric pi_{g,v,j} produced by 03_apply_posthoc_and_robustness.R
# (output/results/specification_analysis/robustness_metrics.rds) - see
# R/robustness_heatmaps.R for the plotting function itself.
#
# With 258 gene sets, this is paginated (ROWS_PER_PAGE below) into a single
# multi-page PDF (save_multi_page_pdf(), R/plot_helpers.R) rather than one
# very tall page, so every page is a normal, consistently-sized page
# instead of requiring heavy zooming. Pages are sized to
# GENESET_PAGE_WIDTH_CM/HEIGHT_CM below, calibrated to the true A4 aspect
# ratio (1:sqrt(2)) at GENESET_PAGE_WIDTH_CM's scale rather than a literal
# 21x29.7cm sheet, so the wide comparison-column layout still has room.
# ROWS_PER_PAGE is scaled to match: at the previous 45cm x 40cm page,
# 45 rows only filled ~60% of the page height (row density ~0.53cm/row);
# ROWS_PER_PAGE below targets ~92% fill of the new, taller page at that
# same row density/text size - full names/descriptions, at the same
# readable row height as before, just more of them per page. Re-tune if a
# generated page over/underfills once you can see the actual PDF.
#
# Safe to run against partial results (e.g. while
# 02_run_raw_specifications.R / 03_apply_posthoc_and_robustness.R are still
# producing more of the specification grid) - it just plots whatever is in
# robustness_metrics.rds at the time.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_btm     <- fs::path("data", "BTM_processed.rds")
p_metrics      <- fs::path("output", "results", "specification_analysis", "robustness_metrics.rds")
out_dir        <- fs::path("output", "figures", "supplementary")
p_fig_genesets <- fs::path(out_dir, "robustness_heatmap_genesets.pdf")

if (!fs::file_exists(p_metrics)) {
  stop(
    "Robustness metrics not found at ", p_metrics, ".\n",
    "Run 01_build_specification_grid.R, 02_run_raw_specifications.R and ",
    "03_apply_posthoc_and_robustness.R first."
  )
}

fs::dir_create(out_dir)

# ── Config ────────────────────────────────────────────────────────────────────

ROWS_PER_PAGE <- 110

GENESET_PAGE_WIDTH_CM  <- 45
GENESET_PAGE_HEIGHT_CM <- 45 * sqrt(2)  # true A4 ratio (1:1.414), ~63.6cm

# ── Load data ─────────────────────────────────────────────────────────────────

BTM                  <- readRDS(p_data_btm)
robustness_metrics   <- readRDS(p_metrics)
robustness_annotated <- join_geneset_aggregates(robustness_metrics, BTM)

# =============================================================================
# SUPPLEMENTARY FIGURE — gene-set-level heatmap
# =============================================================================

p_genesets_pages <- plot_robustness_heatmap_genesets(robustness_annotated, rows_per_page = ROWS_PER_PAGE)

save_multi_page_pdf(
  p_genesets_pages, path = p_fig_genesets,
  width = GENESET_PAGE_WIDTH_CM / 2.54, height = GENESET_PAGE_HEIGHT_CM / 2.54
)
message(sprintf(
  "Saved gene-set-level robustness heatmap (%d pages) to: %s",
  length(p_genesets_pages), p_fig_genesets
))

# rm(list = ls())
