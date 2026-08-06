# =============================================================================
# TEMPORARY: validate the new modular QuSAGE implementation against the old
# =============================================================================
# Compares run_qusage_comparison() (R/dgsa_qusage.R) against the existing,
# monolithic analysis/reanalysis/qusage_dgsa.R for one vaccine x timepoint
# comparison, for the p-values and activation scores that come directly out
# of the QuSAGE algorithm (rawPval, activation.score). Correlation with
# immune response is NOT compared: the new implementation deliberately
# reuses calculate_gs_correlation() (the same helper dearseq uses) instead
# of the old script's bespoke weighted multi-study correlation pipeline -
# see the design notes at the top of R/dgsa_qusage.R.
#
# This file, and analysis/reanalysis/qusage_dgsa.R itself, are meant to be
# deleted once you've confirmed the two implementations agree closely
# enough.
#
# Usage:
#   1. If you have not already, run analysis/reanalysis/qusage_dgsa.R once
#      in full (it caches intermediate results under output/results/qusage/,
#      so this is slow only the first time).
#   2. Edit VAX / DAY below to a comparison you want to check (pick
#      something with a small number of studies/participants for speed).
#   3. source(fs::path("analysis", "reanalysis", "validate_qusage_refactor.R"))
# =============================================================================

library(fs)
library(tidyverse)
library(qusage)

source(fs::path("R", "load_all.R"))

# ── Choose a comparison to validate ─────────────────────────────────────────

VAX <- "Yellow Fever (LV)"
DAY <- 7

# ── Old implementation: read from its cached, processed output ─────────────

p_old_results <- fs::path("output", "results", "qusage", "qusage_dgsa_results_processed.rds")

if (!fs::file_exists(p_old_results)) {
  stop(
    "Old QuSAGE results not found at ", p_old_results, ".\n",
    "Run analysis/reanalysis/qusage_dgsa.R once first (it will cache its ",
    "intermediate and final results under output/results/qusage/)."
  )
}

old_results <- readRDS(p_old_results) %>%
  filter(condition == VAX, as.numeric(as.character(time)) == DAY) %>%
  select(gs.name, rawPval_old = rawPval, activation.score_old = activation.score)

if (nrow(old_results) == 0) {
  stop("No old results found for '", VAX, "' - Day ", DAY, ". Check VAX/DAY match a comparison present in qusage_dgsa_results_processed.rds (see unique(readRDS(p_old_results)$condition)).")
}

# ── New implementation: run directly ────────────────────────────────────────

hipc <- readRDS(fs::path("data", "hipc_merged_young_noNorm.rds"))
BTM  <- readRDS(fs::path("data", "BTM_processed.rds"))
gene_names <- hipc %>% select(a1cf:zzz3) %>% colnames()

new_result <- run_qusage_comparison(VAX, DAY, hipc, BTM, gene_names)

if (is.null(new_result)) {
  stop("New implementation returned NULL for '", VAX, "' - Day ", DAY, " (insufficient data after filtering).")
}

new_results <- tibble(
  gs.name                = BTM$geneset.names,
  rawPval_new             = new_result$pvals$rawPval,
  activation.score_new     = unlist(new_result$score$activation.scores)
)

# ── Compare ──────────────────────────────────────────────────────────────────

comparison <- inner_join(old_results, new_results, by = "gs.name") %>%
  mutate(
    pval_diff = abs(rawPval_old - rawPval_new),
    score_diff = abs(activation.score_old - activation.score_new)
  ) %>%
  arrange(desc(pval_diff))

cat(sprintf("\nComparing '%s' - Day %s (%d gene sets matched)\n\n", VAX, DAY, nrow(comparison)))

print(comparison, n = 20)

cat("\n--- Summary ---\n")
cat(sprintf("rawPval:          max abs diff = %.4g, correlation = %.4f\n",
            max(comparison$pval_diff, na.rm = TRUE),
            cor(comparison$rawPval_old, comparison$rawPval_new, use = "complete.obs")))
cat(sprintf("activation.score: max abs diff = %.4g, correlation = %.4f\n",
            max(comparison$score_diff, na.rm = TRUE),
            cor(comparison$activation.score_old, comparison$activation.score_new, use = "complete.obs")))

rm(list = ls())
