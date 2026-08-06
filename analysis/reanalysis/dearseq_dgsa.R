# =============================================================================
# Differential Gene-Set Analysis with dearseq — baseline specification
# =============================================================================
# Thin driver: runs the baseline dearseq specification (covariates = age,
# sex, study, race; mean-variance weights = local linear regression at the
# observation level - the bolded options in Table 2.1) across every vaccine
# x timepoint comparison, comparing post-vaccination samples against
# self-controls (pre-vaccination baseline). All the actual logic lives in
# R/dgsa_common.R and R/dgsa_dearseq.R.
# Results are saved incrementally to guard against mid-run crashes.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)
library(dearseq)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_expr  <- fs::path("data", "hipc_merged_young_noNorm.rds")
p_data_btm   <- fs::path("data", "BTM_processed.rds")
p_results    <- fs::path("output", "results", "reanalysis", "dearseq_dgsa_results_list.rds")

# ── Load data ─────────────────────────────────────────────────────────────────

hipc <- readRDS(p_data_expr)
BTM  <- readRDS(p_data_btm)

gene_names <- hipc %>% select(a1cf:zzz3) %>% colnames()

# Label pre- (1) and post-vaccination (2) samples; factorise study
hipc <- hipc %>%
  mutate(
    vaccine_code    = as.factor(if_else(time_post_last_vax > 0, 2, 1)),
    study_accession = as.factor(study_accession)
  )

# =============================================================================
# MAIN LOOP
# =============================================================================

# Load any previously saved results so a restart continues rather than reruns
results_list <- if (file_exists(p_results)) readRDS(p_results) else list()

comparisons <- list_valid_comparisons(hipc)

for (i in seq_len(nrow(comparisons))) {

  vax <- comparisons$vaccine_name[i]
  day <- comparisons$day[i]

  comparison_name <- sprintf("%s vs Control - Day %s", vax, day)

  # Skip comparisons already completed in a previous run
  if (comparison_name %in% names(results_list)) {
    message("Skipping (already complete): ", comparison_name)
    next
  }

  message("Running: ", comparison_name)

  result <- tryCatch(
    run_dearseq_comparison(vax, day, hipc, BTM, gene_names),
    error = function(e) {
      message("  ERROR: ", conditionMessage(e))
      NULL
    }
  )

  results_list[[comparison_name]] <- result

  # Save after every comparison so a crash loses at most one result
  saveRDS(object = results_list,
          file = p_results)

  # Free memory before the next iteration
  gc()
}

message("Analysis complete. Results saved to: ", p_results)

rm(list = ls())
