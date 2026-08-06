# =============================================================================
# Synthetic IS2-like dataset for testing
# =============================================================================
# Generates a small dataset mimicking the schema of data/hipc_merged_*.rds
# (clinical + immune response + expression, one row per sample) and
# data/BTM_processed.rds (gene set object), so that the modular DGSA
# pipeline in R/ can be tested without access to the real (private) IS2 data.
# =============================================================================

#' Generate a small synthetic IS2-like dataset for testing
#'
#' Calling this function sets the global RNG seed (via `seed`) for
#' reproducibility.
#'
#' Simulates two vaccines, each administered in two studies, with
#' participants having a pre-vaccination baseline (a few with two baseline
#' measurements, to exercise "most recent baseline" logic) and
#' post-vaccination samples at days 1 and 7. Genes are grouped into four
#' (partially overlapping) synthetic gene sets. If `effect_vaccine` is set, a
#' true post-vaccination effect is injected into one gene set for that
#' vaccine at `effect_day`, so at least one comparison is non-null; all
#' other comparisons are null (noise only).
#'
#' @param seed Integer RNG seed.
#' @param n_participants_per_study Number of participants simulated per
#'   study.
#' @param n_genes Number of synthetic genes to simulate (minimum 40, so all
#'   four gene sets are populated).
#' @param effect_vaccine Name of the vaccine given a synthetic true effect,
#'   or NULL for an entirely null dataset.
#' @param effect_day Post-vaccination day at which the effect above is
#'   injected.
#'
#' @return A list with elements `hipc` (tibble of samples, matching the
#'   hipc_merged_*.rds schema) and `BTM` (gene set list, matching the
#'   BTM_processed.rds schema).
make_synthetic_is2_data <- function(seed = 42,
                                    n_participants_per_study = 12,
                                    n_genes = 40,
                                    effect_vaccine = "Vaccine A",
                                    effect_day = 7) {
  stopifnot(n_genes >= 40)
  set.seed(seed)

  vaccines <- c("Vaccine A", "Vaccine B")
  studies_per_vaccine <- list(
    "Vaccine A" = c("SIM001", "SIM002"),
    "Vaccine B" = c("SIM003", "SIM004")
  )
  post_days <- c(1, 7)

  gene_names <- sprintf("gene%02d", seq_len(n_genes))

  # ── Gene sets (deliberately overlapping, as real gene sets are) ─────────
  geneset_defs <- list(
    geneset1 = gene_names[1:10],
    geneset2 = gene_names[8:18],
    geneset3 = gene_names[19:28],
    geneset4 = gene_names[29:40]
  )

  BTM <- list(
    genesets                    = unname(geneset_defs),
    geneset.names                = sprintf("Geneset %d (S%d)", seq_along(geneset_defs), seq_along(geneset_defs)),
    geneset.descriptions         = sprintf("Synthetic geneset %d", seq_along(geneset_defs)),
    geneset.names.descriptions   = sprintf("Synthetic geneset %d (S%d)", seq_along(geneset_defs), seq_along(geneset_defs)),
    geneset.aggregates           = factor(rep("Synthetic", length(geneset_defs)))
  )

  # ── Clinical + expression rows ──────────────────────────────────────────
  rows <- list()

  for (vax in vaccines) {
    for (study in studies_per_vaccine[[vax]]) {
      for (p in seq_len(n_participants_per_study)) {

        participant_id <- paste0(study, "_SUB", p)
        age             <- round(stats::runif(1, 18, 65), 1)
        gender          <- sample(c("Male", "Female", "Unknown"), 1)
        race            <- sample(c("White", "Asian", "Black or African American", "Unknown"), 1)
        response        <- stats::rnorm(1)

        # Most participants have one pre-vaccination sample; every 4th has
        # two (at different study_time_collected), to exercise the
        # "most recent baseline" logic in filter_paired_samples().
        pre_times <- if (p %% 4 == 0) c(-10, -1) else -1

        for (t in pre_times) {
          rows[[length(rows) + 1]] <- make_synthetic_sample_row(
            participant_id, vax, study, t, t, age, gender, race, NA_real_,
            gene_names
          )
        }

        for (day in post_days) {
          is_effect <- !is.null(effect_vaccine) && vax == effect_vaccine && day == effect_day
          rows[[length(rows) + 1]] <- make_synthetic_sample_row(
            participant_id, vax, study, day, day, age, gender, race, response,
            gene_names,
            effect_genes = if (is_effect) geneset_defs$geneset1 else character(0)
          )
        }
      }
    }
  }

  hipc <- dplyr::bind_rows(rows)

  list(hipc = hipc, BTM = BTM)
}

# Builds one sample row (clinical + immune response + gene expression
# columns). `effect_genes` (if any) are given a fixed upward expression
# shift, simulating a true post-vaccination effect.
make_synthetic_sample_row <- function(participant_id, vaccine_name, study_accession,
                                      time_post_last_vax, study_time_collected,
                                      age, gender, race, response_log2_mfc,
                                      gene_names, effect_genes = character(0)) {

  expr <- stats::rnorm(length(gene_names), mean = 6, sd = 1)
  names(expr) <- gene_names

  if (length(effect_genes) > 0) {
    expr[effect_genes] <- expr[effect_genes] + 1.5
  }

  row <- tibble::tibble(
    participant_id                = participant_id,
    vaccine_name                   = vaccine_name,
    study_accession                = study_accession,
    time_post_last_vax             = time_post_last_vax,
    study_time_collected           = study_time_collected,
    age_imputed                    = age,
    gender                         = gender,
    race                           = race,
    immResp_MFC_anyAssay_log2_MFC  = response_log2_mfc
  )

  dplyr::bind_cols(row, tibble::as_tibble(as.list(expr)))
}
