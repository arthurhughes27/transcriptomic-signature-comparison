test_that("make_synthetic_is2_data() returns a well-formed hipc/BTM pair", {
  d <- make_synthetic_is2_data(seed = 1, n_participants_per_study = 4, n_genes = 40)

  expect_type(d, "list")
  expect_named(d, c("hipc", "BTM"))

  required_cols <- c(
    "participant_id", "vaccine_name", "study_accession",
    "time_post_last_vax", "study_time_collected",
    "age_imputed", "gender", "race", "immResp_MFC_anyAssay_log2_MFC"
  )
  expect_true(all(required_cols %in% colnames(d$hipc)))

  # every participant has at least one pre- and one post-vaccination sample
  coverage <- d$hipc |>
    dplyr::group_by(participant_id) |>
    dplyr::summarise(
      has_pre  = any(time_post_last_vax <= 0),
      has_post = any(time_post_last_vax > 0),
      .groups  = "drop"
    )
  expect_true(all(coverage$has_pre))
  expect_true(all(coverage$has_post))

  # BTM: genesets, names, descriptions, aggregates all align in length
  n_genesets <- length(d$BTM$genesets)
  expect_equal(length(d$BTM$geneset.names), n_genesets)
  expect_equal(length(d$BTM$geneset.descriptions), n_genesets)
  expect_equal(length(d$BTM$geneset.aggregates), n_genesets)
  expect_true(all(unlist(d$BTM$genesets) %in% colnames(d$hipc)))
})

test_that("make_synthetic_is2_data() is reproducible for a fixed seed", {
  d1 <- make_synthetic_is2_data(seed = 7, n_participants_per_study = 3, n_genes = 40)
  d2 <- make_synthetic_is2_data(seed = 7, n_participants_per_study = 3, n_genes = 40)
  expect_equal(d1$hipc, d2$hipc)
})

test_that("make_synthetic_is2_data() injects the requested effect", {
  d <- make_synthetic_is2_data(
    seed = 5, n_participants_per_study = 20, n_genes = 40,
    effect_vaccine = "Vaccine A", effect_day = 7
  )

  effect_genes <- d$BTM$genesets[[1]]  # geneset1

  boosted <- d$hipc |>
    dplyr::filter(vaccine_name == "Vaccine A", time_post_last_vax == 7) |>
    dplyr::pull(effect_genes[1])
  unboosted <- d$hipc |>
    dplyr::filter(vaccine_name == "Vaccine B", time_post_last_vax == 7) |>
    dplyr::pull(effect_genes[1])

  expect_gt(mean(boosted), mean(unboosted) + 1)
})
