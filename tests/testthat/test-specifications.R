test_that("build_raw_specification_grid() produces 64 dearseq + 2 qusage specifications", {
  raw_grid <- build_raw_specification_grid()

  expect_equal(nrow(raw_grid), 66)
  expect_equal(sum(raw_grid$method == "dearseq"), 64)
  expect_equal(sum(raw_grid$method == "qusage"), 2)
  expect_equal(length(unique(raw_grid$raw_spec_id)), 66)
  expect_equal(length(unique(raw_grid$spec_label)), 66)
})

test_that("build_raw_specification_grid() covariate sets are all 16 subsets of the pool", {
  raw_grid <- build_raw_specification_grid()
  covariate_pool <- c("age_imputed", "gender", "study_accession", "race")

  # 16 covariate subsets x 2 weight methods x 2 weight levels = 64 rows
  covariate_lengths <- raw_grid |>
    dplyr::filter(method == "dearseq") |>
    dplyr::pull(covariates) |>
    unique()
  expect_equal(length(covariate_lengths), 16)

  # every subset size is represented, with the right counts (choose(4,k))
  sizes <- vapply(covariate_lengths, length, integer(1))
  expect_equal(sort(as.integer(table(sizes))), sort(choose(length(covariate_pool), 0:length(covariate_pool))))
})

test_that("build_raw_specification_grid() flags exactly one baseline per method", {
  raw_grid <- build_raw_specification_grid()

  baseline <- baseline_raw_specifications(raw_grid)
  expect_equal(nrow(baseline), 2)
  expect_setequal(baseline$method, c("dearseq", "qusage"))

  dearseq_baseline <- dplyr::filter(baseline, method == "dearseq")
  expect_true(setequal(dearseq_baseline$covariates[[1]],
                       c("age_imputed", "gender", "study_accession", "race")))
  expect_equal(dearseq_baseline$which_weights, "loclin")
  expect_false(dearseq_baseline$gene_based_weights)

  qusage_baseline <- dplyr::filter(baseline, method == "qusage")
  expect_false(qusage_baseline$equal_variance)
})

test_that("build_posthoc_specification_grid() produces 540 specifications with one baseline", {
  posthoc_grid <- build_posthoc_specification_grid()

  expect_equal(nrow(posthoc_grid), 540)
  expect_equal(sum(posthoc_grid$is_baseline), 1)

  baseline <- dplyr::filter(posthoc_grid, is_baseline)
  expect_equal(baseline$adjustment_scope, "withinTime")
  expect_equal(baseline$adjustment_method, "BH")
  expect_equal(baseline$alpha, 0.05)
  expect_equal(baseline$fc_threshold, 0)
})

test_that("build_full_specification_grid() produces exactly 35,640 specifications with one baseline", {
  full_grid <- build_full_specification_grid()

  expect_equal(nrow(full_grid), 35640)
  expect_equal(sum(full_grid$is_baseline), 1)
  expect_equal(length(unique(full_grid$spec_id)), 35640)
})
