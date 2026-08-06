test_that("list_valid_comparisons() finds vaccine x timepoint pairs with paired data", {
  hipc <- tibble::tibble(
    vaccine_name       = c("V1", "V1", "V1", "V2", "V2"),
    participant_id     = c("p1", "p1", "p2", "p3", "p3"),
    time_post_last_vax = c(-1,   7,    7,    -1,   3)
  )
  # p1: pre (-1) + day 7  -> V1/day7 valid
  # p2: day 7 only, no pre -> doesn't invalidate V1/day7 (p1 already valid)
  # p3: pre (-1) + day 3  -> V2/day3 valid

  out <- list_valid_comparisons(hipc)

  expect_equal(nrow(out), 2)
  expect_setequal(paste(out$vaccine_name, out$day), c("V1 7", "V2 3"))
})

test_that("filter_paired_samples() keeps only the most recent pre-vaccination sample", {
  hipc <- tibble::tibble(
    vaccine_name          = c("V1",  "V1", "V1", "V1"),
    participant_id        = c("p1",  "p1", "p1", "p2"),
    time_post_last_vax    = c(-10,   -1,   7,    7),
    study_time_collected  = c(-10,   -1,   7,    7)
  )
  # p2 has no pre-vaccination sample at all, so should be excluded entirely

  out <- filter_paired_samples(hipc, "V1", 7)

  expect_equal(nrow(out), 2)
  expect_setequal(out$time_post_last_vax, c(-1, 7))
  expect_true(all(out$participant_id == out$participant_id[1]))  # only p1 remains
})

test_that("build_covariate_matrix() supports an empty covariate set", {
  df <- tibble::tibble(
    vaccine_code = factor(c(1, 1, 2, 2)),
    age_imputed  = c(30, 40, 50, 60)
  )
  x <- build_covariate_matrix(df, covariates = character(0))

  expect_true("(Intercept)" %in% colnames(x))
  expect_false(any(grepl("^vaccine_code", colnames(x))))
  expect_false(any(grepl("age_imputed", colnames(x))))
})

test_that("build_covariate_matrix() drops linearly dependent covariate columns", {
  df <- tibble::tibble(
    vaccine_code = factor(c(1, 1, 2, 2)),
    age_years    = c(30, 40, 50, 60),
    age_months   = c(30, 40, 50, 60) * 12  # perfectly collinear with age_years
  )
  x <- build_covariate_matrix(df, covariates = c("age_years", "age_months"))

  # Both age_years and age_months cannot survive (they're collinear); exactly
  # one adjustment covariate plus the intercept should remain.
  expect_equal(ncol(x), 2)
})

test_that("calculate_gs_correlation() recovers a known correlation structure", {
  set.seed(1)
  response <- rnorm(20)
  y <- rbind(
    gene1 = response * 2 + rnorm(20, sd = 0.01),  # near-perfectly correlated
    gene2 = rnorm(20)                              # unrelated
  )
  colnames(y) <- paste0("s", 1:20)

  GSA <- list(
    genesets      = list(c("gene1"), c("gene2")),
    geneset.names = c("gs1", "gs2")
  )

  out <- calculate_gs_correlation(y, response, GSA)

  expect_gt(out$mean.corr[["gs1"]], 0.9)
  expect_lt(abs(out$mean.corr[["gs2"]]), 0.9)
})

test_that("calculate_gs_correlation() errors on mismatched sample counts", {
  y <- matrix(rnorm(10), nrow = 2, dimnames = list(c("gene1", "gene2"), NULL))
  GSA <- list(genesets = list(c("gene1")), geneset.names = c("gs1"))

  expect_error(calculate_gs_correlation(y, response = rnorm(3), GSA = GSA))
})
