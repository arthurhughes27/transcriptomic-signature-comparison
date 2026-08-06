test_that("count_significant_specifications() matches a brute-force computation", {
  set.seed(1)
  df <- tibble::tibble(
    gs.name             = paste0("gs", 1:6),
    condition           = "V1",
    time                = 7,
    fc.score            = c(0.1, 0.3, 0.5, 0.7, 0.9, NA),
    global.adjPval_BH   = c(0.001, 0.02, 0.04, 0.2, 0.5, 0.9),
    global.adjPval_holm = c(0.2, 0.03, NA, 0.01, 0.4, 0.05)
  )

  alphas        <- c(0.01, 0.05, 0.1)
  fc_thresholds <- c(0.2, 0.6)

  fast <- count_significant_specifications(df, alphas = alphas, fc_thresholds = fc_thresholds)

  adjpval_cols <- c("global.adjPval_BH", "global.adjPval_holm")
  brute <- rep(0L, nrow(df))
  for (col in adjpval_cols) {
    for (a in alphas) {
      for (f in fc_thresholds) {
        sig <- (df[[col]] <= a) & (abs(df$fc.score) >= f)
        sig[is.na(sig)] <- FALSE
        brute <- brute + as.integer(sig)
      }
    }
  }

  expect_equal(fast, brute)
  expect_true(all(fast >= 0 & fast <= length(adjpval_cols) * length(alphas) * length(fc_thresholds)))
})

test_that("count_significant_specifications() errors without adjusted p-value columns", {
  df <- tibble::tibble(gs.name = "gs1", fc.score = 1)
  expect_error(count_significant_specifications(df))
})

test_that("accumulate_robustness_counts() and compute_robustness_metric() combine multiple raw runs correctly", {
  run1 <- tibble::tibble(
    gs.name           = c("gs1", "gs2"),
    condition         = "V1",
    time              = 7,
    fc.score          = c(1, 1),
    global.adjPval_BH = c(0.01, 0.5)
  )
  run2 <- tibble::tibble(
    gs.name           = c("gs1"),  # gs2 not evaluated in this run
    condition         = "V1",
    time              = 7,
    fc.score          = c(1),
    global.adjPval_BH = c(0.01)
  )

  alphas <- c(0.05)
  fc_thresholds <- c(0)

  acc <- accumulate_robustness_counts(NULL, run1, alphas = alphas, fc_thresholds = fc_thresholds)
  acc <- accumulate_robustness_counts(acc, run2, alphas = alphas, fc_thresholds = fc_thresholds)

  out <- compute_robustness_metric(acc)

  gs1 <- dplyr::filter(out, gs.name == "gs1")
  gs2 <- dplyr::filter(out, gs.name == "gs2")

  # gs1: significant in both runs (1 posthoc spec each) -> 2/2 = 1
  expect_equal(gs1$n_evaluated, 2)
  expect_equal(gs1$robustness, 1)

  # gs2: only evaluated in run1, not significant (p=0.5 > 0.05) -> 0/1 = 0
  expect_equal(gs2$n_evaluated, 1)
  expect_equal(gs2$robustness, 0)
})

test_that("accumulate_robustness_counts() is order-independent", {
  run1 <- tibble::tibble(
    gs.name = "gs1", condition = "V1", time = 7,
    fc.score = 1, global.adjPval_BH = 0.01
  )
  run2 <- tibble::tibble(
    gs.name = "gs1", condition = "V1", time = 7,
    fc.score = 1, global.adjPval_BH = 0.9
  )

  acc_ab <- accumulate_robustness_counts(NULL, run1)
  acc_ab <- accumulate_robustness_counts(acc_ab, run2)

  acc_ba <- accumulate_robustness_counts(NULL, run2)
  acc_ba <- accumulate_robustness_counts(acc_ba, run1)

  expect_equal(
    dplyr::arrange(acc_ab, gs.name),
    dplyr::arrange(acc_ba, gs.name)
  )
})
