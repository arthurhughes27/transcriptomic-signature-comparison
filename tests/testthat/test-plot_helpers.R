test_that("assign_colours() uses provided colours when given", {
  out <- assign_colours(c("a", "b"), provided_colours = c("#111111", "#222222"))
  expect_equal(out, c(a = "#111111", b = "#222222"))
})

test_that("assign_colours() generates a colour per label when none provided", {
  out <- assign_colours(c("a", "b", "c"))
  expect_equal(names(out), c("a", "b", "c"))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", out)))
})

test_that("assign_colours() errors if provided_colours length mismatches labels", {
  expect_error(assign_colours(c("a", "b"), provided_colours = "#111111"))
})

test_that("assign_dgsa_colours() adds condition.colour and gs.colour columns", {
  results_df <- tibble::tibble(
    condition    = factor(c("V1", "V1", "V2"), levels = c("V1", "V2")),
    gs.aggregate = factor(c("A", "B", "A"), levels = c("A", "B"))
  )

  out <- assign_dgsa_colours(
    results_df,
    condition_colors = c("#000001", "#000002"),
    aggregate_colors = c("#100001", "#100002")
  )

  expect_true(all(c("condition.colour", "gs.colour") %in% colnames(out)))
  expect_equal(out$condition.colour, c("#000001", "#000001", "#000002"))
  expect_equal(out$gs.colour, c("#100001", "#100002", "#100001"))
})

test_that("default configuration vectors are aligned", {
  expect_equal(length(default_conditions_order()), length(default_condition_colors()))
  expect_equal(length(default_day_order()), length(default_day_colors()))
})

test_that("default_method_colors() gives one colour each for dearseq and qusage", {
  out <- default_method_colors()
  expect_equal(names(out), c("dearseq", "qusage"))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", out)))
})

test_that("assign_day_colours() gives the same colour to a day regardless of which other days are present", {
  colour_full   <- assign_day_colours(c(1, 3, 7, 10, 14, 21))
  colour_subset <- assign_day_colours(c(3, 7))

  expect_equal(colour_full[["Day 3"]], colour_subset[["Day 3"]])
  expect_equal(colour_full[["Day 7"]], colour_subset[["Day 7"]])
})

test_that("assign_day_colours() restricts output to the days present in times", {
  out <- assign_day_colours(c(1, 7))
  expect_equal(names(out), c("Day 1", "Day 7"))
})

test_that("assign_day_colours() handles days outside default_day_order() by ramping extra colours", {
  out <- assign_day_colours(c(1, 3, 100))
  expect_equal(names(out), c("Day 1", "Day 3", "Day 100"))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", out)))
})

test_that("day_facet_strip() builds a strip object without erroring", {
  testthat::skip_if_not_installed("ggh4x")
  strip <- day_facet_strip(c(1, 3, 7))
  expect_false(is.null(strip))
})

test_that("day_highlight_bands() returns one band per highlighted day present, positioned by axis order", {
  testthat::skip_if_not_installed("ggplot2")

  x_levels <- c("1", "3", "7", "10")
  bands <- day_highlight_bands(x_levels, highlight_days = c(1, 7))

  expect_equal(length(bands), 2)
  # Day 1 is at x-position 1, Day 7 at x-position 3.
  expect_equal(bands[[1]]$data$xmin, 1 - 0.5)
  expect_equal(bands[[2]]$data$xmin, 3 - 0.5)
})

test_that("day_highlight_bands() skips highlighted days absent from x_levels", {
  testthat::skip_if_not_installed("ggplot2")

  bands <- day_highlight_bands(c("1", "3"), highlight_days = c(1, 3, 7))

  expect_equal(length(bands), 2)
})

test_that("day_highlight_bands() uses the same colour as assign_day_colours()", {
  testthat::skip_if_not_installed("ggplot2")

  bands <- day_highlight_bands(c("1"), highlight_days = 1)

  expect_equal(bands[[1]]$data$fill, unname(assign_day_colours(1)[["Day 1"]]))
})
