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
})
