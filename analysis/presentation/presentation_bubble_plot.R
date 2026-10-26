# =============================================================================
# Presentation figure — simplified vaccine-level sample bubble plot
# =============================================================================
# A simplified version of the main-text vaccine-level bubble plot
# (analysis/descriptive/is2_descriptive.R) for slides: restricted to days
# 1, 3 and 7 (DAYS_TO_ANALYSE in the reanalysis/specification-analysis
# driver scripts) and, for each of those days, to participants who ALSO
# have a baseline (day <= 0) sample - the exact pairing requirement
# filter_paired_samples() (R/dgsa_common.R) applies before a comparison is
# actually run, reused directly here so the counts shown match what the
# DGSA analysis itself uses. Day 0 itself is not shown (it isn't a
# comparison timepoint; its only role here is as the pairing requirement).
#
# All text sizes are bumped up from the main-text figure for slide
# legibility, and the day-highlight bands (day_highlight_bands(),
# R/plot_helpers.R) are dropped - with only days 1/3/7 shown, every column
# would be highlighted, so the bands no longer distinguish anything.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(dplyr)
library(ggplot2)
library(fs)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

processed_data_folder <- "data"
presentation_figures_folder <- fs::path("output", "figures", "presentation")

fs::dir_create(presentation_figures_folder)

# Must match the DAYS_TO_ANALYSE used throughout the reanalysis/
# specification-analysis driver scripts.
DAYS_TO_SHOW <- c(1, 3, 7)

# ── Load data ─────────────────────────────────────────────────────────────────

hipc_merged_all_norm <- readRDS(
  fs::path(processed_data_folder, "hipc_merged_all_norm.rds")
)

# Vaccine order fixed at preprocessing time (analysis/preprocessing/preprocessing_clinical.R)
vaccine_order <- levels(hipc_merged_all_norm$vaccine_name)

vaccine_colours <- hipc_merged_all_norm |>
  distinct(vaccine_name, vaccine_colour)

# =============================================================================
# PAIRED SAMPLE COUNTS PER VACCINE x DAY (1, 3, 7 ONLY)
# =============================================================================

# Every (vaccine, day) combination with any row at all at that day,
# restricted to DAYS_TO_SHOW - not yet restricted to paired participants;
# that's applied per-combination below via filter_paired_samples(), with
# combinations left with zero paired samples dropped afterwards.
vaccine_day_combos <- hipc_merged_all_norm |>
  filter(time_post_last_vax %in% DAYS_TO_SHOW) |>
  distinct(vaccine_name, time_post_last_vax)

counts <- purrr::pmap_dfr(vaccine_day_combos, function(vaccine_name, time_post_last_vax) {
  paired <- filter_paired_samples(hipc_merged_all_norm, vaccine_name, time_post_last_vax)
  tibble::tibble(
    vaccine_name        = vaccine_name,
    time_post_last_vax  = time_post_last_vax,
    n                    = sum(paired$time_post_last_vax == time_post_last_vax)
  )
}) |>
  filter(n > 0) |>
  left_join(vaccine_colours, by = "vaccine_name") |>
  mutate(
    vaccine_name = factor(vaccine_name, levels = vaccine_order),
    time_post_last_vax = factor(
      as.character(time_post_last_vax),
      levels = as.character(sort(DAYS_TO_SHOW)),
      ordered = TRUE
    ),
    # Same sublinear bubble-size compression as the main-text figure.
    size_var = n ^ (2 / 3)
  )

size_breaks_counts <- c(10, 50, 100, 200, 400)
size_breaks         <- size_breaks_counts ^ (2 / 3)

bubble_plot <- ggplot(counts, aes(x = time_post_last_vax, y = vaccine_name)) +
  geom_point(
    aes(size = size_var, fill = vaccine_name),
    shape = 21, colour = "black", alpha = 0.75, show.legend = c(size = TRUE, fill = FALSE)
  ) +
  geom_text(
    aes(label = n),
    colour = "white", size = 6, vjust = 0.5, show.legend = FALSE
  ) +
  scale_fill_manual(values = setNames(counts$vaccine_colour, counts$vaccine_name)) +
  scale_size_area(
    name     = "Samples",
    max_size = 32,
    breaks   = size_breaks,
    labels   = size_breaks_counts,
    guide    = guide_legend(override.aes = list(fill = "grey80", colour = "black"))
  ) +
  scale_y_discrete(limits = rev(vaccine_order)) +
  labs(
    x     = "Days post vaccination",
    y     = "Vaccine",
    title = "IS2 dataset: transcriptomic samples"
  ) +
  theme_minimal(base_size = 26) +
  theme(
    panel.grid.major.y   = element_line(color = "grey90"),
    panel.grid.minor     = element_blank(),
    axis.title           = element_text(size = 38),
    axis.text            = element_text(size = 24),
    axis.text.x          = element_text(hjust = 1),
    plot.title           = element_text(size = 44, hjust = 0.5, face = "bold"),
    legend.title         = element_text(size = 32, hjust = 0.5, face = "bold"),
    legend.text          = element_text(size = 24),
    legend.key.spacing.y = unit(0.3, "cm"),
    legend.spacing.y     = unit(1.0, "cm")
  )

print(bubble_plot)

ggsave(
  filename = "vaccine_bubble_plot_presentation.pdf",
  path     = presentation_figures_folder,
  plot     = bubble_plot,
  width    = 40, height = 28, units = "cm"
)

rm(list = ls())
