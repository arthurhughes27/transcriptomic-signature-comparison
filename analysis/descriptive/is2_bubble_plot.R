# =============================================================================
# IS2 Dataset — Vaccine x Timepoint Sample Bubble Plot
# =============================================================================
# Produces the main-text figure describing the IS2 dataset (Chapter 2,
# Section 2.3.1): vaccines on the y-axis, days post-vaccination on the
# x-axis, bubble size proportional to the number of transcriptomic samples
# available for that vaccine x timepoint combination (summed across every
# contributing study). The per-study equivalent of this figure has been
# moved to Appendix A (analysis/descriptive/is2_appendix_descriptives.R),
# since a reader wants the vaccine-level picture first.
#
# The days actually used in the DGSA analysis grid (DAYS_TO_HIGHLIGHT,
# below) are highlighted with a light background band, in the same
# per-day colours used throughout the specification-analysis figures
# (R/plot_helpers.R's assign_day_colours()/day_highlight_bands()).
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(dplyr)
library(ggplot2)
library(fs)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

processed_data_folder      <- "data"
descriptive_figures_folder <- fs::path("output", "figures", "descriptive")

fs::dir_create(descriptive_figures_folder)

# Timepoints highlighted on the bubble plots below - matches DAYS_TO_ANALYSE
# in the reanalysis/specification-analysis driver scripts.
DAYS_TO_HIGHLIGHT <- c(1, 3, 7)

# ── Load data ─────────────────────────────────────────────────────────────────

hipc_merged_all_norm <- readRDS(
  fs::path(processed_data_folder, "hipc_merged_all_norm.rds")
)

# Vaccine order fixed at preprocessing time (analysis/preprocessing/preprocessing_clinical.R)
vaccine_order <- levels(hipc_merged_all_norm$vaccine_name)

# =============================================================================
# BUBBLE PLOT: SAMPLE COUNTS PER VACCINE x TIMEPOINT
# =============================================================================

counts <- hipc_merged_all_norm %>%
  group_by(vaccine_name, vaccine_colour, time_post_last_vax) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    vaccine_name = factor(vaccine_name, levels = vaccine_order),
    # Order time points numerically
    time_post_last_vax = factor(
      as.character(time_post_last_vax),
      levels = sort(unique(as.numeric(as.character(time_post_last_vax)))) %>%
        as.character(),
      ordered = TRUE
    ),
    # Compress bubble sizes with a sublinear transform so that a small number
    # of large vaccines do not swamp the size scale; true counts are recovered
    # via the legend breaks/labels below.
    size_var = n ^ (2 / 3)
  )

size_breaks_counts <- c(10, 50, 100, 200, 400)
size_breaks         <- size_breaks_counts ^ (2 / 3)

bubble_plot <- ggplot(counts, aes(x = time_post_last_vax, y = vaccine_name)) +
  day_highlight_bands(levels(counts$time_post_last_vax), DAYS_TO_HIGHLIGHT) +
  geom_point(
    aes(size = size_var, fill = vaccine_name),
    shape = 21, colour = "black", alpha = 0.75, show.legend = c(size = TRUE, fill = FALSE)
  ) +
  geom_text(
    aes(label = n),
    colour = "white", size = 3.5, vjust = 0.5, show.legend = FALSE
  ) +
  scale_fill_manual(values = setNames(counts$vaccine_colour, counts$vaccine_name)) +
  scale_size_area(
    name     = "Count",
    max_size = 28,
    breaks   = size_breaks,
    labels   = size_breaks_counts,
    guide    = guide_legend(override.aes = list(fill = "grey80", colour = "black"))
  ) +
  scale_y_discrete(limits = rev(vaccine_order)) +
  labs(
    x     = "Days post-vaccination",
    y     = "Vaccine",
    title = "Participants with transcriptomic samples per vaccine across time"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major.y   = element_line(color = "grey90"),
    panel.grid.minor     = element_blank(),
    axis.title           = element_text(size = 30),
    axis.text            = element_text(size = 14),
    axis.text.x          = element_text(angle = 45, hjust = 1),
    plot.title           = element_text(size = 35, hjust = 0.5, face = "bold"),
    plot.subtitle        = element_text(size = 15, hjust = 0.5),
    legend.title         = element_text(size = 20, hjust = 0.5),
    legend.key.spacing.y = unit(0.3, "cm"),
    legend.spacing.y     = unit(1.0, "cm")
  )

print(bubble_plot)

ggsave(
  filename = "vaccine_bubble_plot.pdf",
  path     = descriptive_figures_folder,
  plot     = bubble_plot,
  width    = 40, height = 25, units = "cm"
)

rm(list = ls())
