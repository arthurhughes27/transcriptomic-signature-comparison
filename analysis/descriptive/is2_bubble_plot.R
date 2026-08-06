# =============================================================================
# IS2 Dataset — Study x Timepoint Sample Bubble Plot
# =============================================================================
# Produces the main-text figure describing the IS2 dataset (Chapter 2,
# Section 2.3.1): studies on the y-axis (coloured by vaccine), days
# post-vaccination on the x-axis, bubble size proportional to the number of
# transcriptomic samples available for that study x timepoint combination.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(dplyr)
library(ggplot2)
library(fs)

# ── Paths ─────────────────────────────────────────────────────────────────────

processed_data_folder      <- "data"
descriptive_figures_folder <- fs::path("output", "figures", "descriptive")

fs::dir_create(descriptive_figures_folder)

# ── Load data ─────────────────────────────────────────────────────────────────

hipc_merged_all_norm <- readRDS(
  fs::path(processed_data_folder, "hipc_merged_all_norm.rds")
)

# Study order fixed at preprocessing time (grouped by vaccine); reused here so
# that this figure and the Appendix A descriptives (is2_appendix_descriptives.R)
# render studies in the same order without depending on each other.
study_order <- levels(hipc_merged_all_norm$study_accession_unique)

# Vaccine name -> hex colour mapping (fixed at preprocessing time)
fill_values <- hipc_merged_all_norm %>%
  distinct(vaccine_name, vaccine_colour) %>%
  { setNames(.$vaccine_colour, .$vaccine_name) }

# =============================================================================
# BUBBLE PLOT: SAMPLE COUNTS PER STUDY x TIMEPOINT
# =============================================================================

counts <- hipc_merged_all_norm %>%
  group_by(study_accession_unique, vaccine_colour,
           time_post_last_vax, vaccine_name) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    study_accession_unique = factor(study_accession_unique, levels = study_order),
    # Order time points numerically
    time_post_last_vax = factor(
      as.character(time_post_last_vax),
      levels = sort(unique(as.numeric(as.character(time_post_last_vax)))) %>%
        as.character(),
      ordered = TRUE
    ),
    # Compress bubble sizes with a sublinear transform so that a small number
    # of large studies do not swamp the size scale; true counts are recovered
    # via the legend breaks/labels below.
    size_var = n ^ (2 / 3)
  )

size_breaks_counts <- c(10, 50, 100, 200)
size_breaks         <- size_breaks_counts ^ (2 / 3)

bubble_plot <- ggplot(counts, aes(x = time_post_last_vax, y = study_accession_unique)) +
  geom_point(
    aes(size = size_var, fill = vaccine_name),
    shape = 21, colour = "black", alpha = 0.75
  ) +
  geom_text(
    aes(label = n),
    colour = "white", size = 3.5, vjust = 0.5, show.legend = FALSE
  ) +
  scale_fill_manual(
    name   = "Vaccine",
    values = fill_values,
    guide  = guide_legend(override.aes = list(shape = 21, size = 6, colour = "black"))
  ) +
  scale_size_area(
    name     = "Count",
    max_size = 28,
    breaks   = size_breaks,
    labels   = size_breaks_counts,
    guide    = guide_legend(override.aes = list(fill = "grey80", colour = "black"))
  ) +
  scale_y_discrete(limits = rev(study_order)) +
  labs(
    x     = "Days post-vaccination",
    y     = "Study identifier",
    title = "Participants with transcriptomic samples per study across time"
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
  filename = "study_bubble_plot_sequential.pdf",
  path     = descriptive_figures_folder,
  plot     = bubble_plot,
  width    = 45, height = 40, units = "cm"
)

rm(list = ls())
