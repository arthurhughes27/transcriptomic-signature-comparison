# =============================================================================
# IS2 Dataset — Appendix A Descriptives
# =============================================================================
# Produces the Appendix A material describing the IS2 dataset (Chapter 2,
# Section 2.3.1): per-study demographic distributions and a study-level
# sample-size summary table.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(ggnewscale)
library(knitr)
library(fs)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

processed_data_folder      <- "data"
descriptive_figures_folder <- fs::path("output", "figures", "supplementary")
descriptive_tables_folder  <- fs::path("output", "tables", "supplementary")

fs::dir_create(descriptive_figures_folder)
fs::dir_create(descriptive_tables_folder)

# Timepoints highlighted on the study bubble plot below - matches
# DAYS_TO_ANALYSE in the reanalysis/specification-analysis driver scripts.
DAYS_TO_HIGHLIGHT <- c(1, 3, 7)

# ── Load data ─────────────────────────────────────────────────────────────────

hipc_merged_all_norm <- readRDS(
  fs::path(processed_data_folder, "hipc_merged_all_norm.rds")
)

# Study order fixed at preprocessing time (grouped by vaccine); kept in sync
# with is2_descriptive.R by reading the same underlying factor levels rather
# than deriving order from a rendered plot.
study_order <- levels(hipc_merged_all_norm$study_accession_unique)

# Vaccine name -> hex colour mapping (fixed at preprocessing time)
fill_values <- hipc_merged_all_norm %>%
  distinct(vaccine_name, vaccine_colour) %>%
  { setNames(.$vaccine_colour, .$vaccine_name) }

# Per-study colour and vaccine metadata, ordered to match study_order
study_colours_df <- hipc_merged_all_norm %>%
  distinct(study_accession_unique, vaccine_colour, vaccine_name) %>%
  mutate(study_accession_unique = factor(study_accession_unique, levels = study_order)) %>%
  arrange(study_accession_unique)

# Named vector: study -> vaccine hex colour (used for x-axis text colouring)
axis_colours <- setNames(
  study_colours_df$vaccine_colour,
  as.character(study_colours_df$study_accession_unique)
)

# Shared minimal theme for all x-axis plots
# axis_text_colours: named character vector of colours aligned to x-axis levels
base_theme <- function(axis_text_colours) {
  theme_minimal(base_size = 18) +
    theme(
      axis.text.x          = element_text(angle = 45, hjust = 1, size = 14,
                                          colour = axis_text_colours),
      axis.text.y          = element_text(size = 14),
      axis.title           = element_text(size = 30),
      panel.grid.major.x   = element_blank(),
      panel.grid.minor     = element_blank(),
      plot.title           = element_text(size = 35, hjust = 0.5, face = "bold"),
      legend.title         = element_text(size = 20, hjust = 0.5),
      legend.key.spacing.y = unit(0.3, "cm"),
      legend.spacing.y     = unit(1.0, "cm")
    )
}

# ggplot layer list that appends a vaccine colour legend via invisible points.
# df must contain columns: study_accession_unique, vaccine_name.
vaccine_legend_layer <- function(df) {
  list(
    new_scale_fill(),
    geom_point(
      data        = df,
      aes(x = study_accession_unique, y = -Inf, fill = vaccine_name),
      size        = 0,
      inherit.aes = FALSE
    ),
    scale_fill_manual(
      name   = "Vaccine",
      values = fill_values,
      guide  = guide_legend(override.aes = list(shape = 21, size = 6, colour = "black"))
    )
  )
}

# =============================================================================
# BUBBLE PLOT: SAMPLE COUNTS PER STUDY x TIMEPOINT
# =============================================================================
# The per-study companion to the main-text vaccine-level bubble plot
# (analysis/descriptive/is2_descriptive.R): studies on the y-axis (coloured
# by vaccine), days post-vaccination on the x-axis, bubble size
# proportional to the number of transcriptomic samples available for that
# study x timepoint combination. Moved here from the main text so the
# main-text figure gives the vaccine-level picture first. The days
# actually used in the DGSA analysis grid (DAYS_TO_HIGHLIGHT, above) are
# highlighted with a light background band, in the same per-day colours
# used throughout the specification-analysis figures
# (R/plot_helpers.R's assign_day_colours()/day_highlight_bands()).

study_bubble_counts <- hipc_merged_all_norm %>%
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

study_bubble_size_breaks_counts <- c(10, 50, 100, 200)
study_bubble_size_breaks         <- study_bubble_size_breaks_counts ^ (2 / 3)

p_study_bubble <- ggplot(study_bubble_counts, aes(x = time_post_last_vax, y = study_accession_unique)) +
  day_highlight_bands(levels(study_bubble_counts$time_post_last_vax), DAYS_TO_HIGHLIGHT) +
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
    name     = "Samples",
    max_size = 28,
    breaks   = study_bubble_size_breaks,
    labels   = study_bubble_size_breaks_counts,
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

print(p_study_bubble)

ggsave(
  filename = "study_bubble_plot_sequential.pdf",
  path     = descriptive_figures_folder,
  plot     = p_study_bubble,
  width    = 45, height = 40, units = "cm"
)

# =============================================================================
# TABLE: STUDY-LEVEL SAMPLE SIZES
# =============================================================================
# One row per study: vaccine identity, number of unique participants,
# number of transcriptomic samples, and the post-vaccination timepoints
# (days) at which samples were collected.

study_sample_sizes <- hipc_merged_all_norm %>%
  mutate(study_accession_unique = factor(study_accession_unique, levels = study_order)) %>%
  group_by(study_accession_unique) %>%
  summarise(
    pathogen        = dplyr::first(pathogen),
    vaccine_type     = dplyr::first(vaccine_type),
    n_participants   = n_distinct(participant_id),
    n_samples        = n(),
    timepoints_days  = paste(sort(unique(time_post_last_vax)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(study_accession_unique) %>%
  rename(
    `Study` = study_accession_unique,
    `Pathogen` = pathogen,
    `Vaccine type` = vaccine_type,
    `N participants` = n_participants,
    `N samples` = n_samples,
    `Timepoints (days post-vax)` = timepoints_days
  )

# Save as CSV (for records / further processing)
write.csv(
  study_sample_sizes,
  file = fs::path(descriptive_tables_folder, "study_sample_sizes.csv"),
  row.names = FALSE
)

# Save as a LaTeX table ready to \input{} into Appendix A
sample_sizes_tex <- knitr::kable(
  study_sample_sizes,
  format    = "latex",
  booktabs  = TRUE,
  longtable = TRUE,
  caption   = "Study-level sample sizes in the IS2 dataset: number of unique participants, number of transcriptomic samples, and post-vaccination timepoints (days) sampled, per study.",
  label     = "tab:study-sample-sizes"
)

writeLines(
  sample_sizes_tex,
  con = fs::path(descriptive_tables_folder, "study_sample_sizes.tex")
)

# =============================================================================
# DEMOGRAPHIC DISTRIBUTIONS
# =============================================================================

# ── Age (continuous -> violin) ───────────────────────────────────────────────

study_fill <- setNames(
  study_colours_df$vaccine_colour,
  as.character(study_colours_df$study_accession_unique)
)

age_data <- hipc_merged_all_norm %>%
  filter(!is.na(age_imputed)) %>%
  mutate(study_accession_unique = factor(study_accession_unique, levels = study_order))

p_age <- ggplot(age_data,
             aes(x = study_accession_unique, y = age_imputed,
                 fill = study_accession_unique)) +
  geom_violin(
    trim = TRUE, scale = "width", colour = "black",
    linewidth = 0.4, alpha = 0.6, show.legend = FALSE
  ) +
  stat_summary(
    fun = median, geom = "point",
    shape = 21, size = 2.5, fill = "white", colour = "black",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = study_fill) +
  vaccine_legend_layer(study_colours_df) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x     = "Study identifier",
    y     = "Age (years)",
    title = "Age distribution per study"
  ) +
  base_theme(axis_colours[levels(age_data$study_accession_unique)])

print(p_age)

# ── Gender (categorical -> stacked bar) ──────────────────────────────────────

gender_levels  <- c("Male", "Female", "Unknown")

# Muted, mid-saturation palette matching the rest of the repo's figures
# (R/plot_helpers.R's default_condition_colors()/default_aggregate_colors()),
# rather than the previous high-saturation ColorBrewer picks.
gender_colours <- c(Male = "#5B84B1", Female = "#D98880", Unknown = "#A9A9A9")

# ── Race (categorical -> stacked bar) ────────────────────────────────────────

race_levels <- c(
  "American Indian or Alaska Native", "Asian",
  "Black or African American", "White", "Other", "Unknown"
)
# Deliberately disjoint from gender_colours (no shared hues, including for
# "Unknown" beyond the shared neutral grey) so the Gender and Race
# legends/panels remain visually distinct when compared side by side. Same
# muted palette family as gender_colours above.
race_colours <- c(
  "American Indian or Alaska Native" = "#7D9D9C",
  "Asian"                            = "#8E7CC3",
  "Black or African American"        = "#C97064",
  "White"                            = "#8FAE68",
  "Other"                            = "#D9A441",
  "Unknown"                          = "#A9A9A9"
)

# Builds a stacked-bar panel of percentage composition per study for a single
# categorical covariate (e.g. gender or race), with its own vaccine legend
# (via vaccine_legend_layer()) alongside the covariate's own fill legend.
make_categorical_panel <- function(fill_col, factor_levels, colours, legend_title,
                                   plot_title) {

  plot_data <- hipc_merged_all_norm %>%
    mutate(
      !!fill_col := factor(
        if_else(.data[[fill_col]] %in% factor_levels, .data[[fill_col]], "Unknown"),
        levels = factor_levels
      ),
      study_accession_unique = factor(study_accession_unique, levels = study_order)
    ) %>%
    group_by(study_accession_unique, .data[[fill_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(study_accession_unique) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()

  ggplot(plot_data,
         aes(x = study_accession_unique, y = pct, fill = .data[[fill_col]])) +
    geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
    scale_fill_manual(
      name   = legend_title,
      values = colours,
      guide  = guide_legend(override.aes = list(colour = "black", linewidth = 0.3))
    ) +
    vaccine_legend_layer(study_colours_df) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(0, 100), expand = c(0, 0)
    ) +
    labs(
      x     = "Study identifier",
      y     = "Percentage of participants",
      title = plot_title
    ) +
    base_theme(axis_colours[levels(plot_data$study_accession_unique)])
}

p_gender <- make_categorical_panel(
  "gender", gender_levels, gender_colours, "Gender", "Gender distribution per study"
)
p_race <- make_categorical_panel(
  "race", race_levels, race_colours, "Race", "Race distribution per study"
)

print(p_gender)
print(p_race)

# Save individual demographic plots
purrr::walk2(
  list(p_age,                         p_gender,                         p_race),
  c("study_age_distribution.pdf",    "study_gender_distribution.pdf", "study_race_distribution.pdf"),
  ~ ggsave(filename = .y, path = descriptive_figures_folder,
           plot = .x, width = 45, height = 25, units = "cm")
)

rm(list = ls())
