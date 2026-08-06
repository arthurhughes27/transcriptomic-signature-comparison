# =============================================================================
# IS2 Dataset — Appendix A Descriptives
# =============================================================================
# Produces the Appendix A material describing the IS2 dataset (Chapter 2,
# Section 2.3.1): per-study covariate distributions (demographics, immune
# response assay availability/distributions) and a study-level sample-size
# summary table.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(tibble)
library(stringr)
library(ggnewscale)
library(cowplot)
library(knitr)
library(fs)

# ── Paths ─────────────────────────────────────────────────────────────────────

processed_data_folder      <- "data"
descriptive_figures_folder <- fs::path("output", "figures", "descriptive")
descriptive_tables_folder  <- fs::path("output", "tables", "descriptive")

fs::dir_create(descriptive_figures_folder)
fs::dir_create(descriptive_tables_folder)

# ── Load data ─────────────────────────────────────────────────────────────────

hipc_merged_all_norm <- readRDS(
  fs::path(processed_data_folder, "hipc_merged_all_norm.rds")
)

# Study order fixed at preprocessing time (grouped by vaccine); kept in sync
# with is2_bubble_plot.R by reading the same underlying factor levels rather
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
# IMMUNE RESPONSE ASSAY AVAILABILITY
# =============================================================================

assay_cols <- c(
  "immResp_MFC_nAb_pre_value",
  "immResp_MFC_elisa_pre_value",
  "immResp_MFC_hai_pre_value"
)
assay_labels <- c(
  immResp_MFC_nAb_pre_value   = "nAb",
  immResp_MFC_elisa_pre_value = "ELISA",
  immResp_MFC_hai_pre_value   = "HAI"
)

assay_avail <- hipc_merged_all_norm %>%
  select(study_accession_unique, all_of(assay_cols)) %>%
  pivot_longer(cols = all_of(assay_cols), names_to = "assay", values_to = "value") %>%
  group_by(study_accession_unique, assay) %>%
  summarise(available = any(!is.na(value)), .groups = "drop") %>%
  mutate(
    assay                  = factor(assay_labels[assay], levels = rev(assay_labels)),
    study_accession_unique = factor(study_accession_unique, levels = study_order)
  )

p_assay_avail <- ggplot(assay_avail, aes(x = study_accession_unique, y = assay, fill = available)) +
  geom_tile(colour = "grey70", linewidth = 0.4) +
  scale_fill_manual(
    name   = "Assay available",
    values = c("TRUE" = "#2ECC71", "FALSE" = "white"),
    labels = c("TRUE" = "At least one sample", "FALSE" = "No samples"),
    guide  = guide_legend(override.aes = list(colour = "grey70", linewidth = 0.4))
  ) +
  vaccine_legend_layer(study_colours_df) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x     = "Study identifier",
    y     = "Assay",
    title = "Immune response assay availability per study"
  ) +
  base_theme(axis_colours[levels(assay_avail$study_accession_unique)]) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 25, hjust = 0.5)
  )

print(p_assay_avail)

ggsave(
  filename = "study_immResp_availability.pdf",
  path     = descriptive_figures_folder,
  plot     = p_assay_avail,
  width    = 45, height = 20, units = "cm"
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
gender_colours <- c(Male = "#4E79A7", Female = "#F28E2B", Unknown = "#B0B0B0")

# ── Race (categorical -> stacked bar) ────────────────────────────────────────

race_levels <- c(
  "American Indian or Alaska Native", "Asian",
  "Black or African American", "White", "Other", "Unknown"
)
# Deliberately disjoint from gender_colours (no shared hues, including for
# "Unknown") so the Gender and Race legends/panels remain visually distinct
# when shown together in the combined demographics figure.
race_colours <- c(
  "American Indian or Alaska Native" = "#1B9E77",
  "Asian"                            = "#7570B3",
  "Black or African American"        = "#E7298A",
  "White"                            = "#66A61E",
  "Other"                            = "#E6AB02",
  "Unknown"                          = "#A6761D"
)

# Builds a stacked-bar panel of percentage composition per study for a single
# categorical covariate (e.g. gender or race). When `vaccine_legend` is TRUE,
# an additional (invisible) vaccine fill scale is layered on via ggnewscale so
# the plot carries its own vaccine legend when used standalone; this is left
# FALSE for panels destined for the combined figure below, where the vaccine
# legend is shown once (from the age panel) instead of once per panel, and
# keeping a single fill scale per plot makes legend extraction unambiguous.
make_categorical_panel <- function(fill_col, factor_levels, colours, legend_title,
                                   plot_title, vaccine_legend = TRUE) {

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

  p <- ggplot(plot_data,
              aes(x = study_accession_unique, y = pct, fill = .data[[fill_col]])) +
    geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
    scale_fill_manual(
      name   = legend_title,
      values = colours,
      guide  = guide_legend(override.aes = list(colour = "black", linewidth = 0.3))
    )

  if (vaccine_legend) p <- p + vaccine_legend_layer(study_colours_df)

  p +
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
           plot = .x, width = 45, height = 20, units = "cm")
)

# =============================================================================
# COMBINED DEMOGRAPHICS FIGURE (age / gender / race stacked)
# =============================================================================
# Each panel below carries exactly one fill scale (age: vaccine; gender:
# gender; race: race), so legends can be extracted unambiguously with
# cowplot::get_legend() and each shown exactly once in the shared legend
# column, rather than the previous approach of layering a second "Vaccine"
# fill scale onto every panel via ggnewscale, which made get_legend()/guides()
# extraction ambiguous and produced garbled combined legends.

# Bare gender/race panels for the combined figure only: single fill scale,
# no vaccine legend layer (the age panel's vaccine legend already covers this,
# since all three panels share the same vaccine-coloured x-axis text).
p_gender_combined <- make_categorical_panel(
  "gender", gender_levels, gender_colours, "Gender", "Gender distribution per study",
  vaccine_legend = FALSE
)
p_race_combined <- make_categorical_panel(
  "race", race_levels, race_colours, "Race", "Race distribution per study",
  vaccine_legend = FALSE
)

# Strip legends and x-axis labels from the upper two panels
strip_x <- theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank(),
  legend.position = "none"
)

legend_theme <- theme(
  legend.position       = "right",
  legend.title          = element_text(size = 20, hjust = 0.5),
  legend.text           = element_text(size = 14),
  legend.key.spacing.y  = unit(0.3, "cm")
)

p_age_bare    <- p_age             + strip_x
p_gender_bare <- p_gender_combined + strip_x
p_race_bare   <- p_race_combined   + theme(legend.position = "none")

# Extract legends separately for the legend column (each source plot has a
# single fill scale, so each call unambiguously grabs the intended legend)
legend_vaccine <- get_legend(p_age             + legend_theme)
legend_gender  <- get_legend(p_gender_combined + legend_theme)
legend_race    <- get_legend(p_race_combined   + legend_theme)

legend_col <- plot_grid(
  legend_vaccine, legend_gender, legend_race,
  ncol        = 1,
  rel_heights = c(length(fill_values), length(gender_levels), length(race_levels))
)

panels <- plot_grid(
  p_age_bare, p_gender_bare, p_race_bare,
  ncol        = 1,
  align       = "v",
  axis        = "lr",
  rel_heights = c(1.2, 1, 1),
  labels      = c("A", "B", "C"),
  label_size  = 28
)

p_demographics_combined <- ggdraw(
  plot_grid(panels, legend_col, ncol = 2, rel_widths = c(4, 1))
) 

print(p_demographics_combined)

ggsave(
  filename = "study_demographics_combined.pdf",
  path     = descriptive_figures_folder,
  plot     = p_demographics_combined,
  width    = 45, height = 60, units = "cm"
)

# =============================================================================
# IMMUNE ASSAY DISTRIBUTIONS (pre vs post vaccination)
# =============================================================================

assay_display      <- c(nAb = "nAb", hai = "HAI", elisa = "ELISA")
timepoint_colours  <- c("Pre-vaccination" = "#6BAED6", "Post-vaccination" = "#E6550D")

# Builds a dodged violin plot for one assay comparing pre- and post-vaccination
# log-transformed values, restricted to studies with at least one non-NA value.
make_assay_violin <- function(assay_name) {

  pre_col  <- paste0("immResp_MFC_", assay_name, "_pre_value")
  post_col <- paste0("immResp_MFC_", assay_name, "_post_value")

  plot_data <- hipc_merged_all_norm %>%
    select(study_accession_unique, all_of(c(pre_col, post_col))) %>%
    pivot_longer(
      cols      = all_of(c(pre_col, post_col)),
      names_to  = "timepoint",
      values_to = "value"
    ) %>%
    filter(!is.na(value)) %>%
    mutate(
      value     = log(value),
      timepoint = factor(
        if_else(str_detect(timepoint, "_pre_"), "Pre-vaccination", "Post-vaccination"),
        levels = c("Pre-vaccination", "Post-vaccination")
      ),
      # Preserve global study ordering, keeping only studies present in this assay
      study_accession_unique = factor(
        study_accession_unique,
        levels = intersect(study_order, unique(study_accession_unique))
      )
    )

  studies_present   <- levels(plot_data$study_accession_unique)
  study_colours_sub <- study_colours_df %>%
    filter(study_accession_unique %in% studies_present) %>%
    mutate(study_accession_unique = factor(study_accession_unique, levels = studies_present))
  axis_colours_sub  <- axis_colours[studies_present]

  dw <- 0.85

  ggplot(plot_data,
         aes(x = study_accession_unique, y = value, fill = timepoint)) +
    geom_violin(
      aes(group = interaction(study_accession_unique, timepoint)),
      position  = position_dodge(width = dw),
      trim      = TRUE, scale = "width",
      colour    = "black", linewidth = 0.35, alpha = 0.85
    ) +
    stat_summary(
      aes(group = interaction(study_accession_unique, timepoint)),
      fun      = median, geom = "point",
      shape    = 21, size = 2, fill = "white", colour = "black",
      position = position_dodge(width = dw),
      show.legend = FALSE
    ) +
    scale_fill_manual(
      name   = "Timepoint",
      values = timepoint_colours,
      guide  = guide_legend(override.aes = list(colour = "black", linewidth = 0.35))
    ) +
    vaccine_legend_layer(study_colours_sub) +
    scale_x_discrete(drop = FALSE) +
    labs(
      x     = "Study identifier",
      y     = paste0("log(", assay_display[assay_name], " value)"),
      title = paste0(assay_display[assay_name],
                     " assay distribution per study (pre vs post vaccination)")
    ) +
    base_theme(axis_colours_sub[studies_present])
}

p_nab   <- make_assay_violin("nAb")
p_hai   <- make_assay_violin("hai")
p_elisa <- make_assay_violin("elisa")

print(p_nab)
print(p_hai)
print(p_elisa)

# Save assay plots (hai is slightly wider to accommodate more studies)
ggsave(
  filename = "study_nAb_distribution.pdf",
  path =  descriptive_figures_folder,
  plot = p_nab,
  width = 45,
  height = 20,
  units = "cm"
)
ggsave(
  filename = "study_hai_distribution.pdf",
  path =  descriptive_figures_folder,
  plot = p_hai,
  width = 50,
  height = 20,
  units = "cm"
)
ggsave(
  filename = "study_elisa_distribution.pdf",
  path = descriptive_figures_folder,
  plot = p_elisa,
  width = 45,
  height = 20,
  units = "cm"
)

rm(list = ls())
