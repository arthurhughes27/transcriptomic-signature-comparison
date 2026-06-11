# =============================================================================
# HIPC IS2 Dataset — Descriptive Figures
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(tibble)
library(patchwork)
library(stringr)
library(ggnewscale)
library(cowplot)
library(fs)

# ── Paths ─────────────────────────────────────────────────────────────────────

processed_data_folder     <- "data"
descriptive_figures_folder <- fs::path("output", "figures", "descriptive")

# ── Load data ─────────────────────────────────────────────────────────────────

hipc_merged_all_norm <- readRDS(
  fs::path(processed_data_folder, "hipc_merged_all_norm.rds")
)

# =============================================================================
# SHARED AESTHETICS
# =============================================================================

# Vaccine name -> hex colour mapping (used in all fill legends)
fill_values <- hipc_merged_all_norm %>%
  distinct(vaccine_name, vaccine_colour) %>%
  { setNames(.$vaccine_colour, .$vaccine_name) }

# =============================================================================
# P1: SAMPLE COUNTS PER STUDY ACROSS TIME
# =============================================================================

counts <- hipc_merged_all_norm %>%
  group_by(study_accession_unique, vaccine_colour,
           time_post_last_vax, vaccine_name) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    # Order time points numerically
    time_post_last_vax = factor(
      as.character(time_post_last_vax),
      levels = sort(unique(as.numeric(as.character(time_post_last_vax)))) %>%
        as.character(),
      ordered = TRUE
    ),
    # Compress bubble sizes with a sublinear transform
    size_var = n ^ (2 / 3)
  )

size_breaks_counts <- c(10, 50, 100, 200)
size_breaks        <- size_breaks_counts ^ (2 / 3)

p1 <- ggplot(counts, aes(x = time_post_last_vax, y = study_accession_unique)) +
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

print(p1)

ggsave(
  filename = "study_bubble_plot_sequential.pdf",
  path     = descriptive_figures_folder,
  plot     = p1,
  width    = 45, height = 40, units = "cm"
)

# =============================================================================
# SHARED SETUP FOR P2–P9
# =============================================================================

# Study order derived from p1's rendered y-axis (bottom-to-top → left-to-right)
study_order <- rev(ggplot_build(p1)$layout$panel_params[[1]]$y$breaks)

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

# Shared minimal theme for all x-axis plots (p2–p9)
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
# P2: IMMUNE RESPONSE ASSAY AVAILABILITY
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

p2 <- ggplot(assay_avail, aes(x = study_accession_unique, y = assay, fill = available)) +
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

print(p2)

ggsave(
  filename = "study_immResp_availability.pdf",
  path     = descriptive_figures_folder,
  plot     = p2,
  width    = 45, height = 20, units = "cm"
)

# =============================================================================
# P3–P5: DEMOGRAPHIC DISTRIBUTIONS
# =============================================================================

# ── P3: Age (continuous → violin) ────────────────────────────────────────────

study_fill <- setNames(
  study_colours_df$vaccine_colour,
  as.character(study_colours_df$study_accession_unique)
)

age_data <- hipc_merged_all_norm %>%
  filter(!is.na(age_imputed)) %>%
  mutate(study_accession_unique = factor(study_accession_unique, levels = study_order))

p3 <- ggplot(age_data,
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

print(p3)

# ── P4: Gender (categorical → stacked bar) ───────────────────────────────────

gender_levels  <- c("Male", "Female", "Unknown")
gender_colours <- c(Male = "#4E79A7", Female = "#F28E2B", Unknown = "#B0B0B0")

gender_data <- hipc_merged_all_norm %>%
  mutate(
    gender = factor(
      if_else(gender %in% gender_levels, gender, "Unknown"),
      levels = gender_levels
    ),
    study_accession_unique = factor(study_accession_unique, levels = study_order)
  ) %>%
  group_by(study_accession_unique, gender) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(study_accession_unique) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p4 <- ggplot(gender_data,
             aes(x = study_accession_unique, y = pct, fill = gender)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
  scale_fill_manual(
    name   = "Gender",
    values = gender_colours,
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
    title = "Gender distribution per study"
  ) +
  base_theme(axis_colours[levels(gender_data$study_accession_unique)])

print(p4)

# ── P5: Race (categorical → stacked bar) ─────────────────────────────────────

race_levels <- c(
  "American Indian or Alaska Native", "Asian",
  "Black or African American", "White", "Other", "Unknown"
)
race_colours <- c(
  "American Indian or Alaska Native" = "#59A14F",
  "Asian"                            = "#F28E2B",
  "Black or African American"        = "#E15759",
  "White"                            = "#4E79A7",
  "Other"                            = "#B07AA1",
  "Unknown"                          = "#B0B0B0"
)

race_data <- hipc_merged_all_norm %>%
  mutate(
    race = factor(
      if_else(race %in% race_levels, race, "Unknown"),
      levels = race_levels
    ),
    study_accession_unique = factor(study_accession_unique, levels = study_order)
  ) %>%
  group_by(study_accession_unique, race) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(study_accession_unique) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p5 <- ggplot(race_data,
             aes(x = study_accession_unique, y = pct, fill = race)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
  scale_fill_manual(
    name   = "Race",
    values = race_colours,
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
    title = "Race distribution per study"
  ) +
  base_theme(axis_colours[levels(race_data$study_accession_unique)])

print(p5)

# Save individual demographic plots
purrr::walk2(
  list(p3,                           p4,                        p5),
  c("study_age_distribution.pdf",    "study_gender_distribution.pdf", "study_race_distribution.pdf"),
  ~ ggsave(filename = .y, path = descriptive_figures_folder,
           plot = .x, width = 45, height = 20, units = "cm")
)

# =============================================================================
# P6: COMBINED DEMOGRAPHICS FIGURE (p3 / p4 / p5 stacked)
# =============================================================================

# Strip legends and x-axis labels from the upper two panels
strip_x <- theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank(),
  legend.position = "none"
)

p3_bare <- p3 + strip_x
p4_bare <- p4 + strip_x
p5_bare <- p5 + theme(legend.position = "none")

# Extract legends separately for the legend column
legend_vaccine <- get_legend(
  p3 + theme(legend.position = "right",
             legend.title = element_text(size = 20, hjust = 0.5),
             legend.text  = element_text(size = 14),
             legend.key.spacing.y = unit(0.3, "cm"))
)
legend_gender <- get_legend(
  p4 +
    guides(fill = guide_legend(title = "Gender",
                               override.aes = list(colour = "black", linewidth = 0.3))) +
    theme(legend.position = "right",
          legend.title = element_text(size = 20, hjust = 0.5),
          legend.text  = element_text(size = 14),
          legend.key.spacing.y = unit(0.3, "cm"))
)
legend_race <- get_legend(
  p5 +
    guides(fill = guide_legend(title = "Race",
                               override.aes = list(colour = "black", linewidth = 0.3))) +
    theme(legend.position = "right",
          legend.title = element_text(size = 20, hjust = 0.5),
          legend.text  = element_text(size = 14),
          legend.key.spacing.y = unit(0.3, "cm"))
)

legend_col <- plot_grid(
  legend_vaccine, legend_gender, legend_race,
  ncol        = 1,
  rel_heights = c(length(fill_values), length(gender_levels), length(race_levels))
)

panels <- plot_grid(
  p3_bare, p4_bare, p5_bare,
  ncol        = 1,
  align       = "v",
  axis        = "lr",
  rel_heights = c(1.2, 1, 1),
  labels      = c("A", "B", "C"),
  label_size  = 28
)

p6 <- ggdraw(
  plot_grid(panels, legend_col, ncol = 2, rel_widths = c(4, 1))
) +
  draw_label(
    "Participant demographics per study",
    x = 0.5, y = 0.995, vjust = 1, hjust = 0.5,
    fontface = "bold", size = 35
  )

print(p6)

ggsave(
  filename = "study_demographics_combined.pdf",
  path     = descriptive_figures_folder,
  plot     = p6,
  width    = 45, height = 55, units = "cm"
)

# =============================================================================
# P7–P9: IMMUNE ASSAY DISTRIBUTIONS (pre vs post vaccination)
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

p7 <- make_assay_violin("nAb")
p8 <- make_assay_violin("hai")
p9 <- make_assay_violin("elisa")

print(p7)
print(p8)
print(p9)

# Save assay plots (hai is slightly wider to accommodate more studies)
ggsave(
  filename = "study_nAb_distribution.pdf",
  path =  descriptive_figures_folder,
  plot = p7,
  width = 45,
  height = 20,
  units = "cm"
)
ggsave(
  filename = "study_hai_distribution.pdf",
  path =  descriptive_figures_folder,
  plot = p8,
  width = 50,
  height = 20,
  units = "cm"
)
ggsave(
  filename = "study_elisa_distribution.pdf",
  path = descriptive_figures_folder,
  plot = p9,
  width = 45,
  height = 20,
  units = "cm"
)

rm(list = ls())