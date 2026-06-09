# R Script to describe the different studies present in the HIPC IS2 dataset

# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(tibble)

# Directory to access processed data
processed_data_folder = "data"

# Directory to store figures
descriptive_figures_folder = fs::path("output", "figures", "descriptive")

# Path to processed gene-level data
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load merged gene-level data
hipc_merged_all_norm = readRDS(p_load_expr_all_norm)


# PLOT 1: COUNTS ACROSS TIME #

# Summarise counts per vaccine x time
counts <- hipc_merged_all_norm %>%
  group_by(study_accession_unique,
           vaccine_colour,
           time_post_last_vax,
           vaccine_name) %>%
  summarise(n = n(), .groups = "drop")

# Order the time points numerically and make time_post_last_vax an ordered factor
time_levels <- counts %>%
  distinct(time_post_last_vax) %>%
  arrange(as.numeric(time_post_last_vax)) %>%
  pull(time_post_last_vax) %>%
  as.character()   # factor levels must be character

# Order the counts by the study time
counts <- counts %>%
  mutate(time_post_last_vax = factor(
    as.character(time_post_last_vax),
    levels = time_levels,
    ordered = TRUE
  ))

# Since the range of the number of samples is large,
# we make the bubble size proportional to a nonlinear transformation of the count (to the power 2/3)
counts <- counts %>%
  mutate(size_var = n^{
    2 / 3
  })

# Extract study, vaccine and colour info for legend
vaccine_map <- hipc_merged_all_norm %>%
  distinct(study_accession_unique, vaccine_name, vaccine_colour)

fill_values <- vaccine_map %>%
  distinct(vaccine_name, vaccine_colour) %>%
  {
    setNames(.$vaccine_colour, .$vaccine_name)
  }

# Choose discrete values for counts to show in legend
size_breaks_counts <- c(10, 50, 100, 200)

# Convert them to the scale used in the plot
size_breaks <- (size_breaks_counts)^{
  2 / 3
}

# Make the plot
p1 <- ggplot(counts, aes(x = time_post_last_vax, y = study_accession_unique)) +
  geom_point(
    aes(size = size_var, fill = vaccine_name),
    shape = 21,
    colour = "black",
    alpha = 0.75,
    show.legend = TRUE
  ) +
  geom_text(
    aes(label = n),
    colour = "white",
    size = 3.5,
    vjust = 0.5,
    show.legend = FALSE
  ) +
  # Fill legend: vaccine names with their hex colours
  scale_fill_manual(
    name = "Vaccine",
    values = fill_values,
    guide = guide_legend(override.aes = list(
      shape = 21,
      size = 6,
      colour = "black"
    ))
  ) +
  # Size legend: show 10 / 50 / 100 / 200 as legend entries (we pass their sqrt values)
  scale_size_area(
    name = "Count",
    max_size = 28,
    breaks = size_breaks,
    labels = size_breaks_counts,
    guide = guide_legend(override.aes = list(
      fill = "grey80", colour = "black"
    ))
  ) +
  labs(
    x = "Days post-vaccination",
    y = "Study identifier",
    title = "Participants with transcriptomic samples per study across time"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 35, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 15, hjust = 0.5),
    legend.title = element_text(size = 20, hjust = 0.5), 
    legend.key.spacing.y = unit(0.3, 'cm'),
    legend.spacing.y = unit(1.0, 'cm')
  )

print(p1)

ggsave(
  filename = "study_bubble_plot_sequential.pdf",
  path = descriptive_figures_folder,
  plot = p1,
  width = 45,
  height = 40,
  units = "cm"
)

# PLOT 2: ASSAY AVAILABILITY ACROSS STUDIES #

# PLOT 2: ASSAY AVAILABILITY ACROSS STUDIES #

# Define assay names and labels
assay_cols <- c("immResp_MFC_nAb_pre_value",
                "immResp_MFC_elisa_pre_value",
                "immResp_MFC_hai_pre_value")

assay_labels <- c(
  "immResp_MFC_nAb_pre_value"   = "nAb",
  "immResp_MFC_elisa_pre_value" = "ELISA",
  "immResp_MFC_hai_pre_value"   = "HAI"
)

# Use the exact study order from p1's rendered y-axis
study_order <- ggplot_build(p1)$layout$panel_params[[1]]$y$breaks

# Build availability data frame: for each study x assay, check for any non-NA value
assay_avail <- hipc_merged_all_norm %>%
  select(study_accession_unique, all_of(assay_cols)) %>%
  pivot_longer(cols = all_of(assay_cols),
               names_to = "assay",
               values_to = "value") %>%
  group_by(study_accession_unique, assay) %>%
  summarise(available = any(!is.na(value)), .groups = "drop") %>%
  mutate(
    assay = factor(assay_labels[assay], levels = rev(assay_labels)),
    study_accession_unique = factor(study_accession_unique, levels = study_order)
  )

# Extract per-study vaccine colour for x-axis text colouring
study_colours <- hipc_merged_all_norm %>%
  distinct(study_accession_unique, vaccine_colour, vaccine_name) %>%
  mutate(study_accession_unique = factor(study_accession_unique, levels = study_order)) %>%
  arrange(study_accession_unique)

# Build named vector: study -> colour (for axis text colouring)
axis_colours <- setNames(study_colours$vaccine_colour,
                         as.character(study_colours$study_accession_unique))

# Make the plot
p2 <- ggplot(assay_avail,
             aes(x = study_accession_unique, y = assay, fill = available)) +
  geom_tile(colour = "grey70", linewidth = 0.4) +
  scale_fill_manual(
    name = "Assay available",
    values = c("TRUE"  = "#2ECC71",
               "FALSE" = "white"),
    labels = c("TRUE"  = "At least one sample",
               "FALSE" = "No samples"),
    guide = guide_legend(
      override.aes = list(colour = "grey70", linewidth = 0.4)
    )
  ) +
  # Vaccine colour legend (same fill_values as p1)
  ggnewscale::new_scale_fill() +
  geom_point(
    data = study_colours,
    aes(x = study_accession_unique,
        y = -Inf,
        fill = vaccine_name),
    size = 0,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name = "Vaccine",
    values = fill_values,
    guide = guide_legend(
      override.aes = list(shape = 21, size = 6, colour = "black")
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = "Study identifier",
    y = "Assay",
    title = "Immune response assay availability per study"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, size = 14,
      colour = axis_colours[levels(assay_avail$study_accession_unique)]
    ),
    axis.text.y        = element_text(size = 18),
    axis.title         = element_text(size = 30),
    panel.grid         = element_blank(),
    plot.title         = element_text(size = 35, hjust = 0.5, face = "bold"),
    legend.title       = element_text(size = 25, hjust = 0.5),
    legend.key.spacing.y = unit(0.3, "cm"),
    legend.spacing.y   = unit(1.0, "cm")
  )

print(p2)

ggsave(
  filename = "study_immResp_availability.pdf",
  path = descriptive_figures_folder,
  plot = p2,
  width = 45,
  height = 20,
  units = "cm"
)

