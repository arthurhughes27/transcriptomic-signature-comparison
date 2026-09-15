# Master script to run all supplementary analyses in order

# Per-study companion to the main-text vaccine-level bubble plot (Section 2.3.1)
source(fs::path("analysis", "supplementary", "is2_descriptive_supplementary.R"))

# Gene-set-level companion to the main-text aggregate robustness heatmap
source(fs::path("analysis", "supplementary", "specification_heatmap_genesets_supplementary.R"))

# Deduplicated specification-count table, per vaccine x timepoint comparison
source(fs::path("analysis", "supplementary", "effective_specification_totals_table.R"))
