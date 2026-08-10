# Master script to run all scripts in the "reanalysis" section
#
# All four scripts read baseline results from the specification-analysis
# raw output (R/baseline_results.R::load_baseline_results_from_raw()), so
# analysis/specification_analysis/01_build_specification_grid.R and
# 02_run_raw_specifications.R must already have been run first. The
# earlier standalone dearseq_dgsa.R / qusage_dgsa.R baseline-only driver
# scripts (and their process_*_dgsa_results.R post-processors) have been
# removed as redundant: 02_run_raw_specifications.R's raw specification
# grid already includes both methods' baseline specifications.

# Stage 1: dearseq vs QuSAGE baseline comparison
source(fs::path("analysis", "reanalysis", "01_dgsa_comparison_example.R"))

# Circos plots comparing the two methods
source(fs::path("analysis", "reanalysis", "02_plot_circos.R"))

# Heatmap comparison of the two methods
source(fs::path("analysis", "reanalysis", "03_plot_heatmap_comparison.R"))

# Baseline significance comparison (percentage significant, dearseq vs QuSAGE)
source(fs::path("analysis", "reanalysis", "04_plot_baseline_significance_comparison.R"))
