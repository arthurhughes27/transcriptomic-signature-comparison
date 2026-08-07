# Master script to run all scripts in the "reanalysis" section

# Baseline dearseq specification
source(fs::path("analysis", "reanalysis", "dearseq_dgsa.R"))
source(fs::path("analysis", "reanalysis", "process_dearseq_dgsa_results.R"))

# Baseline QuSAGE specification
source(fs::path("analysis", "reanalysis", "qusage_dgsa.R"))
source(fs::path("analysis", "reanalysis", "process_qusage_dgsa_results.R"))

# Stage 1: dearseq vs QuSAGE baseline comparison
# Reads baseline results from the specification-analysis raw output (see
# R/baseline_results.R), so analysis/specification_analysis/01_build_specification_grid.R
# and 02_run_raw_specifications.R must already have been run.
source(fs::path("analysis", "reanalysis", "dgsa_comparison_example.R"))

# Circos plots comparing the two methods
# Also reads baseline results from the specification-analysis raw output
# (see R/baseline_results.R) - same prerequisite as above.
source(fs::path("analysis", "reanalysis", "plot_circos.R"))

# Heatmap comparison of the two methods
# Also reads baseline results from the specification-analysis raw output
# (see R/baseline_results.R) - same prerequisite as above.
source(fs::path("analysis", "reanalysis", "plot_heatmap_comparison.R"))

# Baseline significance comparison (percentage significant, dearseq vs QuSAGE)
# Also reads baseline results from the specification-analysis raw output
# (see R/baseline_results.R) - same prerequisite as above.
source(fs::path("analysis", "reanalysis", "plot_baseline_significance_comparison.R"))
