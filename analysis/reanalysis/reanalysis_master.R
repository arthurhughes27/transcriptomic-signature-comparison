# Master script to run all scripts in the "reanalysis" section

# Baseline dearseq specification
source(fs::path("analysis", "reanalysis", "dearseq_dgsa.R"))
source(fs::path("analysis", "reanalysis", "process_dearseq_dgsa_results.R"))

# Baseline QuSAGE specification
source(fs::path("analysis", "reanalysis", "qusage_dgsa.R"))
source(fs::path("analysis", "reanalysis", "process_qusage_dgsa_results.R"))

# Stage 1: dearseq vs QuSAGE baseline comparison
source(fs::path("analysis", "reanalysis", "dgsa_comparison_example.R"))

# Circos plots comparing the two methods
source(fs::path("analysis", "reanalysis", "plot_circos.R"))
