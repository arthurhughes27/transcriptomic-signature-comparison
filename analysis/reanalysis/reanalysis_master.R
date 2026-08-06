# Master script to run all scripts in the "reanalysis" section

source(fs::path("analysis", "reanalysis", "dearseq_dgsa.R"))

source(fs::path("analysis", "reanalysis", "process_dearseq_dgsa_results.R"))

# QuSAGE driver script (built on R/dgsa_qusage.R's run_qusage_comparison())
# is pending as part of the Stage 1 baseline-comparison driver.