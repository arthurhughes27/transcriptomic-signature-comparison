# Master script to run all scripts in the "reanalysis" section

source(fs::path("analysis", "reanalysis", "dearseq_dgsa.R"))

source(fs::path("analysis", "reanalysis", "process_dearseq_dgsa_results.R"))

source(fs::path("analysis", "reanalysis", "qusage_dgsa.R"))