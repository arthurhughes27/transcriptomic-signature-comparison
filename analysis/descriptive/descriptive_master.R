# Master descriptive analysis script: runs all descriptive analyses in order

# Section 2.3.1 main-text figure: study x timepoint sample bubble plot
source(fs::path("analysis", "descriptive", "is2_bubble_plot.R"))

# Appendix A: covariate distributions and study-level sample size table
source(fs::path("analysis", "descriptive", "is2_appendix_descriptives.R"))