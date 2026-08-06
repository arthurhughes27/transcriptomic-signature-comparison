# Sources every R/*.R file (excluding this loader itself) so that all
# modular analysis functions become available with a single call, e.g. from
# an analysis/ driver script:
#
#   source(fs::path("R", "load_all.R"))
#
# Assumes the working directory is the project root, consistent with every
# other script in this repository.

.load_all_r_dir <- fs::path("R")

.load_all_files <- sort(fs::dir_ls(.load_all_r_dir, regexp = "\\.R$"))
.load_all_files <- .load_all_files[fs::path_file(.load_all_files) != "load_all.R"]

invisible(lapply(.load_all_files, source))

rm(.load_all_r_dir, .load_all_files)
