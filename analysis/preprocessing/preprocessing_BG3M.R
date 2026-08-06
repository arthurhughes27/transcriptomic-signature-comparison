# File to pre-process the Blood Gen 3 Modules for downstream use
## This includes removing empty genesets and genesets with no genes in our data

# Load necessary packages
library(fs)
library(dplyr)
library(stringr)
library(GSA)
library(readr)
library(readxl)

# Specify folder within folder root where the raw data lives
raw_data_folder = "data-raw"

# Specify folder within folder root where the processed data lives
processed_data_folder = "data"

# Use fs::path() to specify the data paths robustly
p_load_BG3M_raw <- fs::path(raw_data_folder, "Suppl_File_1_BIOINF.xls")
# Use the already-processed young_noNorm_expr (the analysis-ready expression
# data feeding hipc_merged_young_noNorm.rds), not the raw all_noNorm_eset, so
# genesets are filtered against the exact set of genes surviving
# preprocessing_expression.R's per-study NA filtering. Filtering against the
# raw, unfiltered gene list let genesets with zero genes actually present in
# the analysis data slip through, silently producing NA p-values in
# downstream DGSA. Requires preprocessing_expression.R to have already run
# (as it does in preprocessing_master.R's ordering).
p_load_expression <- fs::path(processed_data_folder, "young_noNorm_expr.rds")

# Load data
BG3M_raw = read_excel(p_load_BG3M_raw)
if (!fs::file_exists(p_load_expression)) {
  stop(
    "Processed expression data not found at ", p_load_expression, ".\n",
    "Run preprocessing_expression.R first (preprocessing_master.R already does this in the right order)."
  )
}
young_noNorm_expr <- readRDS(p_load_expression)

# First we must preprocess the BG3M_raw raw data into a .gmt format

# Member genes
genesets <- strsplit(BG3M_raw$`Member genes`, ",\\s*")

# Trim spaces and remove duplicates/empties
genesets <- lapply(genesets, function(x)
  unique(trimws(x[x != ""])))

BG3M <- list(
  genesets = genesets,
  geneset.names = BG3M_raw$ID,
  geneset.descriptions = BG3M_raw$`Module title`
)

# Gene names are every young_noNorm_expr column other than the two leading
# identifier columns; already lowercase (see preprocessing_expression.R),
# but re-applied here for robustness
gene_names = young_noNorm_expr %>%
  select(-participant_id, -study_time_collected) %>%
  colnames() %>%
  tolower()

# Now in the geneset data - make the geneset names lowercase
BG3M[["genesets"]] = lapply(
  BG3M[["genesets"]],
  FUN = function(x) {
    tolower(x)
  }
)

# Remove any empty elements
BG3M[["genesets"]]  = lapply(
  BG3M[["genesets"]],
  FUN = function(x) {
    x[x != ""]
  }
)

# Use only genesets which contain at least one of the commonly aligned genes in IS2
## Find the indices of the genesets containing at least one relevant gene
gs_membership <-
  sapply(BG3M[["genesets"]], function(gs)
    any(gs %in% gene_names))

# Filter the BG3M object to remove genesets without IS2 member
BG3M[["genesets"]] <- BG3M[["genesets"]][gs_membership]

BG3M[["geneset.names"]] <- BG3M[["geneset.names"]][gs_membership]

BG3M[["geneset.descriptions"]] <-
  BG3M[["geneset.descriptions"]][gs_membership]

# Rename the column names.descriptions to names
BG3M[["geneset.names.descriptions"]] = paste0(BG3M[["geneset.descriptions"]], " (", BG3M[["geneset.names"]], ")")

names(BG3M)[names(BG3M) == "Cluster"] = "geneset.aggregates"

names(BG3M[["genesets"]]) = NULL

# Remove TBA genesets
BG3M_filter <- !grepl("TBD", BG3M[["geneset.names.descriptions"]])
BG3M[["geneset.names.descriptions"]] <- BG3M[["geneset.names.descriptions"]][BG3M_filter]
BG3M[["geneset.names"]] <- BG3M[["geneset.names"]][BG3M_filter]
BG3M[["geneset.descriptions"]] <- BG3M[["geneset.descriptions"]][BG3M_filter]
BG3M[["geneset.aggregates"]] <- BG3M[["geneset.aggregates"]][BG3M_filter]
BG3M[["genesets"]] <- BG3M[["genesets"]][BG3M_filter]

# Save the processed geneset object

# Use fs::path() to specify the data path robustly
p_save <- fs::path(processed_data_folder, "BG3M_processed.rds")

# Save dataframe
saveRDS(BG3M, file = p_save)

rm(list = ls())
