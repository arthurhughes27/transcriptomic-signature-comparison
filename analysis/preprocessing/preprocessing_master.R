# Master preprocessing script: runs all preprocessing steps in order
# Each script reads from and writes to the data/ directory

# Preprocessing clinical data
source("analysis/preprocessing/preprocessing_clinical.R")

# Preprocessing gene expression data
source("analysis/preprocessing/preprocessing_expression.R")

# Preprocessing immune response data
source("analysis/preprocessing/preprocessing_immuneresponse.R")

# Preprocessing gene set data - BTM and BloodGen3 Modules
source("analysis/preprocessing/preprocessing_BTM.R")
source("analysis/preprocessing/preprocessing_BG3M.R")

# Merging clinical, gene expression, and immune response data
source("analysis/preprocessing/preprocessing_merging.R")
