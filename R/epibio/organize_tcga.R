
library(RnBeads)

source("config.R")
source("common.R")


dir.create(generated_TCGA_folder, recursive=TRUE, showWarnings=FALSE)
tcga_data_folder <- file.path(external_disk_data_path, 'TCGA')

