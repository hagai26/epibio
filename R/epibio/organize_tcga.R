
library(RnBeads)

source("config.R")
source("common.R")

work_on_targets <- function(targets, idat_folder) {
  type <- targets$histological_type[[1]]
  study <- targets$sample_type[[1]]  
  targets <- head(targets, 5) # XXX
  data.source <-list(idat_folder, targets)
  rnb.set <- rnb.execute.import(data.source=data.source, data.type="infinium.idat.dir")
  betas.table <- process_rnb_set_to_betas(rnb.set, FALSE)
  write_beta_values_table(generated_TCGA_folder, tcga_folder_name, study, type, betas.table)  
}

work_on_tcga_folder <- function(tcga_inside_folder) {
  # folder structure looks like this:
  # TCGA\UCS\DNA_Methylation\JHU_USC__HumanMethylation450\Level_1
  tcga_inside_folder <- list.files(tcga_inside_folder, full.names = TRUE)
  tcga_inside_folder <- list.files(tcga_inside_folder, full.names = TRUE)
  samples_filename <- file.path(tcga_inside_folder, 'UCS_sample_annotation.txt')
  idat_folder <- file.path(tcga_inside_folder, 'Level_1')
  
  targets <- read.table(samples_filename, sep='\t', header=TRUE, na.strings=c("NA", "0"), 
                        quote="\"", stringsAsFactors=FALSE)
  targets$barcode <- targets$Array.Data.File
  splited_targets <- split(targets, list(targets$histological_type, targets$sample_type), drop=TRUE)
  ret <- lapply(splited_targets, FUN=work_on_targets, idat_folder)
}


dir.create(generated_TCGA_folder, recursive=TRUE, showWarnings=FALSE)

tcga_folder <- file.path(external_disk_data_path, 'TCGA')
tcga_inside_folders <- list.files(tcga_folder, full.names = TRUE)
lapply(tcga_inside_folders, FUN=work_on_tcga_folder)
