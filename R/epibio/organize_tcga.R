
library(RnBeads)

source("config.R")
source("common.R")


dir.create(generated_TCGA_folder, recursive=TRUE, showWarnings=FALSE)

tcga_folder_name <- 'UCS'
tcga_data_folder <- file.path(external_disk_data_path, 'TCGA', tcga_folder_name, 'DNA_Methylation', 'JHU_USC__HumanMethylation450')
idat_folder <- file.path(tcga_data_folder, 'Level_1')

# targets
samples_filename <- file.path(tcga_data_folder, 'UCS_sample_annotation.txt')
targets <- read.table(samples_filename, sep='\t', header=TRUE, na.strings=c("NA", "0"), 
                quote="\"", stringsAsFactors=FALSE)
targets$barcode <- targets$Array.Data.File

targets <- head(targets, 2) # XXX
data.source <-list(idat_folder, targets)
rnb.set <- rnb.execute.import(data.source=data.source, data.type="infinium.idat.dir")
betas.table <- process_rnb_set_to_betas(rnb.set, FALSE)

name <- paste0(levels(factor(targets$histological_type)), ".", levels(factor(targets$sample_type)))
write_beta_values_table(generated_TCGA_folder, tcga_folder_name, name, betas.table)
