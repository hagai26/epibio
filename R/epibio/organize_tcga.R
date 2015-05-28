
library(RnBeads)

source("config.R")
source("common.R")
source("tcga_utils.R")
args <- commandArgs(trailingOnly = TRUE)

run_organize_tcga <- function() {
	dir.create(generated_TCGA_folder, recursive=TRUE, showWarnings=FALSE)

	tcga_folder <- file.path(external_disk_data_path, 'TCGA')
	stopifnot(file.exists(tcga_folder))
	tcga_inside_folders <- list.dirs(tcga_folder)
	ignore_list <- c("Clinical", "download_2", "download_3", "download_4")
	tcga_inside_folders <- tcga_inside_folders[!(tcga_inside_folders %in% ignore_list)]
	#indices <- c(as.numeric(args[1]))
	indices <- c(3)
	#indices <- seq_along(tcga_inside_folders)
	for (i in indices) {
	  cur <- tcga_inside_folders[[i]]
	  print(sprintf('working on %s (%d/%d)', cur, i, length(tcga_inside_folders)))
	  work_on_tcga_folder(cur, tcga_folder)
	}
	print("DONE")
}

run_organize_tcga()
