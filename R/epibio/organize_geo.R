
library(doParallel)

source("config.R")
source("common.R")
source("geo_utils.R")
args <- commandArgs(trailingOnly = TRUE)

run_organize_geo <- function() {
	dir.create(generated_GEO_folder, recursive=TRUE, showWarnings=FALSE)

	joined_folder <- file.path(data_folder, "global/GEO/joined")
	joined_files <- list.files(joined_folder, full.names = TRUE, pattern="*.txt")

	# == skip serieses ==
	# GEOs which I don't know how to parse:
	# - no l1 signals txt file
	# - different parsing on l1 txt file
	no_l1_list <- c("GSE37965", "GSE39279", "GSE39560", "GSE41169", "GSE53924", 
					"GSE39141", "GSE34777")
	not_released_list <- c("GSE62003", "GSE49064")
	bad_list <- c(no_l1_list, not_released_list,
				  "GSE30338", "GSE37754", "GSE40360", "GSE40279", "GSE41826", 
				  "GSE43976", "GSE49377", "GSE48461", "GSE42882", "GSE46573",
				  "GSE55598", "GSE55438", "GSE56044", "GSE61044", "GSE61380",
				  "GSE42752", "GSE48684", "GSE49542", "GSE42372", "GSE32079",
				  "GSE46168", "GSE47627", "GSE61151", "GSE32146")
	wait_list <- c("GSE62924", "GSE51245", "GSE38266", "GSE46306", "GSE59685")
	ignore_list <- paste0(joined_folder, "/", c(bad_list, wait_list), ".txt")
	#external_disk_data_path <- '/cs/icore/joshua.moss/dor/atlas'
	geo_data_folder <- file.path(external_disk_data_path, 'GEO')
	stopifnot(file.exists(geo_data_folder))
	only_vec <- list.files(geo_data_folder)
	#only_vec <- c("GSE42118") # XXX
	only_list <- paste0(joined_folder, "/", c(only_vec), ".txt")
	joined_files <- joined_files[(joined_files %in% only_list) & !(joined_files %in% ignore_list)]
	stopifnot(length(joined_files) > 0)
	print("joined_files:")
	print(joined_files)

	all.series.info <- do.call("rbind", lapply(joined_files, FUN=read_joined_file))
	# Remove serieses with idats
	all.series.info <- subset(all.series.info, is.na(supplementary_file))
	# get only relevant samples
	relevant.samples.idx <- which(as.numeric(all.series.info$relevant) == 1)
	pheno <- all.series.info[relevant.samples.idx, ]

	# Fix pheno values
	col_vec <- c('series_id', 'title', 'cell_type', 'tissue', 'disease')
	missing_cond <- (is.na(pheno$tissue) | pheno$tissue=='') & (is.na(pheno$cell_type) | pheno$cell_type=='')
	missing_both <- subset(pheno, missing_cond)
	duplicate_cond <- !is.na(pheno$tissue) & pheno$tissue!='' & !is.na(pheno$cell_type) & pheno$cell_type!=''
	duplicate_names <- subset(pheno, duplicate_cond)
	#write.csv(missing_both[,col_vec], file='missing_both.csv')
	#write.csv(duplicate_names[,col_vec], file='duplicate_names.csv')

	pheno <- subset(pheno, !(missing_cond | duplicate_cond))
	pheno$tissue_or_cell_type <- paste3(pheno$tissue, pheno$cell_type)

	splited_targets <- split(pheno, list(pheno$disease, pheno$tissue_or_cell_type), drop=TRUE)
	geo_data_folder <- file.path(external_disk_data_path, 'GEO')
	indices <- get_indices_to_runon(splited_targets, args)
	indices <- c(12)
	logger.start(fname=NA)
	num.cores <- detectCores()/2
	#parallel.setup(num.cores)
	for (i in indices) {
	  print(sprintf('working on %d/%d', i, length(splited_targets)))
	  workOnTargets(splited_targets[[i]], all.series.info, geo_data_folder)
	}
	#parallel.disable()
	print("DONE")  
}

run_organize_geo()
