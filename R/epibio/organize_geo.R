
library(doParallel)
library(stringr)

source("config.R")
source("common.R")
source("geo_utils.R")
source("RnBeadsCommon.R")
args <- commandArgs(trailingOnly = TRUE)

list_series_id_files <- function(series_id_folder) {
  non_relevant_patterns <- c(
    "_[Pp]rocessed[._]", "_Summary_icc_M[.]",
    "upload_Beta[.]","_SampleMethylationProfile[.]",
    "_average_beta[.]", "_betas?[.]",
    "_geo_all_cohorts[.]", "_Results[.]",
    "_dasen[.]", "_NewGSMs[.]",
    "_Normalized_data[.]",
    "_Metrics[.]", "_qc[.]", "_BM_Oligos_samplsheet[.]")
  
  series_id_files <- list.files(series_id_folder, pattern="*.(txt.gz|csv.gz|tsv.gz)$")
  # filter non relevant files
  series_id_files <- series_id_files[!grepl(paste(non_relevant_patterns, 
                                                  collapse="|"), series_id_files)]
  series_id_files
}

# GSE62727 for example
readGeoL1DataWithIdats <- function(series_id_folder, series_id_orig, series_id_files, 
                                   output_filename, targets, all.series.info) {
  # Download idats
  splited_supplementary_file <- strsplit(targets$supplementary_file, ";")
  splited_supplementary_file <- lapply(splited_supplementary_file, trim)
  
  # Each item should be of length two (Green & Red)
  splited_supplementary_file_len_unique <- unique(unlist(lapply(splited_supplementary_file, 
                                                                function(x) length(x))))
  stopifnot(splited_supplementary_file_len_unique[[1]] == 2)
  stopifnot(length(splited_supplementary_file_len_unique) == 1)
  
  idat_folder <- file.path(series_id_folder, "idats")
  dir.create(idat_folder, recursive=TRUE, showWarnings=FALSE)

  targets$idat1_url <- sapply(splited_supplementary_file, "[", 1)
  targets$idat2_url <- sapply(splited_supplementary_file, "[", 2)
  
  targets$idat1_filename <- sapply(strsplit(targets$idat1_url, split='/', fixed=TRUE), tail, 1)
  targets$idat2_filename <- sapply(strsplit(targets$idat2_url, split='/', fixed=TRUE), tail, 1)

  targets$barcode <- gsub("_(Grn|Red).idat.gz", "", targets$idat1_filename)
  rownames(targets) <- targets$barcode
  
  for(i in 1:nrow(targets)) {
    target <- targets[i,]
    destfile <- file.path(idat_folder, target$idat1_filename)
    if(!file.exists(destfile)) {
      print(sprintf('downloading %s', target$idat1_url))
      download.file(target$idat1_url, destfile, "internal")
    }
    destfile <- file.path(idat_folder, target$idat2_filename)
    if(!file.exists(destfile)) {
      print(sprintf('downloading %s', target$idat2_url))
      download.file(target$idat2_url, destfile, "internal")
    }
  }

  # Work on idats
  print("working on idats")
  betas.table <- workOnIdatsFolder(idat_folder, targets, 'barcode')
  write_beta_values_table(output_filename, betas.table)
}

readGeoL1DataWithoutIdats <- function(series_id_folder, series_id_orig, series_id_files, 
                                      output_filename, targets, all.series.info) {
  nrows = 10000 # XXX (should be -1 on production)
  #nrows = -1
  
  this_targets = subset(targets, targets$series_id == series_id_orig)  
  filename_first_level <- levels(factor(this_targets$Filename))[[1]]
  this_all.series.info <- subset(all.series.info, 
                                 all.series.info$Filename == filename_first_level)
  print(series_id_files)
  series_id_fp <- file.path(series_id_folder, series_id_files)
  file_sizes <- sum(file.info(series_id_fp)$size/2**20)
  mem_limits <- tryCatch(memory.limit()/50, warning=function(x) NA)
  if(nrows == -1 & !is.na(mem_limits) && file_sizes > mem_limits)  {
    print(sprintf('GEO file too big for memory. skipping', basename(output_filename)))
  } else {
    p.values <- NULL
    problematic_unmeth_suffixes <- c("[. _-]?[Uu]nmethylated[. _-]?[Ss]ignal$", 
                                     "[_ .]{1,2}[Uu]nmethylated$",
                                     "_Unmethylated[.]Detection$",
                                     "[.]UM$",
                                     "[.]unmeth$")
    unmeth_suffixes <- c(problematic_unmeth_suffixes, "[._: ]Signal[_]?A$")
    meth_suffixes <- c("[. _-]?[Mm]ethylated[. _-]?[Ss]ignal$", 
                       "[_ .]{1,2}[Mm]ethylated$",
                       "_Methylated[.]Detection$",
                       "Signal[_]?B$", 
                       "_ M$", "[.]M$",
                       "[.]meth$",
                       # GSE58218 is strange
                       "[^h]ylated Signal")
    pvalue_suffixes <- c("_[ ]?pValue$",
                         "[. _:-]?Detection[. _-]?P[Vv]al(.\\d+)?$", 
                         "[.]Pval$", "[.]Detection$",
                         "[_ ]detection[_ ]p[-]?value[s]?$",
                         "[.]pval")
    suffixes = c(unmeth_suffixes, meth_suffixes, pvalue_suffixes)
    unmeth_files <- grep("Signal_A.NA|_unmeth|_Non_Methyl_", series_id_files)
    if(length(series_id_files) > 1 && length(series_id_files) - length(unmeth_files) == 1 ) {
      # works for GSE62992, GSE50498
      # => two files of raw signals: signal A and signal B, no pvals
      unmeth_signals <- read_l1_signal_file(series_id_fp[unmeth_files], nrows)
      meth_signals <- read_l1_signal_file(series_id_fp[-unmeth_files], nrows)
      
      colnum <- length(colnames(unmeth_signals))
      # remove unrelevant stuff from colnames
      samples.all <- gsub("[.]Signal_A","", colnames(unmeth_signals))
      relevant_samples <- get_relevant_samples(this_targets, samples.all, this_all.series.info)
      relevant_samples.all <- samples.all[relevant_samples, drop = FALSE]
      # assign unmethylated and methylated
      U <- data.matrix(unmeth_signals)[,relevant_samples, drop = FALSE]
      colnames(U) <- relevant_samples.all
      M <- data.matrix(meth_signals)[,relevant_samples, drop = FALSE]
      colnames(M) <- relevant_samples.all
    } else {
      # raw files with 3 columns for each sample
      # GSE36278 has two raw files for different samples
      if(length(series_id_fp) < 3) {
        signals <- do.call("cbind", lapply(series_id_fp, FUN=read_l1_signal_file, nrows))
      } else {
        stop('too many gz files')
      }
      
      if(grepl("ID_REF$|TargetID$", colnames(signals)[[1]], ignore.case = TRUE)) {
        # GSE46306, GSE48684
        if (colnames(signals)[[2]] == "ProbeID_A" && colnames(signals)[[3]] == "ProbeID_B") {
          # GSE50874
          rownames(signals) <- signals[, 1]
          signals <- signals[,-c(1,2,3)]
        } else {
          # GSE53162 which has two rownames columns
          rownames(signals) <- signals[, 1]
          signals <- signals[,-c(1)]
        }
      }
      
      # Remove AVG_Beta or Intensity columns (as in GSE52576, GSE50874, GSE53816)
      signals <- signals[!grepl("[._]AVG_Beta|[.]Intensity", colnames(signals))]
      
      # if there are two reps we use only the first (as in GSE53816)
      old_ncol <- ncol(signals)
      signals <- signals[!grepl("-rep2", colnames(signals))]
      if(ncol(signals) != old_ncol) {
        # remove the rep1 prefix from colnames
        colnames(signals) <- gsub('-rep1', '', colnames(signals))
      }
      
      # locate relevant samples
      if(length(colnames(signals)) == 0) {
        stop('signals is empty')
      }
      orig <- colnames(signals)
      unmeth_ids = grepl(paste(unmeth_suffixes, collapse="|"), orig)
      stopifnot(sum(unmeth_ids) > 0)
      # remove all unmeth expressions (because meth expressions are included in unmeth sometimes)
      orig <- gsub(paste(problematic_unmeth_suffixes, collapse="|"), "", orig)
      meth_ids =  grepl(paste(meth_suffixes, collapse="|"), orig)
      stopifnot(sum(meth_ids) > 0)
      # TODO - check GSE47627, GSE42118
      if(sum(unmeth_ids) != sum(meth_ids)) {
        print(sprintf("%d %d", sum(unmeth_ids), sum(meth_ids)))
        stop("different unmeth_ids and meth_ids!")
      }
      stopifnot(any(meth_ids & unmeth_ids) == FALSE)
      pval_ids = grepl(paste(pvalue_suffixes, collapse="|"), orig)
      
      # remove suffixes from colnames
      colnames(signals) <- mgsub(suffixes, character(length(suffixes)), colnames(signals))
      samples.all <- colnames(signals)[unmeth_ids]
      
      relevant_samples <- get_relevant_samples(this_targets, samples.all, this_all.series.info)
      # assign  unmethylated, methylated and pvalue matrices
      U <- data.matrix(signals[,unmeth_ids, drop = FALSE])[,relevant_samples, drop = FALSE]
      M <- data.matrix(signals[,meth_ids, drop = FALSE])[,relevant_samples, drop = FALSE]
      if(!is.null(pval_ids) & sum(pval_ids) > 0) {
        if(sum(unmeth_ids) != sum(pval_ids)) {
          print(sprintf("%d %d", sum(unmeth_ids), sum(pval_ids)))
          stop("different unmeth_ids and pval_ids!")
        }
        p.values <- data.matrix(signals[,pval_ids, drop = FALSE])[,relevant_samples, drop = FALSE]
      }
    }
    stopifnot(dim(this_targets)[[1]] == dim(U)[[2]])
    betas.table <- rnbReadL1Betas(this_targets, U, M, p.values)
    write_beta_values_table(output_filename, betas.table)
  }
}

#' Read GEO L1 data of given series id
#' 
#' @param series_id_orig
#' @param targets
#' @param all.series.info
#' @param study
#' @param type
#' @param geo_data_folder
#' @param generated_GEO_folder
#' 
#' @return nothing
readGeoL1Data <- function(series_id_orig, targets, all.series.info, study, type, 
                          geo_data_folder, generated_GEO_folder) {
  cat('\tReading ', series_id_orig, ": ")
  # handle samples which comes from multiple serieses
  series_id_vec <- unlist(strsplit(series_id_orig, ","))
  series_id <- NULL
  
  # check for idat files in the first series id
  idat_targets <- subset(targets, !is.na(supplementary_file))
  stopifnot(length(idat_targets) == length(targets))
  if(nrow(idat_targets) > 0) {
    series_id <- series_id_vec[[1]]
    series_id_folder <- file.path(geo_data_folder, series_id)
    series_id_files <- list_series_id_files(series_id_folder)
  }
  if(is.null(series_id)) {
    for(series_id_tmp in series_id_vec) {
      # check for data files
      series_id_folder <- file.path(geo_data_folder, series_id_tmp)
      series_id_files <- list_series_id_files(series_id_folder)
      if(length(series_id_files) > 0) {
        series_id <- series_id_tmp
        break
      }
    }
    if(is.null(series_id)) {
      stop(paste('no data files found for', series_id_orig))
    }
  }
  output_filename <- get_output_filename(generated_GEO_folder, series_id, study, type)
  if(file.exists(output_filename)) {
    print(sprintf('%s already exists. skipping', basename(output_filename)))
  } else {
    ptime1 <- proc.time()
    if(nrow(idat_targets) > 0 & length(series_id_files) == 0) {
      readGeoL1DataWithIdats(series_id_folder, series_id_orig, series_id_files, 
                             output_filename, targets, all.series.info)
    } else {
      readGeoL1DataWithoutIdats(series_id_folder, series_id_orig, series_id_files, 
                                output_filename, targets, all.series.info)
    }
    stime <- (proc.time() - ptime1)[3]
    cat("   in", stime, "seconds\n")
  }
}

workOnGEOTargets <- function(targets, all.series.info, geo_data_folder) {
  study <- levels(factor(targets$disease))[[1]]
  type <- levels(factor(targets$tissue_or_cell_type))[[1]]
  series_id <- levels(factor(targets$series_id))
  name <- create_name(study, type)
  print(sprintf("Reading %d samples of %s from %d serieses (study=%s, type=%s)", 
                nrow(targets), name, length(series_id), study, type))
  ret <- lapply(series_id, FUN=readGeoL1Data, targets, all.series.info, 
                study, type, geo_data_folder, generated_GEO_folder)
}


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
  # GSE41114 - has problem with the header columns - there is another ID_REF in it
	bad_list <- c(no_l1_list, not_released_list,
				  "GSE30338", "GSE37754", "GSE40360", "GSE40279", "GSE41826", 
				  "GSE43976", "GSE49377", "GSE48461", "GSE42882", "GSE46573",
				  "GSE55598", "GSE55438", "GSE56044", "GSE61044", "GSE61380",
				  "GSE42752", "GSE48684", "GSE49542", "GSE42372", "GSE32079",
				  "GSE46168", "GSE47627", "GSE61151", "GSE32146", "GSE41114")
	wait_list <- c("GSE62924", "GSE51245", "GSE38266")
	ignore_list <- paste0(joined_folder, "/", c(bad_list, wait_list), ".txt")
	#external_disk_data_path <- '/cs/icore/joshua.moss/dor/atlas'
	geo_data_folder <- file.path(external_disk_data_path, 'GEO')
	stopifnot(file.exists(geo_data_folder))
	only_vec <- list.files(geo_data_folder)
	#only_vec <- c("GSE46306") # XXX # TODO - see GSE46306 is working
	# for hai
	#working_vec <- c('GSE36278', 'GSE52556', 'GSE52576', 'GSE61160', 'GSE53924', 'GSE42752', 'GSE30338', 'GSE32283', 'GSE41826', 'GSE42882', 'GSE46573')
	only_vec <- c('GSE54776', 'GSE55712', 'GSE61107', 'GSE61380')
	#only_vec <- c('GSE48472', 'GSE43414', 'GSE50798', 'GSE48461', 'GSE44661', 'GSE53816', 'GSE49576')
	#only_vec <- c('GSE58218', 'GSE59524', 'GSE49377', 'GSE61431', 'GSE62727', 'GSE31848')
	
	only_list <- paste0(joined_folder, "/", c(only_vec), ".txt")
	joined_files <- joined_files[(joined_files %in% only_list) & !(joined_files %in% ignore_list)]
	stopifnot(length(joined_files) > 0)
	print("joined_files:")
	print(joined_files)

	all.series.info <- do.call("rbind", lapply(joined_files, FUN=read_joined_file))
	# get only relevant samples
	relevant.samples.idx <- which(as.numeric(all.series.info$relevant) == 1)
	pheno <- all.series.info[relevant.samples.idx, ]
  
	# Fix pheno values
	col_vec <- c('series_id', 'title', 'cell_type', 'tissue', 'disease')
	missing_cond <- (is.na(pheno$tissue) | pheno$tissue=='') & (is.na(pheno$cell_type) | pheno$cell_type=='')
	duplicate_cond <- !is.na(pheno$tissue) & pheno$tissue!='' & !is.na(pheno$cell_type) & pheno$cell_type!=''
	#missing_both <- subset(pheno, missing_cond)
	#duplicate_names <- subset(pheno, duplicate_cond)
	#write.csv(missing_both[,col_vec], file='missing_both.csv')
	#write.csv(duplicate_names[,col_vec], file='duplicate_names.csv')
	
	pheno <- subset(pheno, !(missing_cond | duplicate_cond))
	pheno$tissue_or_cell_type <- paste3(pheno$tissue, pheno$cell_type)
  
	splited_targets <- split(pheno, list(pheno$disease, pheno$tissue_or_cell_type), drop=TRUE)
	geo_data_folder <- file.path(external_disk_data_path, 'GEO')
	indices <- get_indices_to_runon(splited_targets, args)
	#indices <- c(3)
	logger.start(fname=NA)
	num.cores <- detectCores()/2
	#parallel.setup(num.cores)
	for (i in indices) {
	  print(sprintf('working on %d/%d', i, length(splited_targets)))
	  workOnGEOTargets(splited_targets[[i]], all.series.info, geo_data_folder)
	}
	#parallel.disable()
  
	print("DONE")  
}

run_organize_geo()
