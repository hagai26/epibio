
library(RnBeads)
library(doParallel)

source("config.R")
source("common.R")
source("geo_utils.R")



#' Build new RnbeadsRawset
rnbReadL1Betas <- function(targets, U, M, p.values) {
  pheno <- targets[, c('description','tissue','cell_type','disease')]
  rnb.set <- RnBeadRawSet(pheno, U=U, M=M, p.values=p.values, useff=FALSE)
  betas.table <- process_rnb_set_to_betas(rnb.set, !is.null(p.values))
  betas.table
}

#' get the relevant samples in samples.all
#' 
#' @param this_targets
#' @param samples.all
#' @param this_all.series.info
#' 
#' @return relevant samples boolean vector
get_relevant_samples <- function(this_targets, samples.all, this_all.series.info) {
  
  # GSE43414 has characteristics_ch1 like:
  # GSM1068923 subjectid: NA;\tbarcode: 6057825014_R06C02.1;\tlunnonetal: FALSE;\ttissue_code: NA;\tbraak.stage: NA;\tSex: NA;\tad.disease.status: NA;\tage.brain: NA;\tage.blood: NA;\tsource tissue: cerebellum
  barcode_match <- str_match(this_targets$characteristics_ch1, "barcode: ([^; ]+)")[,c(2)]      
  title_last_word <- gsub(".* ", '', this_targets$title)
  # sometimes its numbers and they add S to each number (as in GSE53816)
  title_last_word2 <- paste0('S', title_last_word)
  try_match_list <- list(
    this_targets$description, 
    this_targets$source_name_ch1,
    # first word
    gsub("[ ;].*", '', this_targets$source_name_ch1), 
    gsub("[ ;].*", '', this_targets$description), 
    # last word
    gsub(".* ", '', this_targets$source_name_ch1), 
    title_last_word, 
    title_last_word2,
    gsub(".*[ ;\t]", '', this_targets$description),
    # middle one barcode
    barcode_match
  )
  
  match_all <- lapply(try_match_list, function(x) match(samples.all, as.character(x)))
  # find most match locations
  match_all_no_na <- lapply(match_all, function(x) x[!is.na(x)])
  most_match_index <- which.max(sapply(match_all_no_na, FUN=length))
  relevant.samples.loc <- match_all[[most_match_index]]
  
  if(all(is.na(relevant.samples.loc))) {
    if(length(samples.all) == dim(this_all.series.info)[[1]]) {
      relevant_samples <- (this_all.series.info$description %in% this_targets$description) & 
                    (this_all.series.info$source_name_ch1 %in% this_targets$source_name_ch1) &
                    (this_all.series.info$title %in% this_targets$title) &
                    (this_all.series.info$extract_protocol_ch1 %in% this_targets$extract_protocol_ch1)
    } else {
      stop('get_relevant_samples failed')
    }
  } else {
    relevant_samples <- !is.na(relevant.samples.loc)
    # each sample should be only once
    # Remove NA's
    relevant.samples.loc_no_na <- relevant.samples.loc[relevant_samples]
    stopifnot(length(unique(relevant.samples.loc_no_na)) == length(relevant.samples.loc_no_na))
  }
  relevant_samples
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
  for(series_id_tmp in series_id_vec) {
    series_id_folder <- file.path(geo_data_folder, series_id_tmp)
    series_id_files <- list.files(series_id_folder, pattern="*.(txt.gz|csv.gz|tsv.gz)$")
    # filter non relevant files
    non_relevant_patterns <- c(
                        "_[Pp]rocessed[._]", 
                        "upload_Beta[.]",
                        "_SampleMethylationProfile[.]",
                        "_average_beta[.]", "_betas?[.]",
                        "_geo_all_cohorts[.]",
                        "_dasen[.]", "_NewGSMs[.]")
    series_id_files <- series_id_files[!grepl(paste(non_relevant_patterns, collapse="|"), 
                                              series_id_files)]
    if(length(series_id_files) > 0) {
      series_id <- series_id_tmp
      break
    }
  }
  if(is.null(series_id)) {
    stop(paste('no data files found for', series_id_orig))
  }
  output_filename <- get_output_filename(generated_GEO_folder, series_id, study, type)
  if(file.exists(output_filename)) {
    print(sprintf('%s already exists. skipping', basename(output_filename)))
  } else {
    ptime1 <- proc.time()
    print(series_id_files)
    this_targets = subset(targets, targets$series_id == series_id_orig)  
    filename_first_level <- levels(factor(this_targets$Filename))[[1]]
    this_all.series.info <- subset(all.series.info, 
                                   all.series.info$Filename == filename_first_level)
    series_id_fp <- file.path(series_id_folder, series_id_files)
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
                      "[._: ]Signal[_]?B$", 
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
    
    nrows = 5000 # XXX (should be -1 on production)
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
      pval_ids = grepl(paste(pvalue_suffixes, collapse="|"), orig)
      # TODO - check GSE47627, GSE42118
      if(sum(unmeth_ids) != sum(meth_ids)) {
        print(sprintf("%d %d", sum(unmeth_ids), sum(meth_ids)))
        stop("different unmeth_ids and meth_ids!")
      }
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
        signals_pval <- data.matrix(signals[,pval_ids, drop = FALSE])[,relevant_samples, drop = FALSE]
        if(length(signals_pval) > 0) {
          # convert the decimal comma into a dot (as in GSE29290)
          x <- lapply(as.data.frame(signals_pval), function(x) gsub(",", ".", x, fixed = TRUE))
          signals_pval <- data.frame(x, row.names=rownames(signals_pval), stringsAsFactors=FALSE)
          p.values <- data.matrix(signals_pval)
        }
      }
    }
    stopifnot(dim(this_targets)[[1]] == dim(U)[[2]])
    if (nrow(this_targets) > 1) {
      betas.table <- rnbReadL1Betas(this_targets, U, M, p.values)
      write_beta_values_table(output_filename, betas.table)
    } else {
      # Error in checkSlotAssignment(object, name, value) : 
      # assignment of an object of class “numeric” is not valid for slot ‘meth.sites’ in an object of class “RnBeadSet”; is(value, "matrixOrff") is not TRUE
      print("Got only one target - rnbeads raises error on these cases - should fix it - TODO")
    }
    
    stime <- (proc.time() - ptime1)[3]
    cat("   in", stime, "seconds\n")
  }
}

workOnTargets <- function(targets, all.series.info, geo_data_folder) {
  study <- levels(factor(targets$disease))[[1]]
  type <- levels(factor(targets$tissue_or_cell_type))[[1]]
  series_id <- levels(factor(targets$series_id))
  name <- create_name(study, type)
  cat("Reading", nrow(targets), "samples of", name, "from", length(series_id), "serieses\n")
  ret <- lapply(series_id, FUN=readGeoL1Data, targets, all.series.info, 
                study, type, geo_data_folder, generated_GEO_folder)
}

dir.create(generated_GEO_folder, recursive=TRUE, showWarnings=FALSE)

joined_folder <- file.path(data_folder, "global/GEO/joined")
joined_files <- list.files(joined_folder, full.names = TRUE, pattern="*.txt")

# == skip serieses ==
# GEOs which I don't know how to parse:
# - no l1 signals txt file
# - different parsing on l1 txt file
no_l1_list <- c("GSE37965", "GSE39279", "GSE39560", "GSE41169", "GSE53924", "GSE39141")
not_released_list <- c("GSE62003")
bad_list <- c(no_l1_list, not_released_list,
              "GSE30338", "GSE37754", "GSE40360", "GSE40279", "GSE41826", 
              "GSE43976", "GSE49377", "GSE48461", "GSE42882", "GSE46573",
              "GSE55598", "GSE55438", "GSE56044", "GSE61044", "GSE61380",
              "GSE42752", "GSE48684", "GSE49542", "GSE42372", "GSE32079",
              "GSE46168")
wait_list <- c("GSE62924")
ignore_list <- paste0(joined_folder, "/", c(bad_list, wait_list), ".txt")

geo_data_folder <- file.path(external_disk_data_path, 'GEO')
stopifnot(file.exists(geo_data_folder))
only_vec <- list.files(geo_data_folder)
#only_vec <- c("GSE62924") # XXX
only_list <- paste0(joined_folder, "/", c(only_vec), ".txt")
joined_files <- joined_files[(joined_files %in% only_list)]
joined_files <- joined_files[!(joined_files %in% ignore_list)]
#joined_files <- head(joined_files, 26) # XXX
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
write.csv(missing_both[,col_vec], file='missing_both.csv')

duplicate_cond <- !is.na(pheno$tissue) & pheno$tissue!='' & !is.na(pheno$cell_type) & pheno$cell_type!=''
duplicate_names <- subset(pheno, duplicate_cond)
write.csv(duplicate_names[,col_vec], file='duplicate_names.csv')

pheno <- subset(pheno, !(missing_cond | duplicate_cond))
pheno$tissue_or_cell_type <- paste3(pheno$tissue, pheno$cell_type)

splited_targets <- split(pheno, list(pheno$disease, pheno$tissue_or_cell_type), drop=TRUE)
geo_data_folder <- file.path(external_disk_data_path, 'GEO')
logger.start(fname=NA)
num.cores <- 2
parallel.setup(num.cores)
for (i in seq_along(splited_targets)) {
  print(sprintf('working on %d/%d', i, length(splited_targets)))
  workOnTargets(splited_targets[[i]], all.series.info, geo_data_folder)
}
parallel.disable()
print("DONE")
