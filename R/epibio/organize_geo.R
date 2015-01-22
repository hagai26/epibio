
library(RnBeads)

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
    non_relevant_patterns <- c("_processed_")
    series_id_files <- series_id_files[!grepl(paste(non_relevant_patterns, collapse="|"), series_id_files)]
    if(length(series_id_files) > 0) {
      series_id <- series_id_tmp
      break
    }
  }
  if(is.null(series_id)) {
    stop(paste('no data files found for', series_id_orig))
  }
  
  ptime1 <- proc.time()
  print(series_id_files)
  this_targets = subset(targets, targets$series_id == series_id_orig)  
  filename_first_level <- levels(factor(this_targets$Filename))[[1]]
  this_all.series.info <- subset(all.series.info, all.series.info$Filename == filename_first_level)
  series_id_fp <- file.path(series_id_folder, series_id_files)
  p.values <- NULL
  
  unmeth_suffixes = c("[. _-]?[Uu]nmethylated[. _-]?[Ss]ignal$", 
                      "[_ .]{1,2}[Uu]nmethylated$",
                      "_Unmethylated[.]Detection$",
                      "[._: ]Signal[_]?A$", 
                      ".UM$")
  meth_suffixes = c("[. _-]?[Mm]ethylated[. _-]?[Ss]ignal$", 
                    "[_ .]{1,2}[Mm]ethylated$",
                    "_Methylated[.]Detection$",
                    "[._: ]Signal[_]?B$", 
                    "_ M$", ".M$")
  pvalue_suffixes = c("_[ ]?pValue$",
                      "[. _:-]?Detection[. _-]?P[Vv]al(.\\d+)?$", 
                      "[.]Pval$", "[.]Detection$",
                      "[_ ]detection[_ ]p[-]?value[s]?$")
  suffixes = c(unmeth_suffixes, meth_suffixes, pvalue_suffixes)
  
  unmeth_files <- grep("Signal_A.NA|_unmeth", series_id_files)
  nrows = 4000 # XXX (should be -1 on production)
  if(length(series_id_files) > 1 && length(series_id_files) - length(unmeth_files) == 1 ) {
    # works for GSE62992
    # => two files of raw signals: signal A and signal B, no pvals
    unmeth_signals <- read_l1_signal_file(series_id_fp[unmeth_files], nrows)
    meth_signals <- read_l1_signal_file(series_id_fp[-unmeth_files], nrows)
    
    colnum <- length(colnames(unmeth_signals))
    samples.all <- gsub("[.]Signal_A","", colnames(unmeth_signals))
    
    if(length(samples.all) == dim(this_all.series.info)[[1]]) {
      v <- (this_targets$description %in% this_all.series.info$description) & (this_targets$source_name_ch1 %in% this_all.series.info$source_name_ch1)
      relevant.samples.loc <- which(v)
    } else {
      stop('try other option 3')
    }
    
    # assign  unmethylated, methylated and pvalue matrices
    U <- data.matrix(unmeth_signals)
    colnames(U) <- samples.all
    M <- data.matrix(meth_signals)
    colnames(M) <- samples.all
  } else {
    # works for GSE32079, GSE29290, GSE57767, GSE61653, ..
    # => one raw file with 3 columns for each sample
    
    # GSE36278 has two raw files
    if(length(series_id_fp) == 2) {
      signals <- do.call("cbind", lapply(series_id_fp, FUN=read_l1_signal_file, nrows))
    } else {
      signals <- read_l1_signal_file(series_id_fp, nrows)
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
    
    # Remove AVG_Beta or Intensity columns (as in GSE52576, GSE50874)
    signals <- signals[!grepl("[.]AVG_Beta|[.]Intensity", colnames(signals))]
    # locate relevant samples
    colnum <- length(colnames(signals))
    if(colnum == 0) {
      stop('signals is empty')
    }
    orig <- colnames(signals)
    
    unmeth_ids = grepl(paste(unmeth_suffixes, collapse="|"), orig)
    meth_ids =  grepl(paste(meth_suffixes, collapse="|"), orig)
    pval_ids = grepl(paste(pvalue_suffixes, collapse="|"), orig)
    # TODO - check GSE47627, GSE42118

    # remove suffixes from colnames
    colnames(signals) <- mgsub(suffixes, character(length(suffixes)), colnames(signals))
    #print (colnames(signals))
    samples.all <- colnames(signals)[unmeth_ids]
    
    try_match_list <- list(this_targets$description, 
                       # first word
                       gsub("[ ;].*", '', this_targets$source_name_ch1), 
                       gsub("[ ;].*", '', this_targets$description), 
                       # last word
                       gsub(".* ", '', this_targets$source_name_ch1), 
                       gsub(".* ", '', this_targets$title), 
                       gsub(".*[ ;\t]", '', this_targets$description)
    )
    
    match_all <- lapply(try_match_list, function(x) match(samples.all, as.character(x)))
    # Remove NA's
    match_all_no_na <- lapply(match_all, function(x) x[!is.na(x)])
    # find most match locations
    most_match_index <- which.max(sapply(match_all_no_na, FUN=length))
    relevant.samples.loc <- match_all_no_na[[most_match_index]]

    if(all(is.na(relevant.samples.loc))) {
      if(length(samples.all) == dim(this_all.series.info)[[1]]) {
        v <- (this_targets$description %in% this_all.series.info$description) & (this_targets$source_name_ch1 %in% this_all.series.info$source_name_ch1)
        relevant.samples.loc <- which(v)
      } else {
        stop('try other option 1')
      }
    }
    
    # assign  unmethylated, methylated and pvalue matrices
    U <- data.matrix(signals[,unmeth_ids, drop = FALSE])[,relevant.samples.loc, drop = FALSE]
    M <- data.matrix(signals[,meth_ids, drop = FALSE])[,relevant.samples.loc, drop = FALSE]
    if(!is.null(pval_ids)) {
      signals_pval <- signals[,pval_ids, drop = FALSE]
      if(length(signals_pval) > 0) {
        # convert the decimal comma into a dot (as in GSE29290)
        signals_pval <- data.frame(lapply(signals_pval, function(x) gsub(",", ".", x, fixed = TRUE)), 
                                 row.names=rownames(signals_pval), stringsAsFactors=FALSE)
        p.values <- data.matrix(signals_pval)[,relevant.samples.loc, drop = FALSE]
      }
    }
  }
  if (nrow(this_targets) > 1) {
    if(dim(this_targets)[[1]] != dim(U)[[2]]) {
      print('different dim!')
    }
    
    betas.table <- rnbReadL1Betas(this_targets, U, M, p.values)
    write_beta_values_table(generated_GEO_folder, series_id, study, type, betas.table)
  } else {
    # Error in checkSlotAssignment(object, name, value) : 
    # assignment of an object of class “numeric” is not valid for slot ‘meth.sites’ in an object of class “RnBeadSet”; is(value, "matrixOrff") is not TRUE
    print("Got only one target - rnbeads raises error on these cases - should fix it - TODO")
  }
  
  stime <- (proc.time() - ptime1)[3]
  cat("   in", stime, "seconds\n")
}

workOnTargets <- function(targets, all.series.info, geo_data_folder) {
  study <- levels(factor(targets$disease))[[1]]
  type <- levels(factor(targets$tissue))[[1]]
  series_id <- levels(factor(targets$series_id))
  name <- create_name(study, type)
  cat("Reading", nrow(targets), "samples of", name, "from", length(series_id), "serieses\n")
  ret <- lapply(series_id, FUN=readGeoL1Data, targets, all.series.info, study, type, geo_data_folder, generated_GEO_folder)
}

dir.create(generated_GEO_folder, recursive=TRUE, showWarnings=FALSE)

joined_folder <- file.path(data_folder, "global/GEO/joined")
joined_files <- list.files(joined_folder, full.names = TRUE, pattern="*.txt")

# == skip serieses ==
# GEOs which I don't know how to parse:
# - no l1 signals txt file
# - different parsing on l1 txt file
no_l1_list <- c("GSE37965", "GSE39279", "GSE39560", "GSE41169", "GSE53924")
not_released_list <- c("GSE62003")
bad_list <- c(no_l1_list, not_released_list,
              "GSE30338", "GSE37754", "GSE40360", "GSE40279", "GSE41826", 
              "GSE43976", "GSE49377", "GSE48461", "GSE42882", "GSE46573",
              "GSE55598", "GSE55438", "GSE56044", "GSE61044", "GSE61380",
              "GSE42752")
# GEOs which I still don't have
too_big <- c("GSE53816", "GSE54882", "GSE58218", "GSE59685", "GSE61151")
wait_list <- c(too_big)
# working GEOs
working_list <- c("GSE32079", "GSE38266", "GSE35069", "GSE32283", "GSE36278", 
                  "GSE29290", "GSE32146", "GSE37362", "GSE40853", "GSE59157",
                  "GSE39958", "GSE41114", "GSE42372", "GSE41273", "GSE43091", 
                  "GSE44661", "GSE43293", "GSE43298", "GSE46394", "GSE44684",
                  "GSE42118", "GSE42119", "GSE43414", "GSE47512", "GSE62640",
                  "GSE45187", "GSE48325", "GSE49031", "GSE49656", "GSE49576", 
                  "GSE31803", "GSE49542", "GSE44667", "GSE45353", "GSE51758",
                  "GSE50498", "GSE50774", "GSE50759", "GSE52826", "GSE53128",
                  "GSE52576", "GSE50874", "GSE52731", "GSE52401", "GSE50798", 
                  "GSE49393", "GSE47627", "GSE53740", "GSE52113", "GSE53162",
                  "GSE46306", "GSE48684", "GSE54399", "GSE54503", "GSE54415",
                  "GSE54880", "GSE55571", "GSE54776", "GSE54670", "GSE57767",
                  "GSE55712", "GSE57831", "GSE55734", "GSE56420", "GSE53840",
                  "GSE58280", "GSE61653", "GSE63499", "GSE62992", "GSE61256",
                  "GSE60753", "GSE61256", "GSE61259", "GSE61257", "GSE61431",
                  "GSE61258", "GSE58651", "GSE38268")
ignore_list <- paste0(joined_folder, "/", c(bad_list, wait_list), ".txt")

geo_data_folder <- file.path(external_disk_data_path, 'GEO')
only_vec <- list.files(geo_data_folder)
#only_vec <- c("GSE42752") # XXX
only_list <- paste0(joined_folder, "/", c(only_vec), ".txt")
joined_files <- joined_files[(joined_files %in% only_list)]
joined_files <- joined_files[!(joined_files %in% ignore_list)]
#joined_files <- head(joined_files, 26) # XXX
print("joined_files:")
print(joined_files)

all.series.info <- do.call("rbind", lapply(joined_files, FUN=read_joined_file))
# Remove serieses with idats
all.series.info <- subset(all.series.info, is.na(supplementary_file))

# get only relevant samples
relevant.samples.idx <- which(as.numeric(all.series.info$relevant) == 1)
pheno <- all.series.info[relevant.samples.idx, ]
splited_targets <- split(pheno, list(pheno$disease, pheno$tissue), drop=TRUE)

geo_data_folder <- file.path(external_disk_data_path, 'GEO')
ret <- lapply(splited_targets, FUN=workOnTargets, all.series.info, geo_data_folder)

print("DONE")
