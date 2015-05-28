
library(RnBeads)
library(stringr)  # str_count

source("common.R")


read_joined_file <- function(filename) {
  t <- read.table(filename, sep='\t', header=TRUE, row.names=1, fill=TRUE, 
                  na.strings=c("NA", "0"), quote="\"", stringsAsFactors=FALSE,
                  strip.white=TRUE)
  if(length(rownames(t)) > 0) {
    df <- data.frame(Filename=filename, t)
  } else {
    # on case that t is empty
    df <- data.frame(t)
  }
  df
}

read_l1_signal_file <- function(filename, nrows) {
  # skip comments starts with: [#"]
  con <- gzfile(filename)
  lines <- readLines(con, n=30)
  close(con)
  comment_lines <- grepl('^[#"].*', lines)
  empty_lines <- grepl('^\\s*$', lines)
  skip <- which.min(comment_lines | empty_lines) - 1
  lines <- lines[skip+1:length(lines)]
  
  # choose correct sep & dec
  sep <- NULL
  MIN_COLS = 4
  sep_list <- c('\t', ',', ' ')
  for(sep_it in sep_list) {
    sep_count <- str_count(lines, sep_it)
    if(mean(sep_count) >= MIN_COLS) {
      sep <- sep_it
      break
    }
  }
  if(is.null(sep)) {
    stop("Can't figure out correct sep")
  }
  dec <- NULL
  dec_list <- c(',')
  for(dec_it in dec_list) {
    dec_count <- str_count(lines[2:length(lines)], dec_it)
    if(mean(dec_count) >= MIN_COLS) {
      dec <- dec_it
      break
    }
  }
  if(is.null(dec)) {
    dec = "."
  }
  
  # skip comments which doesn't look like comments by checking for sep count inside them
  cols_count <- sapply(strsplit(lines, sep), function(x) sum(x!=""))
  more_skip <- which.max(cols_count > MIN_COLS) - 1
  skip <- skip + more_skip
  lines <- lines[skip+1:length(lines)]
  
  # handle GSE50759 - which has first row with only: "unmeth" "meth"   "pval" column names
  first_row_names <- Filter(function(x) x != "", strsplit(lines[[1]], sep)[[1]])
  first_row_names_unique <- unique(first_row_names)
  colnames_suffixes = NULL
  if(length(first_row_names_unique) < MIN_COLS && "meth" %in% first_row_names_unique) {
    skip <- skip + 1
    colnames_suffixes <- paste0(".", first_row_names)
  }
  
  id_ref_on_other_line <- grepl(paste0("^ID_REF", sep), lines[[2 + skip]])
  if(id_ref_on_other_line) {
    skip <- skip + 1
  }

  # turn off the interpretation of comments
  # because there are samples names with # sometimes (as in GSE58280)
  con <- gzfile(filename)
  t <- read.table(con, header=TRUE, row.names=1, skip=skip, 
                  sep=sep, dec=dec, nrows=nrows, check.names=FALSE, 
                  stringsAsFactors=FALSE, comment.char="")
  # remove columns which doesn't have labels on header (like in GSE32146)
  good_cols <- colnames(t)[colnames(t) != ""]
  out_table <- t[good_cols]
  
  if(!is.null(colnames_suffixes)) {
    colnames(out_table) <- paste0(colnames(out_table), colnames_suffixes)
  }
  
  # convert the decimal comma into a dot (even after using dec=, - should fix the case of .0 -> 0
  if(dec == ',') {
    x <- lapply(out_table, function(x) gsub(",", ".", x, fixed = TRUE))
    out_table <- data.frame(x, row.names=rownames(out_table), stringsAsFactors=FALSE)
  }
  
  out_table
}


#' Build new RnbeadsRawset
rnbReadL1Betas <- function(targets, U, M, p.values) {
  pheno <- targets[, c('description','tissue','cell_type','disease')]
  #next three lines added to avoid error in cases of no pvalues ##josh##
  if (is.null(p.values)){
      p.values <- matrix(0,nrow=nrow(M),ncol=ncol(M),dimnames=list(rownames(M),colnames(M)))
  }
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
  title_last_word_no_leading_zeros <- gsub("(?<![0-9])0+", "", title_last_word, perl = TRUE)
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
    title_last_word_no_leading_zeros,
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
  non_relevant_patterns <- c(
    "_[Pp]rocessed[._]", "_Summary_icc_M[.]",
    "upload_Beta[.]","_SampleMethylationProfile[.]",
    "_average_beta[.]", "_betas?[.]",
    "_geo_all_cohorts[.]", "_Results[.]",
    "_dasen[.]", "_NewGSMs[.]")
  
  for(series_id_tmp in series_id_vec) {
    series_id_folder <- file.path(geo_data_folder, series_id_tmp)
    series_id_files <- list.files(series_id_folder, pattern="*.(txt.gz|csv.gz|tsv.gz)$")
    # filter non relevant files
    series_id_files <- series_id_files[!grepl(paste(non_relevant_patterns, 
                                                    collapse="|"), series_id_files)]
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
    nrows = 50000 # XXX (should be -1 on production)
    #nrows = -1
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
