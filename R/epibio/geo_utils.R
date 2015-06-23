
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
  print(sprintf('read_l1_signal_file called on %s for %d rows', filename, nrows))
  # skip comments starts with: [#"]
  con <- gzfile(filename)
  lines <- readLines(con, n=40)
  close(con)
  comment_lines <- grepl('^[#"].*', lines)
  empty_lines <- grepl('^\\s*$', lines)
  skip <- which.min(comment_lines | empty_lines) - 1
  lines <- lines[skip+1:length(lines)]
  lines <- lines[!is.na(lines)]
  
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
    stop("read_l1_signal_file - Can't figure out correct sep")
  }
  dec <- NULL
  dec_list <- c(',')
  lines_without_header <- lines[2:length(lines)]
  for(dec_it in dec_list) {
    dec_count <- str_count(lines_without_header, dec_it)
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
  barcode_match_with_leading_X <- paste0("X", barcode_match)
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
    barcode_match, barcode_match_with_leading_X
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
