
library(stringr)  # str_count

source("common.R")


read_joined_file <- function(filename) {
  t <- read.table(filename, sep='\t', header=TRUE, row.names=1, fill=TRUE, 
                  na.strings=c("NA", "0"), quote="\"", stringsAsFactors=FALSE)
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
  
  # choose correct sep
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
  
  # skip comments which doesn't look like comments by checking for sep count inside them
  cols_count <- sapply(strsplit(lines, sep), function(x) sum(x!=""))
  more_skip <- which.max(cols_count > MIN_COLS) - 1
  skip <- skip + more_skip
  lines <- lines[skip+1:length(lines)]
  
  # handle GSE50759 - which has first row with only: "unmeth" "meth"   "pval" column names
  first_row_names <- Filter(function(x) x!="", unique(strsplit(lines[[1]], sep)[[1]]))
  if(length(first_row_names) < MIN_COLS) {
    skip <- skip + 1
  }
  
  id_ref_on_other_line <- grepl(paste0("^ID_REF", sep), lines[[2 + skip]])
  if(id_ref_on_other_line) {
    skip <- skip + 1
  }
  

  # turn off the interpretation of comments
  # because there are samples names with # sometimes (as in GSE58280)
  t <- read.table(gzfile(filename), header=TRUE, row.names=1, skip=skip, sep=sep, dec='.', nrows=nrows,
                  check.names=FALSE, stringsAsFactors=FALSE, comment.char="")
  # remove columns which doesn't have labels on header (like in GSE32146)
  good_cols <- colnames(t)[colnames(t) != ""]
  t[good_cols]
}
