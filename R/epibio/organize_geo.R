
library(RnBeads)
library(stringr)

source("config.R")
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

read_l1_signal_file <- function(filename) {
  nrows = 250

  # skip comments starts with: [#"]
  lines <- readLines(filename, n=30)
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
  
  t <- read.table(filename, header=TRUE, row.names=1, skip=skip, sep=sep, dec = ".",
                  nrows=nrows,
                  check.names=FALSE, stringsAsFactors=FALSE)
  # remove columns which doesn't have labels on header (like in GSE32146)
  good_cols <- colnames(t)[colnames(t) != ""]
  t[good_cols]
}

rnb_read_l1_betas <- function(targets, U, M, p.values) {
  pheno <- targets[, c('description','tissue','cell_type','disease')]
  rnb.set <- new('RnBeadRawSet', pheno, U=U, M=M, p.values=p.values, useff=FALSE)
  
  logger.start(fname=NA)
  rnb.options(disk.dump.big.matrices=TRUE)
  rnb.options(enforce.memory.management=TRUE)
  
  tryCatch({
    rnb.set <- rnb.execute.snp.removal(rnb.set)$dataset
  }, error = function(err) {
    # TODO
    
    # on GSE36278 this causes stops with error:
    # Error in checkSlotAssignment(object, name, value) : 
    # assignment of an object of class “numeric” is not valid for slot ‘M’ in an object of class “RnBeadRawSet”; is(value, "matrixOrffOrNULL") is not TRUE
    
    # on GSE42118 (and GSE52576, GSE44667):
    # <simpleError in checkSlotAssignment(object, name, value): 
    # assignment of an object of class “integer” is not valid for slot ‘M’ in an object of class “RnBeadRawSet”; is(value, "matrixOrffOrNULL") is not TRUE>
    
    # currently, we ignore this error and doesn't call snp.removal
    print(err)
  })
  
  betas.table <- meth(rnb.set, row.names=TRUE)
  if(!is.null(p.values)) {
    pvalue.high <- which(dpval(rnb.set) > 0.05, arr.ind=TRUE)
    betas.table[pvalue.high[,'row'], pvalue.high[,'col']] <- NA
  }
  destroy(rnb.set)  
  betas.table
}


read_geo_l1_data <- function(series_id_orig, targets, all.series.info, name) {
  cat('Reading ', series_id_orig, ": ")
  # handle samples which comes from multiple serieses
  series_id_vec <- unlist(strsplit(series_id_orig, ","))
  
  series_id <- NULL
  for(series_id_tmp in series_id_vec) {
    series_id_folder <- file.path(big_data_folder, "GEO", series_id_tmp)
    series_id_files <- list.files(series_id_folder, pattern="*.(txt|csv|tsv)$")
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
  
  if(length(grep("Signal_A.NA|_unmeth", series_id_files)) > 0 ) {
    # works for GSE62992
    # => two files of raw signals: signal A and signal B, no pvals
    signals <- lapply(series_id_fp, FUN=read_l1_signal_file)
    colnum <- length(colnames(signals[[1]]))
    samples.all <- gsub("[.]Signal_A","", colnames(signals[[1]]))
    
    if(length(samples.all) == dim(this_all.series.info)[[1]]) {
      v <- (this_targets$description %in% this_all.series.info$description) & (this_targets$source_name_ch1 %in% this_all.series.info$source_name_ch1)
      relevant.samples.loc <- which(v)
    } else {
      stop('try other option 3')
    }
    
    # assign  unmethylated, methylated and pvalue matrices
    U <- data.matrix(signals[[1]])
    colnames(U) <- samples.all
    M <- data.matrix(signals[[2]])
    colnames(M) <- samples.all
  } else {
    # works for GSE32079, GSE29290, GSE57767, GSE61653
    # => one raw file with 3 columns for each sample
    
    # GSE36278 has two raw files
    if(length(series_id_fp) == 2) {
      signals <- do.call("cbind", lapply(series_id_fp, FUN=read_l1_signal_file))
    } else {
      signals <- read_l1_signal_file(series_id_fp)
    }
    
    # locate relevant samples

    # remove suffixes from colnames
    unmeth_suffixes = c("[. _]?[Uu]nmethylated[. _][Ss]ignal$", "[_ ]{1,2}Unmethylated$",
                        "[.]Signal_A$", 
                        "_Unmethylated[.]Detection$")
    meth_suffixes = c("[. _]?[Mm]ethylated[. _][Ss]ignal$", "[_ ]{1,2}Methylated$",
                      "[.]Signal_B$", 
                      "_Methylated[.]Detection$")
    pvalue_suffixes = c("_[ ]?pValue$",
                        "[. _]Detection[. ]?Pval$", "[.]Pval$", "[.]Detection$",
                        "_detection_pvalue$")
    other_suffixes = c("_ M$")
    suffixes = c(unmeth_suffixes, meth_suffixes, pvalue_suffixes, other_suffixes)

    if(colnames(signals)[[1]] == "ID_Ref") {
      # GSE53162 which has two rownames columns
      rownames(signals) <- signals[, 1]
      signals <- signals[,-c(1)]
    } else if (colnames(signals)[[1]] == "TargetID" 
               && colnames(signals)[[2]] == "ProbeID_A" && colnames(signals)[[3]] == "ProbeID_B") {
      # GSE50874
      rownames(signals) <- signals[, 1]
      signals <- signals[,-c(1,2,3)]
      
    }
    
    # GSE52576 - has 5 columns per sample: AVG_Beta, Intensity
    # GSE50874 - has 4 columns per sample: AVG_Beta
    signals <- signals[!grepl("[.]AVG_Beta|[.]Intensity", colnames(signals))]
    colnum <- length(colnames(signals))
    orig <- colnames(signals)
    
    unmeth_ids = grepl(paste(unmeth_suffixes, collapse="|"), orig)
    meth_ids =  grepl(paste(meth_suffixes, collapse="|"), orig)
    pval_ids = grepl(paste(pvalue_suffixes, collapse="|"), orig)
    # check GSE47627, GSE42118

    colnames(signals) <- mgsub(suffixes, character(length(suffixes)), colnames(signals))
    print (colnames(signals))
    samples.all <- colnames(signals)[unmeth_ids]
    
    try_match_list <- list(this_targets$description, 
                       # first word
                       gsub("[ ;].*", '', this_targets$source_name_ch1), gsub("[ ;].*", '', this_targets$description), 
                       # last word
                       gsub(".* ", '', this_targets$source_name_ch1), gsub(".* ", '', this_targets$title), gsub(".*[ ;\t]", '', this_targets$description)
    )
    
    for(try_match in try_match_list) {
      relevant.samples.loc <- match(as.character(try_match), samples.all)  
      if(!all(is.na(relevant.samples.loc))) {
        break;
      }
    }
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
      # convert the decimal comma into a dot (as in GSE29290)
      signals_pval <- data.frame(lapply(signals_pval, function(x) gsub(",", ".", x, fixed = TRUE)), 
                                 row.names=rownames(signals_pval), stringsAsFactors=FALSE)
      p.values <- data.matrix(signals_pval)[,relevant.samples.loc, drop = FALSE]
    }
  }
  betas.table <- rnb_read_l1_betas(this_targets, U, M, p.values)
  write.table(betas.table, file.path(generated_GEO_folder, paste0(series_id, '_', name, '.txt')), 
              sep='\t', col.names=NA, quote=FALSE)
  
  stime <- (proc.time() - ptime1)[3]
  cat("   in", stime, "seconds\n")
}

work_on_targets <- function(targets, all.series.info) {
  series_id <- levels(factor(targets$series_id))
  name <- paste0(levels(factor(targets$disease)), ".", levels(factor(targets$tissue)))
  cat("Reading", nrow(targets), "samples of", name, "from", length(series_id), "serieses\n")
  ret <- lapply(series_id, FUN=read_geo_l1_data, targets, all.series.info, name)
}

dir.create(generated_GEO_folder, recursive=TRUE, showWarnings=FALSE)
folder <- file.path(data_folder, "global/GEO/joined")
joined_files <- list.files(folder, full.names = TRUE, pattern="*.txt")
joined_files <- joined_files[1:130]

# == skip serieses ==
# GEOs which I don't know how to parse:
# - no l1 signals txt file
# - different parsing on l1 txt file
bad_list <- c("GSE30338", "GSE37754", "GSE37965", "GSE39279", "GSE40360", 
              "GSE39560", "GSE40279", "GSE41826", "GSE41169", "GSE43976", 
              "GSE49377", "GSE48461", "GSE42882", "GSE45529", "GSE46573", 
              "GSE46306", "GSE48684", "GSE31803")
# GEOs which I still don't have
wait_list <- c()
# working GEOs
working_list <- c("GSE32079", "GSE38266", "GSE35069", "GSE32283", "GSE36278", 
                  "GSE29290", "GSE32146", "GSE37362", "GSE38268", "GSE40853", 
                  "GSE39958", "GSE41114", "GSE42372", "GSE41273", "GSE43091", 
                  "GSE44661", "GSE43293", "GSE43298", "GSE46394", "GSE44684",
                  "GSE42118", "GSE42119", "GSE43414", "GSE47512", "GSE42752",
                  "GSE45187", "GSE48325", "GSE49031", "GSE49656", "GSE49576", 
                  "GSE31803", "GSE49542", "GSE44667", "GSE45353", "GSE51758",
                  "GSE50498", "GSE50774", "GSE50759", "GSE53162",
                  "GSE52826", "GSE52576", "GSE50874")
ignore_list <- paste0("../../data/global/GEO/joined/", c(bad_list, wait_list, working_list), ".txt")
joined_files <- joined_files[!(joined_files %in% ignore_list)]

all.series.info <- do.call("rbind", lapply(joined_files, FUN=read_joined_file))
# Remove serieses with idats
all.series.info <- subset(all.series.info, is.na(supplementary_file))

# get only relevant samples
relevant.samples.idx <- which(as.numeric(all.series.info$relevant) == 1)
pheno <- all.series.info[relevant.samples.idx, ]
splited_targets <- split(pheno, list(pheno$disease, pheno$tissue), drop=TRUE)

ret <- lapply(splited_targets, FUN=work_on_targets, all.series.info)

write_nrow_per_group(splited_targets, file.path(generated_GEO_folder, 'GEO_all_kinds.csv'))
print("DONE")

# = TODO =

# why GSE35069 output doesn't have tissue name: GSE35069_Healthy..txt

# GSE29290 have same samples in GSE29290_Healthy.Breast.txt and GSE29290_Breast cancer.Breast.txt:
#   Sample_1	Sample_2	Sample_3	Sample_4	Sample_5	Sample_6	Sample_7	Sample_8

# GSE32283_Glioblastoma.Brain.txt has lots of NAs

# GSE38266 has hard columns names - should figure how to resolve them against joiner table

# GSE32146 last columns doesn't have name: 509 Unmethylated Signal, 509 Methylated Signal, Detection Pval

# GSE40360 gets: Error in read.table(filename, header = TRUE, row.names = 1, skip = 0,  : more columns than column names

# GSE40279 has different column names (4 per sample):
# "5815284007_R01C01.AVG_Beta"  "5815284007_R01C01.Intensity" "5815284007_R01C01.SignalA" "5815284007_R01C01.SignalB"

# GSE47627 has two samples: 
#  GSM1180517   B cells [CD19_Fer] (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1180517)
#  GSM1180518 	B cells [CD19_Javi] (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1180518)
# which have raw data on GEO site to download

# GSE46306
# GSE48684 (which also has two different parsing style files)
# has five columns: SignalA, SignalB, Intensity, AVG_Beta & Detection Pval
# moreover, it has first column of numeric indexing: 1, 2, 3 ...
