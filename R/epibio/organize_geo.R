
library(RnBeads)

source("config.R")
source("common.R")
#source("geo_l1_reader.R")

read_l1_signal_file <- function(filename) {
  t <- read.table(filename, header=TRUE, row.names=1, skip=0, sep='\t', dec = ".",
                  nrows=200, 
                  check.names=FALSE, stringsAsFactors=FALSE)
  t
}

read_joined_file <- function(filename) {
  t <- read.table(filename, sep='\t', header=TRUE, row.names=1, fill=TRUE, 
                  na.strings=c("NA", "0"), quote="\"", stringsAsFactors=FALSE)
  df <- data.frame(Filename=filename, t)
  df
}

rnb_read_l1_betas <- function(targets, U, M, p.values) {
  pheno <- targets[, c('description','tissue','cell_type','disease')]
  rnb.set <- new('RnBeadRawSet', pheno, U=U, M=M, p.values=p.values, useff=FALSE)
  
  logger.start(fname=NA)
  rnb.options(disk.dump.big.matrices=TRUE)
  rnb.options(enforce.memory.management=TRUE)
  rnb.set <- rnb.execute.snp.removal(rnb.set)$dataset
  
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
  # when sample is from multiple serieses - use the first only
  series_id <- sub(",.*", "", series_id_orig)
  this_targets = subset(targets, targets$series_id == series_id_orig)  
  series_id_folder <- file.path(big_data_folder, "GEO", series_id)
  series_id_files <- list.files(series_id_folder, pattern="*.txt$")
  filename_first_level <- levels(factor(this_targets$Filename))[[1]]
  this_all.series.info <- subset(all.series.info, all.series.info$Filename == filename_first_level)
  
  if(length(series_id_files) == 0) {
    msg <- paste('no txt files inside', series_id_folder)
    stop(msg)
  } else {
    ptime1 <- proc.time()
    print(series_id_files)
    series_id_fp <- file.path(series_id_folder, series_id_files)
    if(length(grep("Signal_A.NA", series_id_files)) > 0 ) {
      # works for GSE62992
      # => two files of raw signals, signal A and signal B, no pvals      
      signals <- lapply(series_id_fp, FUN=read_l1_signal_file)
      colnum <- length(colnames(signals[[1]]))
      samples.all <- gsub("[.]Signal_A","", colnames(signals[[1]]))
      
      if(length(samples.all) == dim(this_all.series.info)[[1]]) {
        v <- (this_targets$description %in% this_all.series.info$description) & (this_targets$source_name_ch1 %in% this_all.series.info$source_name_ch1)
        relevant.samples.loc <- which(v)
      } else {
        stop('try other option')
      }
      
      # assign  unmethylated, methylated and pvalue matrices
      U <- data.matrix(signals[[1]])
      colnames(U) <- samples.all
      M <- data.matrix(signals[[2]])
      colnames(M) <- samples.all
      p.values <- NULL
    } else {
      # works for GSE32079, GSE29290, GSE57767, GSE61653
      # => one raw file with 3 columns for each sample
      signals <- read_l1_signal_file(series_id_fp)
      # locate relevant samples
      colnum <- length(colnames(signals))
      unmeth_ids = seq(1, colnum-2, 3)
      meth_ids = seq(2, colnum-1, 3)
      pval_ids = seq(3, colnum, 3)
      
      # remove suffixes from colnames
      suffixes = c("[.]Unmethylated[. ]Signal$", "[.]Methylated[. ]Signal$", 
                   "_Methylated signal$",
                   "[.]Signal_A$", "[.]Signal_B$", 
                   "_Unmethylated[.]Detection$", "_Methylated[.]Detection$",
                   "[.]Detection[.]Pval$", "[.]Detection Pval", "[.]Pval$", "[.]Detection$")
      colnames(signals) <- mgsub(suffixes, character(length(suffixes)), colnames(signals))
      samples.all <- colnames(signals)[unmeth_ids]
      relevant.samples.loc <- match(as.character(this_targets$description), samples.all)
      if(all(is.na(relevant.samples.loc))) {
        first_word <- gsub(" .*", '', this_targets$source_name_ch1)
        relevant.samples.loc <- match(as.character(first_word), samples.all)
        if(all(is.na(relevant.samples.loc))) {
          if(length(samples.all) == dim(this_all.series.info)[[1]]) {
            v <- (this_targets$description %in% this_all.series.info$description) & (this_targets$source_name_ch1 %in% this_all.series.info$source_name_ch1)
            relevant.samples.loc <- which(v)
          } else {
            stop('try other option')
          }
        }
      }
      
      # assign  unmethylated, methylated and pvalue matrices
      U <- data.matrix(signals[,unmeth_ids])[,relevant.samples.loc]
      M <- data.matrix(signals[,meth_ids])[,relevant.samples.loc]
      signals_pval <- signals[,pval_ids]
      # convert the decimal comma into a dot (as in GSE29290)
      signals_pval <- data.frame(lapply(signals_pval, 
                                        function(x) gsub(",", ".", x, fixed = TRUE)), 
                                        row.names=rownames(signals_pval),
                                        stringsAsFactors=FALSE)
      p.values <- data.matrix(signals_pval)[,relevant.samples.loc]
    }
    betas.table <- rnb_read_l1_betas(this_targets, U, M, p.values)
    write.table(betas.table, file.path(generated_GEO_folder, paste0(series_id, '_', name, '.txt')), 
                sep='\t', col.names=NA, quote=FALSE)
    
    stime <- (proc.time() - ptime1)[3]
    cat(" in", stime, "seconds\n")
  }
}

work_on_targets <- function(targets, all.series.info) {
  series_id <- levels(factor(targets$series_id))
  name <- paste0(levels(factor(targets$disease)), ".", levels(factor(targets$tissue)))
  cat("Reading", nrow(targets), "samples of", name, "from", length(series_id), "serieses\n")
  ret <- lapply(series_id, FUN=read_geo_l1_data, targets, all.series.info, name)
}

dir.create(generated_GEO_folder, recursive=TRUE, showWarnings=FALSE)
folder <- file.path(data_folder, "global/GEO/joined")
joined_files <- file.path(folder, list.files(folder, pattern="*.txt"))
joined_files <- head(joined_files, 40)
all.series.info <- do.call("rbind", lapply(joined_files, FUN=read_joined_file))
# get only relevant samples
relevant.samples.idx <- which(as.numeric(all.series.info$relevant) == 1)
pheno <- all.series.info[relevant.samples.idx, ]
splited_targets <- split(pheno, list(pheno$disease, pheno$tissue), drop=TRUE)

ret <- lapply(splited_targets, FUN=work_on_targets, all.series.info)

write_nrow_per_group(splited_targets, file.path(generated_GEO_folder, 'GEO_all_kinds.csv'))
print("DONE")
