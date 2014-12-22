
library(RnBeads)

source("config.R")
source("common.R")
#source("geo_l1_reader.R")

read_l1_signal_file <- function(filename) {
  t <- read.table(filename, nrows=200, 
             header=TRUE, row.names=1, skip=0, sep='\t', dec = ".",
             check.names=FALSE)
  return(t)
}


read_geo_l1_data <- function(series_id_orig, targets, all.series.info) {
  cat('Reading ', series_id_orig, ": ")
  # when sample is from multiple serieses - use the first only
  series_id <- sub(",.*", "", series_id_orig)
  this_targets = subset(targets, targets$series_id == series_id_orig)
  
  series_id_folder <- file.path(big_data_folder, "GEO", series_id)
  series_id_files <- list.files(series_id_folder, pattern="*.txt$")
  if(length(series_id_files) == 0) {
    print('not found')
  } else {
    ptime1 <- proc.time()
    print(series_id_files)
    if(length(grep("Signal_A.NA", series_id_files)) > 0 ) {
      # works for GSE62992
      # two files of raw signals, signal A and signal B, no pvals
      series_id_fp <- file.path(series_id_folder, series_id_files)
      print(series_id_fp)
      signals <- lapply(series_id_fp, FUN=read_l1_signal_file)
      colnum <- length(colnames(signals[[1]]))
      samples.all <- gsub(".Signal_A","", colnames(signals[[1]]))
      #relevant.samples.loc <- match(as.character(rownames(relevant_targets)), samples.all)
      # should fix the description vs samples names
      
      # assign  unmethylated, methylated and pvalue matrices
      U <- data.matrix(signals[[1]])
      colnames(U) <- samples.all
      M <- data.matrix(signals[[2]])
      colnames(M) <- samples.all
      p.values <- NULL
      
    } else {
      # works for GSE32079
      # one raw file with 3 columns for each sample
      series_id_fp <- file.path(series_id_folder, series_id_files[[1]])
      print(series_id_fp)
      signals <- read_l1_signal_file(series_id_fp)
      
      # locate relevant samples
      colnum <- length(colnames(signals))
      unmeth_ids = seq(1, colnum-2, 3)
      meth_ids = seq(2, colnum-1, 3)
      pval_ids = seq(3, colnum, 3)
      
      # remove suffixes from colnames
      suffixes = c("[.]Unmethylated[.]Signal$", "[.]Methylated[.]Signal$", 
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
          fn <- levels(factor(this_targets$Filename))[[1]]
          this_all.series.info <- subset(all.series.info, all.series.info$Filename == fn)
          if(length(samples.all) == dim(this_all.series.info)[[1]]) {
            relevant.samples.loc <- match(as.character(this_targets$description), this_all.series.info$description)
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
                                        stringsAsFactors=FALSE)
      p.values <- data.matrix(signals_pval)[,relevant.samples.loc]
    }
    # run rnbeads preprecossing
    pheno <- this_targets[, c('description','tissue','cell_type','disease')]
    rnb.set <- new('RnBeadRawSet', pheno, U=U, M=M, p.values=p.values, useff=FALSE)
    
    logger.start(fname=NA)
    rnb.options(disk.dump.big.matrices=TRUE)
    rnb.options(enforce.memory.management=TRUE)
    #rnb.set <- rnb.execute.greedycut(rnb.set)
    rnb.set <- rnb.execute.snp.removal(rnb.set)$dataset
    #rnb.set <- rnb.execute.normalization(rnb.set, 
    #                                          method="bmiq", 
    #                                          bgcorr.method="methylumi.lumi")
    #rnb.set <- rnb.execute.sex.removal(rnb.set)$dataset
    
    betas.table <- meth(rnb.set, row.names=TRUE)
    pvalue.high <- which(dpval(rnb.set)>0.05, arr.ind=TRUE)
    betas.table[pvalue.high[,'row'], pvalue.high[,'col']] <- NA
    destroy(rnb.set)
    filename <- file.path(generated_GEO_folder, paste0(series_id, '_', "type", '_', "study", '.txt'))
    write.table(betas.table, filename, sep='\t', col.names=NA, quote=FALSE)
    #rnb.execute.export.csv(rnb.set.sexrem, NA)
    
    stime <- (proc.time() - ptime1)[3]
    cat(" in", stime, "seconds\n")
  }
}


work_on_targets <- function(targets, all.series.info) {
  series_id <- levels(factor(targets$series_id))
  cat("Reading", nrow(targets), "samples", "from", length(series_id), "serieses\n")
  ret <- lapply(series_id, FUN=read_geo_l1_data, targets, all.series.info)
}

read_joined_file <- function(filename) {
  t <- read.table(filename, sep='\t', row.names=1, header=TRUE, fill=TRUE, 
                  na.strings=c("NA", "0"), quote="\"", stringsAsFactors=FALSE)
  df <- data.frame(Filename=filename, t)
  return(df)
}


dir.create(generated_GEO_folder, recursive=TRUE, showWarnings=FALSE)
folder <- file.path(data_folder, "global/GEO/joined")
joined_files <- file.path(folder, list.files(folder, pattern="*.txt"))

f1 <- "../../data/global/GEO/joined/GSE32079.txt"
f2 <- "../../data/global/GEO/joined/GSE32148.txt"
f3 <- "../../data/global/GEO/joined/GSE62992.txt"
f4 <- "../../data/global/GEO/joined/GSE62640.txt"
f5 <- "../../data/global/GEO/joined/GSE61653.txt"
f6 <- "../../data/global/GEO/joined/GSE61446.txt"
f7 <- "../../data/global/GEO/joined/GSE57894.txt"
f8 <- "../../data/global/GEO/joined/GSE57767.txt"
f9 <- "../../data/global/GEO/joined/GSE32146.txt"
f10 <- "../../data/global/GEO/joined/GSE30870.txt"
f11 <- "../../data/global/GEO/joined/GSE29290.txt"
#joined_files <- c(f1, f2, f3, f4, f5, f6, f7, f8, f10, f11)
joined_files <- c(f11)
all.series.info <- do.call("rbind", lapply(joined_files, FUN=read_joined_file))
# get only relevant samples
relevant.samples.idx <- which(as.numeric(all.series.info$relevant) == 1)
pheno <- all.series.info[relevant.samples.idx, ]
splited_targets <- split(pheno, list(pheno$disease, pheno$tissue), drop=TRUE)

ret <- lapply(splited_targets, FUN=work_on_targets, all.series.info)

write_nrow_per_group(splited_targets, file.path(generated_GEO_folder, 'GEO_all_kinds.csv'))
print("DONE")
