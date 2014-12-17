
library(RnBeads)

source("config.R")
source("common.R")
#source("geo_l1_reader.R")


read_geo_l1_data <- function(series_id, targets, type, study) {
  cat('Reading ', series_id, ": ")
  series_id_folder <- file.path(big_data_folder, "GEO", series_id)
  series_id_files <- list.files(series_id_folder, pattern="*_small.txt$") # XXX
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
      signals <- lapply(series_id_fp, function(x) read.table(x, nrows=-1, header=TRUE, row.names=1, skip=0, sep='\t', dec = "."))
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
      signals <- read.table(series_id_fp, nrows=-1, 
                            header=TRUE, row.names=1, skip=0, sep='\t', dec = ".",
                            check.names=FALSE)
      
      # locate relevant samples
      colnum <- length(colnames(signals))
      unmeth_ids = seq(1, colnum-2, 3)
      meth_ids = seq(2, colnum-1, 3)
      pval_ids = seq(3, colnum, 3)
      
      # remove suffixes from colnames
      colnames(signals) <- gsub(".Signal_A", "", colnames(signals))
      colnames(signals) <- gsub(".Signal_B", "", colnames(signals))
      colnames(signals) <- gsub(".Unmethylated.Signal", "", colnames(signals))
      colnames(signals) <- gsub(".Methylated.Signal", "", colnames(signals))
      colnames(signals) <- gsub(".Detection.Pval", "", colnames(signals))
      colnames(signals) <- gsub(".Pval", "", colnames(signals))
      
      samples.all <- colnames(signals)[unmeth_ids]
      relevant.samples.loc <- match(as.character(targets$description), samples.all)
      if(all(is.na(relevant.samples.loc))) {
        print('trying first word as sample name')
        first_word <- gsub(" .*", '', targets$source_name_ch1)
        relevant.samples.loc <- match(as.character(first_word), samples.all)
        if(all(is.na(relevant.samples.loc))) {
          print('try other option')
        }
      }
      
      # assign  unmethylated, methylated and pvalue matrices
      U <- data.matrix(signals[,unmeth_ids])[,relevant.samples.loc]
      M <- data.matrix(signals[,meth_ids])[,relevant.samples.loc]
      p.values <- data.matrix(signals[,pval_ids])[,relevant.samples.loc]
    }
    # run rnbeads preprecossing
    pheno <- targets[, c('description','tissue','cell_type','disease')]
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
    filename <- file.path(generated_GEO_folder, paste0(series_id, '_', type, '_', study, '.txt'))
    write.table(betas.table, filename, sep='\t', col.names=NA, quote=FALSE)
    #rnb.execute.export.csv(rnb.set.sexrem, NA)
    
    stime <- (proc.time() - ptime1)[3]
    cat(" in", stime, "seconds\n")
  }
}


work_on_targets <- function(targets) {
  print("work_on_targets called")
  ptime1 <- proc.time()
  series_id <- levels(factor(targets$series_id))
  cat("Reading", nrow(targets), "samples", "from", length(series_id), "serieses")
  print("")
  # when sample is from multiple serieses - use the first only
  series_id <- sub(",.*", "", series_id) 
  ret <- lapply(series_id, FUN=read_geo_l1_data, targets, type, study)
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
joined_files <- c(f8, f10, f11)
#joined_files <- head(joined_files, 4) # XXX
series.info <- do.call("rbind", lapply(joined_files, FUN=read_joined_file))
# get only relevant samples
relevant.samples.idx <- which(as.numeric(series.info$relevant) == 1)
pheno <- series.info[relevant.samples.idx, c('series_id', 'tissue', 'cell_type', 'disease')]
splited_targets <- split(pheno, list(pheno$disease, pheno$tissue), drop=TRUE)

write_nrow_per_group(splited_targets, file.path(generated_GEO_folder, 'GEO_all_kinds.csv'))

ret <- lapply(splited_targets, FUN=work_on_targets)
print("DONE")
