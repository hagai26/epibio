
library(RnBeads)

source("config.R")
source("common.R")


read_geo_l1_data <- function(one_series_id, relevant_targets) {
  cat('Reading ', series_id, ": ")
  series_id_folder <- file.path(big_data_folder, "GEO", series_id)
  series_id_files <- list.files(series_id_folder, pattern="*.txt$")
  if(length(series_id_files) == 0) {
    print('not found')
  } else {
    print(series_id_files)
    series_id_fp <- file.path(series_id_folder, series_id_files[[1]])
    print(series_id_fp)
    signals <- read.table(series_id_fp, nrows=-1, header=TRUE, row.names=1, skip=0, sep='\t', dec = ".")
    
    # locate relevant samples
    colnum <- length(colnames(signals))
    samples.all <- gsub(".Signal_A","", colnames(signals)[seq(1, (colnum-2), 3)])
    relevant.samples.loc <- c("Q09_MayoCC027T", "Q10_MayoCC049T")
    
    # remove suffixes from colnames
    colnames(signals) <- gsub(".Signal_A", "", colnames(signals))
    colnames(signals) <- gsub(".Signal_B", "", colnames(signals))
    colnames(signals) <- gsub(".Detection.Pval", "", colnames(signals))
    
    # assign  unmethylated, methylated and pvalue matrices
    unmeth_ids = seq(1, colnum-2, 3)
    meth_ids = seq(2, colnum-1, 3)
    pval_ids = seq(3, colnum, 3)
    U <- data.matrix(signals[,unmeth_ids])[,relevant.samples.loc]
    M <- data.matrix(signals[,meth_ids])[,relevant.samples.loc]
    p.values <- data.matrix(signals[,pval_ids])[,relevant.samples.loc]
    
    # run rnbeads preprecossing
    rnb.raw.set <- new('RnBeadRawSet', pheno, U=U, M=M, p.values=p.values, useff=FALSE)
    
    #logger.start(fname=NA)
    #rnb.raw.set.greedy <- rnb.execute.greedycut(rnb.raw.set)
    #rnb.raw.set.greedy.snprem <- rnb.execute.snp.removal(rnb.raw.set)$dataset
    #rnb.set.norm <- rnb.execute.normalization(rnb.raw.set.greedy.snprem,method="bmiq",bgcorr.method="methylumi.lumi")
    #rnb.set.sexrem <- rnb.execute.sex.removal(rnb.set.norm)$dataset
    
    #meth.beta <- meth(rnb.set.sexrem)
  }
}

