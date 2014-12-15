
library(RnBeads)

source("config.R")
source("common.R")

read_geo_l1_data <- function(series_id) {
  cat('Reading ', series_id, ": ")
  series_id_folder <- file.path(big_data_folder, "GEO", series_id)
  series_id_files <- list.files(series_id_folder, pattern="*.txt$")
  if(length(series_id_files) == 0) {
    print('not found')
  } else {
    print(series_id_files)
    series_id_fp <- file.path(series_id_folder, series_id_files[[1]])
    print(series_id_fp)
    series_id_fp <- "../../data/big/GEO/GSE32079/GSE32079_non-normalized_small.txt"
    print(series_id_fp)
    #signals <- read.table(series_id_fp, nrows=-1, header=TRUE, row.names=1, skip=0, sep='\t', dec = ".")
    
    # locate relevant samples
    #print(colnames(signals))
    #samples.all <- gsub(".Signal_A","", colnames(signals)[seq(1,(3*num_samples-2),3)])
    #relevant.samples.loc <- match(as.character(pheno$description),samples.all)
    
    #remove suffixes from colnames
    #colnames(signals) <- gsub(".Signal_A","",colnames(signals))
    #colnames(signals) <- gsub(".Signal_B","",colnames(signals))
    #colnames(signals) <- gsub(".Detection","",colnames(signals))
    
    # assign  unmethylated, methylated and pvalue matrices
    #U <- data.matrix(signals[,seq(1,(3*num_samples-2),3)])[,relevant.samples.loc]
    #M <- data.matrix(signals[,seq(2,(3*num_samples-1),3)])[,relevant.samples.loc]
    #p.values <- data.matrix(signals[,seq(3,(3*num_samples),3)])[,relevant.samples.loc]
    
    #run rnbeads preprecossing
    #rnb.raw.set <- new('RnBeadRawSet',pheno,U=U,M=M,p.values=p.values,useff=FALSE)
    
    #logger.start(fname=NA)
    #rnb.raw.set.greedy <- rnb.execute.greedycut(rnb.raw.set)
    #rnb.raw.set.greedy.snprem <- rnb.execute.snp.removal(rnb.raw.set)$dataset
    #rnb.set.norm <- rnb.execute.normalization(rnb.raw.set.greedy.snprem,method="bmiq",bgcorr.method="methylumi.lumi")
    #rnb.set.sexrem <- rnb.execute.sex.removal(rnb.set.norm)$dataset
    
    #meth.beta <- meth(rnb.set.sexrem)
  }  
}

print('start')
series_id <- "GSE32079"
read_geo_l1_data(series_id)
print('done')
