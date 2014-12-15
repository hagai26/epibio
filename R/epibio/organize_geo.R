
library(RnBeads)

source("config.R")
source("common.R")
#source("geo_l1_reader.R")


read_geo_l1_data <- function(series_id, relevant_targets) {
  cat('Reading ', series_id, ": ")
  series_id_folder <- file.path(big_data_folder, "GEO", series_id)
  series_id_files <- list.files(series_id_folder, pattern="*.txt$")
  if(length(series_id_files) == 0) {
    print('not found')
  } else {
    ptime1 <- proc.time()
    print(series_id_files)
    series_id_fp <- file.path(series_id_folder, series_id_files[[1]])
    series_id_fp <- "../../data/big/GEO/GSE32079/GSE32079_non-normalized_small.txt"
    print(series_id_fp)
    signals <- read.table(series_id_fp, nrows=-1, header=TRUE, row.names=1, skip=0, sep='\t', dec = ".")
    
    # locate relevant samples
    colnum <- length(colnames(signals))
    samples.all <- gsub(".Signal_A","", colnames(signals)[seq(1, (colnum-2), 3)])
    relevant.samples.loc <- match(as.character(relevant_targets$description), samples.all)
    
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
    pheno <- relevant_targets[, c('description','tissue','cell_type','disease')]
    rnb.raw.set <- new('RnBeadRawSet', pheno, U=U, M=M, p.values=p.values, useff=FALSE)
    
    logger.start(fname=NA)
    rnb.raw.set.greedy <- rnb.execute.greedycut(rnb.raw.set)
    rnb.raw.set.greedy.snprem <- rnb.execute.snp.removal(rnb.raw.set)$dataset
    #rnb.set.norm <- rnb.execute.normalization(rnb.raw.set.greedy.snprem, 
    #                                          method="bmiq", 
    #                                          bgcorr.method="methylumi.lumi")
    #rnb.set.sexrem <- rnb.execute.sex.removal(rnb.set.norm)$dataset
    rnb.set.sexrem <- rnb.raw.set.greedy.snprem
    
    meth.beta <- meth(rnb.set.sexrem)
    print(head(meth.beta))
    stime <- (proc.time() - ptime1)[3]
    cat(" in", stime, "seconds\n")
  }
}


work_on_targets <- function(targets) {
  print("work_on_targets called")
  study_levels <- levels(factor(targets$disease))
  type_levels <- levels(factor(targets$tissue))
  
  ptime1 <- proc.time()
  cat("Reading", nrow(targets), "samples")
  print("")
  # work on these targets
  for(type in type_levels) {
    for(study in study_levels) {
      is_relevant_targets <- targets$disease==study & targets$tissue==type
      if(sum(is_relevant_targets) > 0) {
        relevant_targets <- targets[is_relevant_targets,, drop = FALSE]
        cat('currently reading:', type, study, "(", sum(is_relevant_targets), "samples) :")
        series_id <- levels(factor(targets$series_id))
        series_id <- sub(",.*", "", series_id) # read each sample only once
        for(one_series_id in series_id) {
          read_geo_l1_data(one_series_id, relevant_targets)
        }
      }
      
    }
  }
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
joined_files <- c(f1, f2, f3, f4, f5, f6, f7, f8, f10, f11)
joined_files <- head(joined_files, 1)
series.info <- do.call("rbind", lapply(joined_files, function(fn) 
  data.frame(
    Filename=fn, 
    read.table(fn, sep='\t', row.names=1, header=TRUE, fill=TRUE, na.strings=c("NA", "0"), quote="\"", stringsAsFactors=FALSE)
    )
  ))
# get only relevant samples
relevant.samples.idx <- as.numeric(series.info$relevant)
relevant.samples.idx[is.na(relevant.samples.idx)] <- 0
relevant.samples.idx <- relevant.samples.idx == 1

pheno <- series.info[relevant.samples.idx, c('series_id', 'description','tissue','cell_type','disease')]
num_samples <- length(series.info[,1])

result <- chunked_group_by(pheno, list(pheno$disease, pheno$tissue), 3)

all_kinds = data.frame(number = sapply(result$splited, FUN=nrow))
all_kinds_filename <- file.path(generated_GEO_folder, 'GEO_all_kinds.csv')
write.csv(cbind(kind=rownames(all_kinds), all_kinds), file = all_kinds_filename, row.names=FALSE, quote=FALSE)

ret <- lapply(result$grouped, FUN=work_on_targets)
print("DONE")
