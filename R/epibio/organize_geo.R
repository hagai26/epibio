
library(RnBeads)

source("config.R")
source("common.R")


work_on_targets <- function(targets) {
  # work on thses targets using RnBeads
  print(targets)
}


dir.create(generated_GEO_folder, recursive=TRUE, showWarnings=FALSE)
folder <- file.path(data_folder, "global/GEO/joined")
joined_files <- file.path(folder, list.files(folder, pattern="*.txt"))

f1 <- "../../data/global/GEO/joined/GSE32079.txt"
f2 <- "../../data/global/GEO/joined/GSE32148.txt" # bad - more columns than column names
f3 <- "../../data/global/GEO/joined/GSE62992.txt"
f4 <- "../../data/global/GEO/joined/GSE62640.txt"
f5 <- "../../data/global/GEO/joined/GSE61653.txt"
f6 <- "../../data/global/GEO/joined/GSE61446.txt"
f7 <- "../../data/global/GEO/joined/GSE57894.txt"
f8 <- "../../data/global/GEO/joined/GSE57767.txt"
f9 <- "../../data/global/GEO/joined/GSE32146.txt"
f10 <- "../../data/global/GEO/joined/GSE30870.txt"
f11 <- "../../data/global/GEO/joined/GSE29290.txt"
joined_files <- c(f1, f3, f4, f5, f6, f7, f8, f10, f11)
series.info <- do.call("rbind", lapply(joined_files, function(fn) 
  data.frame(Filename=fn, read.table(fn, sep='\t', row.names=1, header=TRUE))
  ))
# get only relevant samples
relevant.samples.idx <- as.numeric(series.info$relevant)
relevant.samples.idx[is.na(relevant.samples.idx)] <- 0
relevant.samples.idx <- relevant.samples.idx == 1

pheno <- series.info[relevant.samples.idx, c('description','tissue','cell_type','disease')]
num_samples <- length(series.info[,1])

result <- chunked_group_by(pheno, list(pheno$disease, pheno$tissue), 2)

all_kinds = data.frame(number = sapply(result$splited, FUN=nrow))
all_kinds_filename <- file.path(generated_GEO_folder, 'GEO_all_kinds.csv')
write.csv(cbind(kind=rownames(all_kinds), all_kinds), file = all_kinds_filename, row.names=FALSE, quote=FALSE)

ret <- lapply(result$grouped, FUN=work_on_targets)
print("DONE")
