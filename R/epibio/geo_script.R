
# folders
data_folder <- "..\\..\\data"
generated_folder  <- file.path(data_folder, "generated")
geo_folder <- file.path(generated_folder, "GEO")
txt_folder <- file.path(data_folder, "static\\GSE32148_matrix_signal_peripheralBlood.txt")
dir.create(geo_folder, recursive=TRUE, showWarnings=FALSE)

table_filename <- file.path(txt_folder, "GSE32148_matrix_signal_peripheralBlood_small.txt")
targets <- read.table(table_filename, sep="\t", header=TRUE, row.names=1)

relevant_target_names <- c("X10B", "X264")

names_len <- length(names(targets))
unmeth_ids = seq(1, names_len, 3)
meth_ids = seq(2, names_len, 3)
pval_ids = seq(3, names_len, 3)


splited_names <- strsplit(names(targets), '[.]')
target_names <- sapply(splited_names, function(x) x[[1]])

  
#MethylSet(targets$)