## summarize merged gz files files to 3 files:
## mean, 10th percentile, 90th percentile

library(data.table)

source("config.R")
source("common.R")

read_merged_table <- function(filepath, col_name, index, num.files) {
  print(sprintf('Summarizing %s %d/%d', col_name, index, num.files))
  data <- data.table(read.table(filepath, sep='\t', header=TRUE, row.names=1), keep.rownames=TRUE, key='rn')
  setnames(data, c("mean", "quantile_0.1", "quantile_0.9"), rep(col_name, 3))
  data
}

summary_files <- function(input_folder, output_folder) {
  mean.col <- 2
  tenth.col <- 4
  ninetieth.col <- 5
  
  gz.files <- list.files(input_folder, full.names=TRUE)
  gz.files_basenames <- basename(gz.files)
  new.col.names <- substr(gz.files_basenames, 1, nchar(gz.files_basenames) - nchar('.txt.gz'))
  num.files <- length(gz.files)
  # initiate variables using first file
  data <- read_merged_table(gz.files[1], new.col.names[1], 1, num.files)
  
  means <- data[, c(1, mean.col), with=FALSE]
  perc.10 <- data[, c(1, tenth.col), with=FALSE]
  perc.90 <- data[, c(1, ninetieth.col), with=FALSE]
  
  for (i in 2:num.files) {
    data <- read_merged_table(gz.files[i], new.col.names[i], i, num.files)
    
    means.new <- data[, c(1, mean.col), with = FALSE]
    perc.10.new <- data[, c(1, tenth.col), with = FALSE]
    perc.90.new <- data[, c(1, ninetieth.col), with = FALSE]
    
    means <- merge(means, means.new, all=TRUE)
    perc.10 <- merge(perc.10, perc.10.new, all=TRUE)
    perc.90 <- merge(perc.90, perc.90.new, all=TRUE)
  }
  write.table(means, file.path(output_folder,'means.txt'), sep = '\t', col.names=TRUE, row.names=FALSE)
  write.table(perc.10, file.path(output_folder,'perc_10.txt'), sep = '\t', col.names=TRUE, row.names=FALSE)
  write.table(perc.90, file.path(output_folder,'perc_90.txt'), sep = '\t', col.names=TRUE, row.names=FALSE)
}


summary_merged_folder <- function(folder_name) {
  print(sprintf('- Summarizing %s -', folder_name))
  input_folder <- file.path(generated_merged_folder, folder_name)
  output_folder <- file.path(summary_folder, folder_name)
  
  dir.create(output_folder, showWarnings = FALSE)
  summary_files(input_folder, output_folder)  
}

dir.create(file.path(summary_folder), showWarnings = FALSE)

summary_merged_folder('GEO')
summary_merged_folder('TCGA')
summary_merged_folder('lab_data')

print("DONE")
