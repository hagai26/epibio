## summarize merged gz files files to 3 files:
## mean, std, 10th percentile, 90th percentile

library(data.table)

source("config.R")
source("common.R")

summary_files <- function(input_folder, output_folder) {
  mean.col <- 2
  std.col <- 3
  tenth.col <- 4
  ninetieth.col <- 5
  
  gz.files <- list.files(input_folder)
  new.col.names <- substr(basename(gz.files), 1, nchar(gz.files) - nchar('.txt.gz'))
  num.files <- length(gz.files)
  # initiate variables using first file
  print(paste('Summarizing ', new.col.names[1], ' ', 1, '/', num.files, sep=''))
  data.file <- file.path(input_folder, gz.files[1])
  data <-
    data.table(
      read.table(data.file, sep='\t', header=TRUE, row.names=1), 
      keep.rownames=TRUE, key='rn')
  
  setnames(data, c("mean", "quantile_0.1", "quantile_0.9"), rep(new.col.names[1], 3))
  means <- data[, c(1,mean.col), with=FALSE]
  perc.10 <- data[, c(1,tenth.col), with=FALSE]
  perc.90 <- data[, c(1,ninetieth.col), with=FALSE]
  
  for (i in 2:num.files) {
    print(paste('Summarizing ', new.col.names[i],' ',i,'/',num.files,sep = ''))
    data.file <- file.path(input_folder,gz.files[i])
    data <-
      data.table(
        read.table(
          data.file,sep = '\t',header = T,row.names = 1
        ),keep.rownames = T,key = 'rn'
      )
    setnames(data, c("mean","quantile_0.1","quantile_0.9"),rep(new.col.names[i],3))
    means.new <- data[, c(1,mean.col), with = FALSE]
    perc.10.new <- data[, c(1, tenth.col),with = FALSE]
    perc.90.new <- data[, c(1, ninetieth.col),with = FALSE]
    means <- merge(means, means.new, all=TRUE)
    perc.10 <- merge(perc.10, perc.10.new, all=TRUE)
    perc.90 <- merge(perc.90, perc.90.new, all=TRUE)
  }
  write.table(means, file.path(output_folder,'means.txt'), sep = '\t', col.names=TRUE, row.names=FALSE)
  write.table(perc.10, file.path(output_folder,'perc_10.txt'), sep = '\t', col.names=TRUE, row.names=FALSE)
  write.table(perc.90, file.path(output_folder,'perc_90.txt'), sep = '\t', col.names=TRUE, row.names=FALSE)
}

summary_merged_folder <- function(folder_name) {
  print(sprintf('Summarizing %s', folder_name))
  input_folder <- file.path(generated_merged_folder, folder_name)
  output_folder <- file.path(summary_folder, folder_name)
  dir.create(output_folder, showWarnings = FALSE)
  summary_files(input_folder, output_folder)  
}

dir.create(file.path(summary_folder), showWarnings = FALSE)

#summary_merged_folder('GEO')
#summary_merged_folder('TCGA')
summary_merged_folder('lab_data')

print("DONE")
