## saved summarized files to smaller ones, using these methods:
##  - filter some data
##  - round floats

library(data.table)

source("config.R")
source("common.R")


make_summary_smaller <- function(folder_name) {
  print(sprintf('- Making %s smaller -', folder_name))
  input_folder <- file.path(summary_folder, folder_name)
  output_folder <- file.path(summary_folder, folder_name)
  
  csv_filename <- file.path(input_folder, 'means.txt')
  data <- read.table(csv_filename, sep='\t', header = TRUE)
  data[,-1] <- round(data[,-1], 2) # the "-1" excludes column 1
  write.csv(data, file.path(output_folder, 'means_small.csv'), row.names=FALSE, quote=FALSE)
}

make_summary_smaller('GEO')
make_summary_smaller('TCGA')
make_summary_smaller('lab_data')

print("DONE")
