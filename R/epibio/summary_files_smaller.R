## saved summarized files to smaller ones, using these methods:
##  - filter some data
##  - round floats

library(data.table)

source("config.R")
source("common.R")


make_summary_smaller <- function(folder_name, good_blood_rows) {
  print(sprintf('- Making %s smaller -', folder_name))
  input_folder <- file.path(summary_folder, folder_name)
  output_folder <- file.path(summary_folder, folder_name)
  
  csv_filename <- file.path(input_folder, 'means.txt')
  data <- read.table(csv_filename, sep='\t', header = TRUE)
  old_row_num <- nrow(data)
  if(is.null(good_blood_rows)) {
    good_blood_rows <- data[!is.na(data$healthy.whole_blood) & data$healthy.whole_blood>=0.7,]$rn
  }
  striped_data <- data[good_blood_rows,]
  print(sprintf('new nrow=%d, old nrow=%d', nrow(striped_data), old_row_num))
  striped_data[,-1] <- round(striped_data[,-1], 2) # the "-1" excludes column 1
  write.csv(striped_data, file.path(output_folder, 'means_small.csv'), 
            row.names=FALSE, quote=FALSE)
  good_blood_rows
}

# we use GEO to get good_blood_rows
good_blood_rows <- make_summary_smaller('GEO', NULL)
make_summary_smaller('TCGA', good_blood_rows)
make_summary_smaller('lab_data', good_blood_rows)

print("DONE")
