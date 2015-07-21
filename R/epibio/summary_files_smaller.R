
library(data.table)
#library(openxlsx)

source("config.R")
source("common.R")

folder_name <- 'lab_data'
input_folder <- file.path(summary_folder, folder_name)
output_folder <- file.path(summary_folder, folder_name)

csv_filename <- file.path(input_folder, 'means.txt')
data <- read.table(csv_filename, sep='\t', header = TRUE)
data[,-1] <- round(data[,-1], 3) # the "-1" excludes column 1
#write.xlsx(data, file.path(output_folder,'means.xlsx'))
write.table(data, file.path(output_folder,'means_small.csv'), 
            sep = '\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

print("DONE")
