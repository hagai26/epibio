
library(data.table)
library(openxlsx)

source("config.R")
source("common.R")

input_folder <- file.path(summary_folder, 'GEO')
csv_filename <- file.path(input_folder, 'means.txt')
data <- read.table(csv_filename, sep='\t', header = TRUE)
print(head(data))
#write.xlsx(data, file.path(output_folder,'means.xlsx'))

print("DONE")
