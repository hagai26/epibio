
library(data.table)
library(openxlsx)

source("config.R")
source("common.R")

input_folder <- file.path(generated_merged_folder, 'GEO')
data <- read.table(csv_filename, sep='\t', header = TRUE)
print(head(data))
#write.xlsx(data, file.path(output_folder,'means.xlsx'))

print("DONE")
