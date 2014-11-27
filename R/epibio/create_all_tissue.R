
library(tools)


myReadSheet <- function(filename) {
  data <- read.csv(file.path(TCGA_folder, filename), stringsAsFactors = FALSE)
  filename_no_ext = file_path_sans_ext(filename)
  colnames(data) <- c(names(data)[1], paste0(names(data)[-1], ".", filename_no_ext))
  data
}

# folders
data_folder <- "..\\..\\data"
generated_folder  <- file.path(data_folder, "generated")
all_folder <- file.path(generated_folder, "all")
dir.create(all_folder, recursive=TRUE, showWarnings=FALSE)
TCGA_folder <- file.path(generated_folder, "TCGA")

csv_files = list.files(TCGA_folder)
all_data_list <- lapply(csv_files, myReadSheet)
merged <- Reduce(function(...) merge(..., all=T), all_data_list)
even_indexes <- seq(2, length(merged), 2)
odd_indexes <- seq(3, length(merged), 2)
merged_mean = merged[c(1, even_indexes)]
merged_std = merged[c(1, odd_indexes)]

write.csv(merged_mean, file = file.path(all_folder, "all_tissues_mean.csv"), row.names=FALSE, quote=FALSE)
write.csv(merged_std, file = file.path(all_folder, "all_tissues_std.csv"), row.names=FALSE, quote=FALSE)
