
source("config.R")
source("common.R")

normalize_names <- function(name) {
  name <- tolower(name)
  # removes
  remove_re <- c('.*___', '.txt.gz')
  name <- mgsub(remove_re, character(length(remove_re)), name)
  # converts
  name <- gsub('healthy_obese', 'obese', name)
  name <- gsub('healthy_control', 'healthy', name)
  name <- gsub('healthy_(leiomyoma_control)', 'healthy', name)
  name
}

readOrganizedFile <- function(filename, folder) {
  fd <- gzfile(file.path(folder, filename))
  betas.table <- read.table(fd, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE)
  print(sprintf('-> Read file %s with %d samples', filename, ncol(betas.table)))
  betas.table
}

workOnKind <- function(group, folder, output_folder) {
  name <- group$normalized[[1]]
  print(sprintf('Working on %s (%d files)', name, length(group$filename)))
  betas.table <- do.call("cbind", lapply(group$filename, FUN=readOrganizedFile, folder))
  samples_num <- ncol(betas.table)
  mean <- rowMeans(betas.table)
  std <- apply(betas.table, 1, sd)
  merged_data <- data.frame(mean, std, samples_num)
  
  output_filename <- file.path(output_folder, paste0(name, '.csv'))
  write.table(merged_data, output_filename, sep='\t', col.names=NA, quote=FALSE)
}

dir.create(generated_merged_folder, recursive=TRUE, showWarnings=FALSE)
#beta_files_tcga <- list.files(generated_TCGA_folder, pattern="*.txt.gz")
beta_files <- list.files(generated_GEO_folder, pattern="*.txt.gz")
#beta_files <- beta_files[c(6,16,17, 88,214, 87,213,209)] # XXX
df <- data.frame(filename=beta_files)
df$normalized <- sapply(beta_files, FUN=normalize_names)
s <- split(df, df$normalized, drop=TRUE)
lapply(s, FUN=workOnKind, generated_GEO_folder, generated_merged_folder)

