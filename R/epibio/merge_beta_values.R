
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
  name <- gsub('temporal_lobe', 'temporal_cortex', name)
  name <- gsub('peripheral_blood_mononuclear_cells', 'pbmc', name)
  name
}

readOrganizedFile <- function(filename, folder) {
  print(sprintf('-> Read file: %s', filename))
  fd <- gzfile(file.path(folder, filename))
  betas.table <- read.table(fd, sep='\t', header=TRUE, row.names=1, 
                            stringsAsFactors=FALSE, comment.char="")
  print(sprintf('-> has %d samples', ncol(betas.table)))
  betas.table
}

workOnKind <- function(group, folder, output_folder) {
  name <- group$normalized[[1]]
  output_filename <- file.path(output_folder, paste0(name, '.txt.gz'))
  print(sprintf('Working on %s (%d files)', name, length(group$filename)))
  if (file.exists(output_filename)) {
    print(sprintf('file exists - skipping'))
  } else {
    betas.table <- do.call("cbind", lapply(group$filename, FUN=readOrganizedFile, folder))
    row.has.na <- apply(betas.table, 1, function(x) any(is.na(x)) )
    betas.table <- betas.table[!row.has.na,]
    samples_num <- ncol(betas.table)
    mean <- rowMeans(betas.table)
    std <- apply(betas.table, 1, sd)
    quantile_0.1 <- apply(betas.table, 1, quantile, probs=c(0.1))
    quantile_0.9 <- apply(betas.table, 1, quantile, probs=c(0.9))
    merged_data <- data.frame(mean, std, quantile_0.1, quantile_0.9, samples_num)
    fd <- gzfile(output_filename)
    write.table(merged_data, fd, sep='\t', col.names=NA, quote=FALSE)
  }
}

dir.create(generated_merged_folder, recursive=TRUE, showWarnings=FALSE)
#beta_files_tcga <- list.files(generated_TCGA_folder, pattern="*.txt.gz")
beta_files <- list.files(generated_GEO_folder, pattern="*.txt.gz")
#beta_files <- beta_files[c(6,16,17, 88,214, 87,213,209)] # XXX
df <- data.frame(filename=beta_files)
df$normalized <- sapply(beta_files, FUN=normalize_names)
kinds <- split(df, df$normalized, drop=TRUE)
c <- 1
for (kind in kinds) {
  print(sprintf("%d/%d", c, length(kinds)))
  workOnKind(kind, generated_GEO_folder, generated_merged_folder)
  c <- c+1
}

