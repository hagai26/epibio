
source("config.R")
source("common.R")
args <- commandArgs(trailingOnly = TRUE)

normalize_names <- function(name) {
  name <- tolower(name)
  # removes
  remove_re <- c('gse.*___', '.txt.gz')
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
  #nrows = 10 # XXX (should be -1 on production)
  nrows = -1
  fd <- gzfile(file.path(folder, filename))
  betas.table <- read.table(fd, sep='\t', header=TRUE, row.names=1, 
                            stringsAsFactors=FALSE, comment.char="", 
                            nrow=nrows)
  print(sprintf('-> has %d samples', ncol(betas.table)))
  betas.table
}

workOnKind <- function(group, folder, output_folder) {
  name <- group$normalized[[1]]
  output_filename <- file.path(output_folder, paste0(name, '.txt.gz'))
  print(sprintf('Working on %s (%d files)', name, length(group$filename)))
  todo <- c()
  if(name %in% todo) {
    print("currently skipping this - TODO - fix BUGS on these")
  } else {
    if (file.exists(output_filename)) {
      print(sprintf('file exists - skipping'))
    } else {
      betas.table_list <- lapply(group$filename, FUN=readOrganizedFile, folder)
      # rename all columns to be with sample index in list (so we wont have same column names)
      for (i in 1:length(betas.table_list)) {
        colnames(betas.table_list[[i]]) <- paste0(colnames(betas.table_list[[i]]), '__', i)
      }
      # merge this list by row names
      betas.table_list_with_rn_col <- lapply(betas.table_list, 
                                         function(x) data.frame(x, Row.names = row.names(x)))
      merged_betas.table_with_rn_col <- Reduce(merge, betas.table_list_with_rn_col)
      betas.table <- transform(merged_betas.table_with_rn_col, row.names=Row.names, Row.names=NULL)
      
      row.has.na <- apply(betas.table, 1, function(x) any(is.na(x)) )
      betas.table <- betas.table[!row.has.na,, drop=FALSE]
      if(nrow(betas.table) == 0) {
        stop("nrow(betas.table) == 0")
      }
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
}

merge_beta_values <- function(generated_folder, output_folder) {
  dir.create(output_folder, recursive=TRUE, showWarnings=FALSE)
  beta_files <- list.files(generated_folder, pattern="*.txt.gz")
  if(length(beta_files) > 0) {
    df <- data.frame(filename=beta_files)
    df$normalized <- sapply(beta_files, FUN=normalize_names)
    groups <- split(df, df$normalized, drop=TRUE)
    indices <- get_indices_to_runon(groups, args)
    for (i in indices) {
      # prevent errors in middle of run
      # we could give indices which are good for geo or tcga and the other will skip them
      if(i <= length(groups)) {
  	    group <- groups[[i]]
  	    print(sprintf("%d/%d", i, length(groups)))
  	    workOnKind(group, generated_folder, output_folder)
      } else {
  	    print(sprintf("skipping %d/%d, index out of range", i, length(groups)))
  	  }
    }
  }
}

dir.create(generated_merged_folder, recursive=TRUE, showWarnings=FALSE)

print("merging GEO")
output_folder <- file.path(generated_merged_folder, 'GEO')
merge_beta_values(generated_GEO_folder, output_folder)

print("merging TCGA")
output_folder <- file.path(generated_merged_folder, 'TCGA')
merge_beta_values(generated_TCGA_folder, output_folder)

print("merging lab_data")
output_folder <- file.path(generated_merged_folder, 'lab_data')
merge_beta_values(generated_lab_data_folder, output_folder)

print("DONE")
