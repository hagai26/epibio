
require(minfi)


work_on_targets <- function(targets, baseDir) {
  study_levels <- levels(factor(targets$Study))
  type_levels <- levels(factor(targets$Type))
  
  ptime1 <- proc.time()
  cat("Reading", nrow(targets), "pairs of idat files...")
  RGSet <- read.450k.exp(base = baseDir, targets = targets)
  MSet <- preprocessIllumina(RGSet, bg.correct = TRUE, normalize = "controls")
  ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
  beta <- head(getBeta(ratioSet), 2000) #XXX
  beta <- beta[order(row.names(beta)),, drop = FALSE]
  ptime2 <- proc.time()
  stime <- (ptime2 - ptime1)[3]
  cat(" in", stime, "seconds\n")

  for(type in type_levels) {
    for(study in study_levels) {
      relevant_targets <- targets$Study==study & targets$Type==type
      if(sum(relevant_targets) > 0) {
        relevant_beta <- beta[,relevant_targets, drop = FALSE]
        mean <- rowMeans(relevant_beta)
        std <- apply(relevant_beta, 1, sd)
        level_data <- data.frame(mean, std)
        
        kind_name <- paste0(study, ".", type)
        kind_output_filename <- paste0(gsub(" ","_", kind_name , fixed=TRUE), ".csv")
        full_filename <- file.path(new_TCGA_folder, kind_output_filename)
        write.csv(cbind(CpG=rownames(level_data), level_data), file = full_filename, row.names=FALSE)
      }
    }
  }
}

# folders
data_folder <- file.path("..", "..", "data_small")
generated_folder  <- file.path("..", "..", "generated")
new_TCGA_folder <- file.path(generated_folder, "TCGA")
idat_folder <- file.path(data_folder, "TCGA_L1", "450K_idat")
dir.create(new_TCGA_folder, recursive=TRUE, showWarnings=FALSE)

# targets
targets <- read.csv(file.path(idat_folder, "samples.csv"), stringsAsFactors = FALSE)
targets$Basename <- file.path(idat_folder, paste(targets$Sentrix.ID, targets$Sentrix.position, sep="_"))
splited_targets <- split(targets,list(targets$Study, targets$Type), drop=TRUE)
splited_targets <- lapply(splited_targets, FUN=function(x) head(x, 20))  #XXX

all_kinds = data.frame(number = sapply(splited_targets, FUN=nrow))
print("Reading all kinds:")
print(all_kinds)
print("")
chunk_size <- 3
chunked_targets <- split(splited_targets, ceiling(seq_along(splited_targets)/chunk_size))
grouped_chunked_targets <- lapply(chunked_targets, function(x) do.call("rbind", x))
ret <- lapply(grouped_chunked_targets, FUN=work_on_targets, NULL)

all_kinds_filename = file.path(generated_folder, 'TCGA_all_kinds.csv')
write.csv(cbind(kind=rownames(all_kinds), all_kinds), file = all_kinds_filename, row.names=FALSE)
print("DONE")
