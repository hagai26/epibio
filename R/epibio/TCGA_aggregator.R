
require(minfi)

source("config.R")
source("common.R")


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
        full_filename <- file.path(generated_TCGA_folder, kind_output_filename)
        write.csv(cbind(CpG=rownames(level_data), level_data), file = full_filename, row.names=FALSE)
      }
    }
  }
}


dir.create(generated_TCGA_folder, recursive=TRUE, showWarnings=FALSE)

# targets
targets <- read.csv(samples_filename, stringsAsFactors = FALSE)
targets$Basename <- file.path(idat_folder, paste(targets$Sentrix.ID, targets$Sentrix.position, sep="_"))

result <- chunked_group_by(targets, list(targets$Study, targets$Type), 2)

all_kinds = data.frame(number = sapply(result$splited, FUN=nrow))
print("Reading all kinds:")
print(all_kinds)
print("")
ret <- lapply(result$grouped, FUN=work_on_targets, NULL)

all_kinds_filename <- file.path(generated_TCGA_folder, 'TCGA_all_kinds.csv')
write.csv(cbind(kind=rownames(all_kinds), all_kinds), file = all_kinds_filename, row.names=FALSE, quote=FALSE)
print("DONE")
