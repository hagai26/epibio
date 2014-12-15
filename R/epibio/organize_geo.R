
library(RnBeads)

source("config.R")
source("common.R")
source("geo_l1_reader.R")


work_on_targets <- function(targets) {
  print("work_on_targets called")
  study_levels <- levels(factor(targets$disease))
  type_levels <- levels(factor(targets$tissue))
  
  ptime1 <- proc.time()
  cat("Reading", nrow(targets), "samples")
  print("")
  # work on these targets
  for(type in type_levels) {
    for(study in study_levels) {
      is_relevant_targets <- targets$disease==study & targets$tissue==type
      if(sum(is_relevant_targets) > 0) {
        relevant_targets <- targets[is_relevant_targets,, drop = FALSE]
        cat('currently reading:', type, study, "(", sum(is_relevant_targets), "samples) :")
        
        series_id <- levels(factor(targets$series_id))
        series_id <- sub(",.*", "", series_id) # read each sample only once
        print(series_id)
        
        for(one_series_id in series_id) {
          read_geo_l1_data(one_series_id)
        }
      }
      
    }
  }
}


dir.create(generated_GEO_folder, recursive=TRUE, showWarnings=FALSE)
folder <- file.path(data_folder, "global/GEO/joined")
joined_files <- file.path(folder, list.files(folder, pattern="*.txt"))

f1 <- "../../data/global/GEO/joined/GSE32079.txt"
f2 <- "../../data/global/GEO/joined/GSE32148.txt"
f3 <- "../../data/global/GEO/joined/GSE62992.txt"
f4 <- "../../data/global/GEO/joined/GSE62640.txt"
f5 <- "../../data/global/GEO/joined/GSE61653.txt"
f6 <- "../../data/global/GEO/joined/GSE61446.txt"
f7 <- "../../data/global/GEO/joined/GSE57894.txt"
f8 <- "../../data/global/GEO/joined/GSE57767.txt"
f9 <- "../../data/global/GEO/joined/GSE32146.txt"
f10 <- "../../data/global/GEO/joined/GSE30870.txt"
f11 <- "../../data/global/GEO/joined/GSE29290.txt"
joined_files <- c(f1, f2, f3, f4, f5, f6, f7, f8, f10, f11)
joined_files <- head(joined_files, 1)
series.info <- do.call("rbind", lapply(joined_files, function(fn) 
  data.frame(
    Filename=fn, 
    read.table(fn, sep='\t', row.names=1, header=TRUE, fill=TRUE, na.strings=c("NA", "0"), quote="\"", stringsAsFactors=FALSE)
    )
  ))
# get only relevant samples
relevant.samples.idx <- as.numeric(series.info$relevant)
relevant.samples.idx[is.na(relevant.samples.idx)] <- 0
relevant.samples.idx <- relevant.samples.idx == 1

pheno <- series.info[relevant.samples.idx, c('series_id', 'description','tissue','cell_type','disease')]
num_samples <- length(series.info[,1])

result <- chunked_group_by(pheno, list(pheno$disease, pheno$tissue), 3)

all_kinds = data.frame(number = sapply(result$splited, FUN=nrow))
all_kinds_filename <- file.path(generated_GEO_folder, 'GEO_all_kinds.csv')
write.csv(cbind(kind=rownames(all_kinds), all_kinds), file = all_kinds_filename, row.names=FALSE, quote=FALSE)

ret <- lapply(result$grouped, FUN=work_on_targets)
print("DONE")
