
library(RnBeads)

source("config.R")

dir.create(new_TCGA_folder, recursive=TRUE, showWarnings=FALSE)

# targets
targets <- read.csv(samples_filename, stringsAsFactors = FALSE)
targets$Basename <- file.path(idat_folder, paste(targets$Sentrix.ID, targets$Sentrix.position, sep="_"))

report.dir <- file.path("reports")

data.source <-c(idat_folder, samples_filename)
result <- rnb.execute.import(data.source=data.source, data.type="infinium.idat.dir")

rnb.set <- result$rnb.set

head(meth(rnb.set, type="promoters"))
