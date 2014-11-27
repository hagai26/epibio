
library(RnBeads)

data_folder <- "..\\..\\data"
generated_folder  <- file.path(data_folder, "generated")
TCGA_folder <- file.path(generated_folder, "TCGA")
idat.dir <- file.path(data_folder, "static\\L3_data_small\\idat")
dir.create(TCGA_folder, recursive=TRUE, showWarnings=FALSE)

sample.annotation <- file.path(idat.dir, "samples.csv")

# targets
targets <- read.csv(file.path(idat.dir, "samples.csv"), stringsAsFactors = FALSE)
targets$Basename <- file.path(idat.dir, paste(targets$Sentrix.ID, targets$Sentrix.position, sep="_"))

report.dir <- file.path("reports")

data.source <-c(idat.dir, sample.annotation)
result <- rnb.execute.import(data.source=data.source, data.type="infinium.idat.dir")

rnb.set <- result$rnb.set

head(meth(rnb.set, type="promoters"))
