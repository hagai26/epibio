
library(GEOquery)

source("config.R")

dir.create(generated_GEO_folder, recursive=TRUE, showWarnings=FALSE)

geo_path <- file.path(data_folder_big, "GEO")
geo_folders <- file.path(geo_path, list.files(geo_path))
soft_filenames <- file.path(geo_folders, list.files(geo_folders, pattern = ".soft$"))
soft_filenames <- soft_filenames[2] # XXX

soft_filename <-  soft_filenames[[1]]

# XXX  - GSElimits = c(1, 10)
gse <- GEOquery::getGEO(filename = soft_filename, GSElimits = c(1, 10))
characteristics <- lapply(GSMList(gse), FUN=function(x) x@header$characteristics_ch1)
gender <- unlist(lapply(characteristics, FUN=function(x) x[[1]]), use.names = FALSE)
type <- unlist(lapply(characteristics, FUN=function(x) x[[2]]), use.names = FALSE)
t <- data.frame(name=paste0(gse@header$geo_accession, ".", names(GSMList(gse))), gender=gender, type=type)

splited_targets <- split(t, list(t$type, t$gender), drop=TRUE)
#gse@header$geo_accession
