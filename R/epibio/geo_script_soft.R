
require(minfi)

source("config.R")


path <- file.path(data_folder_big, "GEO\\GSE52621_family.soft")

g <- read.450k.GEO(path=path)

gset <- GEOquery::getGEO(filename = file.path(path, list.files(path, pattern = ".soft")))
