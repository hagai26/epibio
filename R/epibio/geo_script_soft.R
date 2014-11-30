
require(minfi)

path <- file.path("..\\..\\data\\GEO\\GSE52621_family.soft")

#g <- read.450k.GEO(path=path)

gset <- GEOquery::getGEO(filename = file.path(path, list.files(path, pattern = ".soft")))
