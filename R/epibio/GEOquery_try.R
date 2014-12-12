
library(GEOquery)

source("config.R")

geo_path <- file.path(data_folder_big, "GEO")
geo_folders <- list.files(geo_path)
soft_filenames <- file.path(geo_folders, list.files(path, pattern = ".soft.gz$"))
soft_filenames <- head(soft_filenames, 1) # XXX
# when reading from soft.gz it can make these warnings: "seek on a gzfile connection returned an internal error"
# which can be ignored
# https://stat.ethz.ch/pipermail/bioconductor/2011-September/040957.html
gse <- GEOquery::getGEO(filename = soft_filename)

# GPL8490   is Illumina HumanMethylation27 BeadChip
# GPL13534  is Illumina HumanMethylation450 BeadChip (HumanMethylation450_15017482)
print (names(GPLList(gsm)))
print(names(GSMList(gsm)))

print(GSMList(gsm)[[1]])
