library(RnBeads)
source('common.R')

data_folder <- '/cs/icore/joshua.moss/dor/lab_data'
idat_dir <- file.path(data_folder,'idats')
pheno_file <- file.path(data_folder,'samples_cluster.csv')
output_dir <- '/cs/icore/joshua.moss/dor/hagaic/epibio/generated/lab_data'
#output_file <- file.path(output_dir,'lab_betas.txt')
rnb.options(identifiers.column='Sample_Name')
rnb.set <- rnb.execute.import(list(idat_dir,pheno_file))
betas.table <- process_rnb_set_to_betas(rnb.set,TRUE)

pheno <- pheno(rnb.set)
sample_groups <- unique(pheno$Sample_Group)

for (i in 1:length(sample_groups)){
	cols <- which(pheno$Sample_Group == sample_groups[i])
	output_file <- gzfile(file.path(output_dir,paste(sample_groups[i],'.txt.gz',sep='')))
	write.table(betas.table[,cols],output_file,sep='\t',col.names=NA)
	}

