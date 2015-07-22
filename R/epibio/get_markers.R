library(data.table)
BLOOD_METH_THRESHOLD <- 0.8
BLOOD_UNMETH_THRESHOLD <- 0.2
OTHER_TISSUES_METH_THRESHOLD <- 0.8
OTHER_TISSUES_UNMETH_THRESHOLD <- 0.2
TISSUE_UNMETH_THRESHOLD <- 0.5
TISSUE_METH_THRESHOLD <- 0.5
NUM_ALLOWED_BLOOD <- 0
nrows=-1

get_rows_methylated_in_blood <- function(blood.table.10,BLOOD_METH_THRESHOLD){
	num.blood <- ncol(blood.table.10)
	blood.methylated <- blood.table.10 > BLOOD_METH_THRESHOLD
	blood.methylated.count <- apply(blood.methylated,1,sum)
	rows.methylated.in.blood <- which(blood.methylated.count >= (num.blood-NUM_ALLOWED_BLOOD))
	return(rows.methylated.in.blood)
}

get_rows_unmethylated_in_blood <- function(blood.table.90,BLOOD_UNMETH_THRESHOLD){
	num.blood <- ncol(blood.table.90)
	blood.unmethylated <- blood.table.90 < BLOOD_UNMETH_THRESHOLD
	blood.unmethylated.count <- apply(blood.unmethylated,1,sum)
	rows.unmethylated.in.blood <- which(blood.unmethylated.count >= (num.blood-NUM_ALLOWED_BLOOD))
	return(rows.unmethylated.in.blood)
  }
get_unmeth_markers <- function(tissue,blood.table.10.meth_blood,tissues.table.10.meth_blood,tissues.table.90.meth_blood,TISSUE_UNMETH_THRESHOLD,OTHER_TISSUES_METH_THRESHOLD){
	other.tissues <- tissues.table.10.meth_blood[,setdiff(names(tissues.table.10.meth_blood),tissue)]
	num.tissues <- ncol(other.tissues)
	target.tissue <- tissues.table.90.meth_blood[,tissue,drop=FALSE]
	target.tissue.unmeth <- target.tissue < TISSUE_UNMETH_THRESHOLD
	
	other.tissues.meth <- tissues.table.10.meth_blood > OTHER_TISSUES_METH_THRESHOLD
	other.tissues.meth.count <- data.frame(apply(other.tissues.meth,1,sum,na.rm=T,drop=FALSE))
	other.tissues.meth.many <- other.tissues.meth.count > (num.tissues/2)

	markers <- data.frame(target.tissue,blood.table.10.meth_blood,tissues.table.10.meth_blood,other.tissues.meth.count)
	markers <- markers[target.tissue.unmeth & other.tissues.meth.many,]
	#markers <- markers[target.tissue.unmeth,]
	names(markers)[ncol(markers)] <- 'num_methylated_tissues'
	markers <- markers[order(markers[,ncol(markers)],decreasing=T),]
	return(markers)
}
get_meth_markers <- function(tissue,blood.table.90.unmeth_blood,tissues.table.10.unmeth_blood,tissues.table.90.unmeth_blood,TISSUE_METH_THRESHOLD,OTHER_TISSUES_UNMETH_THRESHOLD){
	other.tissues <- tissues.table.90.unmeth_blood[,setdiff(names(tissues.table.90.unmeth_blood),tissue)]
	num.tissues <- ncol(other.tissues)
	target.tissue <- tissues.table.10.unmeth_blood[,tissue,drop=FALSE]
	target.tissue.meth <- target.tissue > TISSUE_METH_THRESHOLD
	
	other.tissues.unmeth <- tissues.table.90.unmeth_blood < OTHER_TISSUES_UNMETH_THRESHOLD
	other.tissues.unmeth.count <- data.frame(apply(other.tissues.unmeth,1,sum,na.rm=T,drop=FALSE))
	other.tissues.unmeth.many <- other.tissues.unmeth.count > (num.tissues/2)
	markers <- data.frame(target.tissue,blood.table.90.unmeth_blood,tissues.table.90.unmeth_blood,other.tissues.unmeth.count)
	markers <- markers[target.tissue.meth & other.tissues.unmeth.many,]
	#markers <- markers[target.tissue.meth,]
	names(markers)[ncol(markers)] <- 'num_unmethylated_tissues'
	markers <- markers[order(markers[,ncol(markers)],decreasing=T),]
	return(markers)
}
geos.10 <- read.table('../../generated/summary/GEO/perc_10_healthy.txt',sep='\t',header=T,row.names=1,nrows=nrows)
tcga.10 <- read.table('../../generated/summary/TCGA/perc_10_healthy.txt',sep='\t',header=T,row.names=1,nrows=nrows)
geos.90 <- read.table('../../generated/summary/GEO/perc_90_healthy.txt',sep='\t',header=T,row.names=1,nrows=nrows)
tcga.90 <- read.table('../../generated/summary/TCGA/perc_90_healthy.txt',sep='\t',header=T,row.names=1,nrows=nrows)
lab_data.10 <-  read.table('../../generated/summary/lab_data/perc_10.txt',sep='\t',header=T,row.names=1,nrows=nrows)
lab_data.90 <-  read.table('../../generated/summary/lab_data/perc_90.txt',sep='\t',header=T,row.names=1,nrows=nrows)



geo.blood.cols <- c("healthy.cd14._monocytes","healthy.cd19._b.cells","healthy.cd19._cells",
"healthy.cd25.cd127._regulatory_t.cells","healthy.cd4._t.cells","healthy.cd45ro.ra._memory_t.cells",
"healthy.cd45ro.ra._naive_t.cells","healthy.cd56._nk_cells","healthy.cmp_cells_.bone_marrow.",
"healthy.eosinophils","healthy.gmp_cells_.bone_marrow.","healthy.immortalized_b.cells",
"healthy.lymphocyte","healthy.naive_b.cells","healthy.whole_blood","healthy_new_born.whole_blood",
"healthy_with_folate_supplementation.cd16.","healthy_with_folate_supplementation.whole_blood")
#geo.non.blood.cols <- setdiff(names(geos),geo.blood.cols)
lab_data.blood.cols <- 'leukocyte'
blood.cols <- c(geo.blood.cols,lab_data.blood.cols)
#convert to data.table
geos.10 <- data.table(geos.10,keep.rownames=T,key="rn")
tcga.10 <- data.table(tcga.10,keep.rownames=T,key="rn")
geos.90 <- data.table(geos.90,keep.rownames=T,key="rn")
tcga.90 <- data.table(tcga.90,keep.rownames=T,key="rn")
lab_data.10 <- data.table(lab_data.10,keep.rownames=T,key="rn")
lab_data.90 <- data.table(lab_data.90,keep.rownames=T,key="rn")

# merge and back to data frame
#all.table.10 <- data.frame(merge(geos.10,tcga.10),row.names=1)
#all.table.90 <- data.frame(merge(geos.90,tcga.90),row.names=1)
all.table.10 <- data.frame(merge(merge(geos.10,tcga.10),lab_data.10),row.names=1)
all.table.90 <- data.frame(merge(merge(geos.90,tcga.90),lab_data.90),row.names=1)
blood.cols <- intersect(blood.cols, colnames(all.table.10))

blood.table.10 <- all.table.10[,blood.cols]
blood.table.90 <- all.table.90[,blood.cols]

tissues.cols <- setdiff(names(all.table.10),blood.cols)
tissues.table.10 <- all.table.10[,tissues.cols]
tissues.table.90 <- all.table.90[,tissues.cols]

# get markers
dir.create('../../generated/markers',showWarnings=FALSE)
setwd('../../generated/markers')

rows.meth.blood <- get_rows_methylated_in_blood(blood.table.10,BLOOD_METH_THRESHOLD)
rows.unmeth.blood <- get_rows_unmethylated_in_blood(blood.table.10,BLOOD_UNMETH_THRESHOLD)
blood.table.10.meth_blood <- blood.table.10[rows.meth.blood,]
blood.table.90.unmeth_blood <- blood.table.90[rows.unmeth.blood,]

tissues.table.10.meth_blood <- tissues.table.10[rows.meth.blood,]
tissues.table.90.meth_blood <- tissues.table.90[rows.meth.blood,]
tissues.table.10.unmeth_blood <- tissues.table.10[rows.unmeth.blood,]
tissues.table.90.unmeth_blood <- tissues.table.90[rows.unmeth.blood,]

for (tissue in tissues.cols){
	unmeth.markers <- get_unmeth_markers(tissue,blood.table.10.meth_blood,tissues.table.10.meth_blood,tissues.table.90.meth_blood,TISSUE_UNMETH_THRESHOLD,OTHER_TISSUES_METH_THRESHOLD)
	unmeth.markers.file <- paste(tissue,'_unmeth.txt',sep='')
	write.table(unmeth.markers,unmeth.markers.file,sep='\t',col.names=NA)
	num.markers <- nrow(unmeth.markers)
	print(paste('Found',num.markers,'unmethylated markers for',tissue))

	meth.markers <- get_meth_markers(tissue,blood.table.90.unmeth_blood,tissues.table.10.unmeth_blood,tissues.table.90.unmeth_blood,TISSUE_METH_THRESHOLD,OTHER_TISSUES_UNMETH_THRESHOLD)
	meth.markers.file <- paste(tissue,'_meth.txt',sep='')
	write.table(meth.markers,meth.markers.file,sep='\t',col.names=NA)
	num.markers <- nrow(meth.markers)
	print(paste('Found',num.markers,'methylated markers for',tissue))
}




