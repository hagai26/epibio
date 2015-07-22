summary.dir <- '/cs/icore/joshua.moss/dor/hagaic/epibio/generated/summary'

nrow=-1
geo.dir <- file.path(summary.dir,'GEO')
tcga.dir <- file.path(summary.dir,'TCGA')
lab_data.dir <- file.path(summary.dir,'lab_data')
geo.perc_10.file <- file.path(geo.dir,'perc_10.txt')
geo.perc_10.data <- read.table(geo.perc_10.file,sep='\t', header=T, row.names=1,nrow=nrow)
geo.healthy.cols <- names(geo.perc_10.data)[substr(names(geo.perc_10.data),1,7)=='healthy']

geo.perc_10.data.healthy <- geo.perc_10.data[,geo.healthy.cols]

geo.perc_90.file <- file.path(geo.dir,'perc_90.txt')
geo.perc_90.data <- read.table(geo.perc_90.file,sep='\t', header=T, row.names=1,nrow=nrow)
geo.perc_90.data.healthy <- geo.perc_90.data[,geo.healthy.cols]

geo.means.file <- file.path(geo.dir,'means.txt')
geo.means.data <- read.table(geo.means.file,sep='\t', header=T, row.names=1,nrow=nrow)
geo.means.data.healthy <- geo.means.data[,geo.healthy.cols]

write.table(geo.perc_10.data.healthy,file.path(geo.dir,'perc_10_healthy.txt'),sep='\t',col.names=NA)
write.table(geo.perc_90.data.healthy,file.path(geo.dir,'perc_90_healthy.txt'),sep='\t',col.names=NA)
write.table(geo.means.data.healthy,file.path(geo.dir,'means_healthy.txt'),sep='\t',col.names=NA)
print('Done GEO!')


tcga.perc_10.file <- file.path(tcga.dir,'perc_10.txt')
tcga.perc_10.data <- read.table(tcga.perc_10.file,sep='\t', header=T, row.names=1,nrow=nrow)
#tcga.healthy.cols <- names(tcga.perc_10.data)[substr(names(tcga.perc_10.data),1,19)=='solid_tissue_normal']
tcga.healthy.cols <- names(tcga.perc_10.data)[grepl('*solid_tissue_normal*',names(tcga.perc_10.data))]
#tcga.healthy.cols <- setdiff(tcga.healthy.cols,c('solid_tissue_normal.other__specify','solid_tissue_normal..not_available.','solid_tissue_normal.mixed_histology_.please_specify.'))

tcga.perc_10.data.healthy <- tcga.perc_10.data[,tcga.healthy.cols]

tcga.perc_90.file <- file.path(tcga.dir,'perc_90.txt')
tcga.perc_90.data <- read.table(tcga.perc_90.file,sep='\t', header=T, row.names=1,nrow=nrow)
tcga.perc_90.data.healthy <- tcga.perc_90.data[,tcga.healthy.cols]

tcga.means.file <- file.path(tcga.dir,'means.txt')
tcga.means.data <- read.table(tcga.means.file,sep='\t', header=T, row.names=1,nrow=nrow)
tcga.means.data.healthy <- tcga.means.data[,tcga.healthy.cols]

#write results
write.table(tcga.perc_10.data.healthy,file.path(tcga.dir,'perc_10_healthy.txt'),sep='\t',col.names=NA)
write.table(tcga.perc_90.data.healthy,file.path(tcga.dir,'perc_90_healthy.txt'),sep='\t',col.names=NA)
write.table(tcga.means.data.healthy,file.path(tcga.dir,'means_healthy.txt'),sep='\t',col.names=NA) 
print('Done TCGA!') 
