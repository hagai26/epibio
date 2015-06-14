## summarize .gz files to 3 files: mean, 10th percentile, 90th percentile
library(data.table)

summary_files <- function(input_folder,output_folder){
	  mean.col <- 2
  tenth.col <- 4
    ninetieth.col <- 5
    gz.files <- list.files(input_folder)
	  new.col.names <- substr(basename(gz.files),1,nchar(gz.files)-nchar('.txt.gz'))
	  num.files <- length(gz.files)
	    #initiate variables using first file
		print(paste('Summarizing ', new.col.names[1],' ',1,'/',num.files,sep=''))
		data.file <- file.path(input_folder,gz.files[1])
	    data <- data.table(read.table(data.file,sep='\t',header=T,row.names=1),keep.rownames=T,key='rn')
	    
	    setnames(data,c("mean","quantile_0.1","quantile_0.9"),rep(new.col.names[1],3))
		  means <- data[,c(1,mean.col),with=F]
		  perc.10 <- data[,c(1,tenth.col),with=F]
		    perc.90 <- data[,c(1,ninetieth.col),with=F]
		    
		    for (i in 2:num.files){
				print(paste('Summarizing ', new.col.names[i],' ',i,'/',num.files,sep=''))
				data.file <- file.path(input_folder,gz.files[i])
				data <- data.table(read.table(data.file,sep='\t',header=T,row.names=1),keep.rownames=T,key='rn')
			    setnames(data,c("mean","quantile_0.1","quantile_0.9"),rep(new.col.names[i],3))
				    means.new <- data[,c(1,mean.col),with=F]
				    perc.10.new <- data[,c(1,tenth.col),with=F]
					    perc.90.new <- data[,c(1,ninetieth.col),with=F]
					    means <- merge(means,means.new,all=T)
						    perc.10 <- merge(perc.10,perc.10.new,all=T)
						    perc.90 <- merge(perc.90,perc.90.new,all=T)
							  }
			  write.table(means,file.path(output_folder,'means.txt'),sep='\t',col.names=T,row.names=F)
			  write.table(perc.10,file.path(output_folder,'perc_10.txt'),sep='\t',col.names=T,row.names=F)
			    write.table(perc.90,file.path(output_folder,'perc_90.txt'),sep='\t',col.names=T,row.names=F)
}

dir.create('/cs/icore/joshua.moss/dor/hagaic/epibio/generated/summary',showWarnings=F)
print('Summarizing lab_data')
input_folder <- '/cs/icore/joshua.moss/dor/hagaic/epibio/generated/merged/lab_data'
output_folder <- '/cs/icore/joshua.moss/dor/hagaic/epibio/generated/summary/lab_data'
dir.create(output_folder,showWarnings=F)
summary_files(input_folder,output_folder)

