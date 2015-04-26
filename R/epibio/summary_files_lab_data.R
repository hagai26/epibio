ll .gz files to 3 files: mean, 10th percentile, 90th percentile
library(data.table)

summary_files <- function(input_folder,output_folder){
	  mean.col <- 2
  tenth.col <- 4
    ninetieth.col <- 5
    gz.files <- list.files(input_folder)
	  new.col.names <- substr(basename(gz.files),1,nchar(gz.files)-nchar('.txt.gz'))
	  num.files <- length(gz.files)
	    #initiate variables using first file
	    data <- data.table(read.table(gz.files[1],sep='\t',header=T,row.names=1),keep.rownames=T,key='rn')
	    
	    setnames(data,c("mean","quantile_0.1","quantile_0.9"),rep(new.col.names[1],3))
		  means <- data[,c(1,mean.col),with=F]
		  perc.10 <- data[,c(1,tenth.col),with=F]
		    perc.90 <- data[,c(1,ninetieth.col),with=F]
		    
		    for (i in 2:num.files){
				    data <- data.table(read.table(gz.files[i],sep='\t',header=T,row.names=1),keep.rownames=T,key='rn')
			    setnames(data,c("mean","quantile_0.1","quantile_0.9"),rep(new.col.names[i],3))
				    means.new <- data[,c(1,mean.col),with=F]
				    perc.10.new <- data[,c(1,tenth.col),with=F]
					    perc.90.new <- data[,c(1,ninetieth.col),with=F]
					    means <- merge(means,means.new,all=T)
						    perc.10 <- merge(perc.10,perc.10.new,all=T)
						    perc.90 <- merge(perc.90,perc.90.new,all=T)
							  }
			  write.table(means,file.path(output_dir,'means.txt'),sep='\t',col.names=NA)
			  write.table(perc.10,file.path(output_dir,'perc_10.txt'),sep='\t',col.names=NA)
			    write.table(perc.90,file.path(output_dir,'perc_90.txt'),sep='\t',col.names=NA)
}

print('Summarizing GEO')
input_folder <- '/cs/icore/joshua.moss/dor/hagaic/epibio/generated/merged/lab_data')
output_folder <- '/cs/icore/joshua.moss/dor/hagaic/epibio/generated/summary/lab_data')
summary_files(input_folder,output_folder)

