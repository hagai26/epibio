## Example differential methylation analysis using minfi and SWAN in R
## More info @ http://www.bioconductor.org/packages/2.10/bioc/vignettes/minfi/inst/doc/minfi.pdf
require(minfi)
require(minfiData)

## directory that contains the idat files
idatDir = "/path/to/idat_files" 
## read the sample sheet that describes the experiment
targets = read.450k.sheet("/path/to/sample_sheet",pattern="Sample_Sheet.csv") 
## read idat files specified in targets file
rgSet = read.450k.exp(base = idatDir, targets = targets)

## rename samples with meaningful identifiers from targets file
targets = targets[order(targets$Sample_Name),]
rgSet = rgSet[,order(sampleNames(rgSet))]
sampleNames(rgSet) = targets$Sample_Name

## calculate the detection p-values
pVals = detectionP(rgSet)

## plot mean detection p-values for all samples
barplot(apply(pValsDave,2,mean),las=2,cex.names=0.5)
abline(h=0.05,col="red")

## generate a quality control report for the raw data
## samples with large mean detecton p-values will generally also look poor in one or more
## of the quality control plots and may need to be excluded from further analysis
qcReport(rgSet,sampNames=targets$Sample_Name,sampGroups=targets$Sample_Group,pdf="qcReport.pdf")

## SWAN can be used in 3 different ways
## 1) process the raw data to create a methylation set and normalize using SWAN
mSetSw = preprocessSWAN(rgSet)
## 2) the above command is equivalent to these two commands:
mSet = preprocessRaw(rgSet)
mSetSw = preprocessSWAN(rgSet, mSet)
## 3) you can first use Illumina's processing method for background correction and color normalization:
mSet = preprocessIllumina(rgSet, bg.correct = TRUE, normalize = "controls")
mSetSw = preprocessSWAN(rgSet, mSet)

## MDS plots are a quick way to get a sense of the relationship between samples 
mdsPlot(mSetSw,numPositions=1000,sampGroups=targets$Sample_Groups)

## filter out probes that have failed in one or more samples based on detection p-value
mSetSwFlt = mSetSw[rowSums(pVals <= 0.01) == ncol(pVals),]

## perform differential methylation analysis
dmp = dmpFinder(mSetSwFlt,pheno=as.vector(targets$Sample_Group, type="categorical", shrinkVar=TRUE)
