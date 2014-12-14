
# Import from GSE52621
meth.file1 <- "/Users/joshuamoss/Desktop/breast cancer serum analysis/GSE52621/GSE52621_meth.txt"
meth.table1 <- read.table(meth.file1,header = TRUE,row.names = 1,sep = "\t")
unmeth.file1 <- "/Users/joshuamoss/Desktop/breast cancer serum analysis/GSE52621/GSE52621_unmeth.txt"
unmeth.table1 <- read.table(unmeth.file1,header = TRUE,row.names = 1,sep = "\t")
pval.file1 <- "/Users/joshuamoss/Desktop/breast cancer serum analysis/GSE52621/GSE52621_pval.txt"
pval.table1 <- read.table(pval.file1,header = TRUE,row.names = 1,sep = "\t")
pheno.file1 <- "/Users/joshuamoss/Desktop/breast cancer serum analysis/GSE52621/GSE52621_sample.txt"
pheno.table1 <- read.table(pheno.file1,header = TRUE,row.names = 1,sep = "\t")

M = as.matrix(meth.table1)
U = as.matrix(unmeth.table1)
P = as.matrix(pval.table1)
GSE52621.rnb.set <- new("RnBeadRawSet",pheno = pheno.table1, M = M, U = U,platform = "27k", useff = FALSE, p.values = P)

GSE52621.rnb.set.mls <- as(GSE52621.rnb.set,"MethyLumiSet")