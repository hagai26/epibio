
library(RnBeads)


workOnIdatsFolder <- function(idat_folder, targets, identifiers_column) {
  data.source <- list(idat_folder, targets)
  rnb.options(identifiers.column = identifiers_column)
  rnb.set <- rnb.execute.import(data.source=data.source, data.type="infinium.idat.dir")
  betas.table <- process_rnb_set_to_betas(rnb.set, TRUE)
  betas.table
}
# 

process_rnb_set_to_betas <- function(rnb.set, has_pvalues) {
  #print('process_rnb_set_to_betas called')
  
  #tryCatch({
   # rnb.set <- rnb.execute.snp.removal(rnb.set)$dataset
  #}, error = function(err) {
    # TODO
    
    # on GSE36278 this causes stops with error:
    # Error in checkSlotAssignment(object, name, value) : 
    # assignment of an object of class “numeric” is not valid for slot ‘M’ in an object of class “RnBeadRawSet”; is(value, "matrixOrffOrNULL") is not TRUE
    
    # on GSE42118 (and GSE52576, GSE44667, GSE53740):
    # <simpleError in checkSlotAssignment(object, name, value): 
    # assignment of an object of class “integer” is not valid for slot ‘M’ in an object of class “RnBeadRawSet”; is(value, "matrixOrffOrNULL") is not TRUE>
    
    # currently, we ignore this error and doesn't call snp.removal
 #   print(err)
 # })
 
  # XXX
  #rnb.set <- rnb.execute.normalization(rnb.set, method="bmiq",bgcorr.method="methylumi.lumi")
  betas.table <- meth(rnb.set, row.names=TRUE)
  if(has_pvalues) {
    pvalue.high <- which(dpval(rnb.set) > 0.05, arr.ind=TRUE)
    betas.table[pvalue.high[,'row'], pvalue.high[,'col']] <- NA
  }
  destroy(rnb.set)  
  betas.table
}
