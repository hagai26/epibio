
library(data.table)


chunked_group_by <- function(targets, group_by_list, chunk_size) {
  orig_colnames <- colnames(targets)
  targets$tmp_row_name <- rownames(targets)
  splited_targets <- split(targets, group_by_list, drop=TRUE)
  chunked_targets <- split(splited_targets, ceiling(seq_along(splited_targets)/chunk_size))
  grouped_chunked_targets <- lapply(chunked_targets, 
                                    function(x) as.data.frame(rbindlist(x)))
  t <- lapply(grouped_chunked_targets, function(x) as.data.frame(x, row.names=x$tmp_row_name))
  grouped <- lapply(t, function(x) subset(x, select=orig_colnames))
  result <- list("grouped"=grouped, "splited"=splited_targets)
  return(result)
}

write_nrow_per_group <- function(splited_targets, filepath) {
  all_kinds = data.frame(number = sapply(splited_targets, FUN=nrow))
  write.csv(cbind(kind=rownames(all_kinds), all_kinds), file = filepath, row.names=FALSE, quote=FALSE)
}

mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern) != length(replacement)) {
    stop("mgsub: pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}


process_rnb_set_to_betas <- function(rnb.set, has_pvalues) {
  logger.start(fname=NA)
  rnb.options(disk.dump.big.matrices=TRUE)
  rnb.options(enforce.memory.management=TRUE)
  
  tryCatch({
    rnb.set <- rnb.execute.snp.removal(rnb.set)$dataset
  }, error = function(err) {
    # TODO
    
    # on GSE36278 this causes stops with error:
    # Error in checkSlotAssignment(object, name, value) : 
    # assignment of an object of class “numeric” is not valid for slot ‘M’ in an object of class “RnBeadRawSet”; is(value, "matrixOrffOrNULL") is not TRUE
    
    # on GSE42118 (and GSE52576, GSE44667, GSE53740):
    # <simpleError in checkSlotAssignment(object, name, value): 
    # assignment of an object of class “integer” is not valid for slot ‘M’ in an object of class “RnBeadRawSet”; is(value, "matrixOrffOrNULL") is not TRUE>
    
    # currently, we ignore this error and doesn't call snp.removal
    print(err)
  })
  
  
  #rnb.set <- rnb.execute.normalization(rnb.set, 
  #                                     method="bmiq",bgcorr.method="methylumi.lumi")
  betas.table <- meth(rnb.set, row.names=TRUE)
  if(has_pvalues) {
    pvalue.high <- which(dpval(rnb.set) > 0.05, arr.ind=TRUE)
    betas.table[pvalue.high[,'row'], pvalue.high[,'col']] <- NA
  }
  destroy(rnb.set)  
  betas.table
}

create_name <- function(study, type) {
  name <- paste0(study, ".", type)  
  name
}

write_beta_values_table <- function(folder, fn_prefix, study, type, betas.table) {
  name <- create_name(study, type)  
  fn <- file.path(folder, paste0(fn_prefix, '__', mgsub(c(" ", "/"), rep(c("_"), 2), name , fixed=TRUE), '.txt'))
  write.table(betas.table, fn, sep='\t', col.names=NA, quote=FALSE)
}
