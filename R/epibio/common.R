
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
