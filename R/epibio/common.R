
library(reshape)

mymerge <- function(x) {
  if(length(x) == 1) {
    return(x[[1]])
  } else {
    return (merge_recurse(x))
  }
}

chunked_group_by <- function(targets, group_by_list, chunk_size) {
  splited_targets <- split(targets, group_by_list, drop=TRUE)
  chunked_targets <- split(splited_targets, ceiling(seq_along(splited_targets)/chunk_size))
  grouped_chunked_targets <- lapply(chunked_targets, function(x) mymerge(x))
  
  result <- list("grouped"=grouped_chunked_targets, "splited"=splited_targets)
  return(result)
}
