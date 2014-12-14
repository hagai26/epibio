

chunked_group_by <- function(targets, group_by_list, chunk_size) {
  splited_targets <- split(targets, group_by_list, drop=TRUE)
  splited_targets <- lapply(splited_targets, FUN=function(x) head(x, 20))  #XXX
  
  chunked_targets <- split(splited_targets, ceiling(seq_along(splited_targets)/chunk_size))
  grouped_chunked_targets <- lapply(chunked_targets, function(x) do.call("rbind", x))
  
  result <- list("grouped"=grouped_chunked_targets, "splited"=splited_targets)
  return(result)
}
