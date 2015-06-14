
library(data.table)

trim <- function (x) {
  gsub("^\\s+|\\s+$", "", x)
}

get_indices_to_runon <- function(vec, args) {
  if(length(args) > 0) {
    indices <- c(as.numeric(args[1]))
  } else {
    indices <- seq_along(vec)  
  }
  indices
}

list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
                    full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names)) {
    dirs
  } else {
    basename(dirs)
  }
}

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
  result
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

create_name <- function(study, type) {
  name <- paste0(study, ".", type)
  rep_pat <- c(" ", "/", ":", ">", "<")
  name <- mgsub(rep_pat, rep(c("_"), length(rep_pat)), name , fixed=TRUE)
  name
}

get_output_filename <- function(folder, fn_prefix, study, type) {
  name <- create_name(study, type)  
  fn <- file.path(folder, paste0(fn_prefix, '___', name, '.txt.gz'))
  fn
}

write_beta_values_table <- function(output_filename, betas.table) {
  fd <- gzfile(output_filename)
  write.table(betas.table, fd, sep='\t', col.names=NA, quote=FALSE)
}

#' paste which suppress NAs
#' based on http://stackoverflow.com/questions/13673894/suppress-nas-in-paste
paste3 <- function(...,sep="") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}
