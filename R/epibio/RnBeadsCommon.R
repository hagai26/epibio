
library(RnBeads)

workOnIdatsFolder <- function(idat_folder, targets, output_filename) {
  data.source <-list(idat_folder, targets)
  rnb.options(identifiers.column = 'barcode')
  rnb.set <- rnb.execute.import(data.source=data.source, data.type="infinium.idat.dir")
  betas.table <- process_rnb_set_to_betas(rnb.set, TRUE)
  write_beta_values_table(output_filename, betas.table)
}

