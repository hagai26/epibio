
data_folder <- file.path("..", "..", "data")

generated_folder  <- file.path("..", "..", "generated")
generated_TCGA_folder <- file.path(generated_folder, "TCGA")
generated_GEO_folder <- file.path(generated_folder, "GEO")

idat_folder <- file.path(data_folder, "small", "TCGA_L1", "450K_idat")
samples_filename <- file.path(idat_folder, "samples.csv")

big_data_folder <- file.path(data_folder, "big")
small_data_folder <- file.path(data_folder, "small")

external_disk_data_path <- file.path(data_folder, 'external_disk')
