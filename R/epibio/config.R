
data_folder <- file.path("..", "..", "data", "small")
generated_folder  <- file.path("..", "..", "generated")
idat_folder <- file.path(data_folder, "TCGA_L1", "450K_idat")
samples_filename <- file.path(idat_folder, "samples.csv")

generated_TCGA_folder <- file.path(generated_folder, "TCGA")
generated_GEO_folder <- file.path(generated_folder, "GEO")

data_folder_big <- file.path("..", "..", "data", "big")
