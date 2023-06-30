# Download and import original HORUS data ####

data_dir <- tempfile("HORUS_Ita_Catalog")
dir.create(data_dir)
data_file <- file.path(data_dir, "HORUS_Ita_Catalog.zip")

download.file("http://horus.bo.ingv.it/DataFolder/HORUS_Ita_Catalog.zip",
              destfile = file.path(data_dir, "HORUS_Ita_Catalog.zip"))

system(paste0("unzip -d ", data_dir, " ", data_file))

# load ISIDE catalogue
horus_orig <- read.table(
  file = file.path(data_dir, "HORUS_Ita_Catalog.txt"),
  header = TRUE,
  sep = "\t"
)

unlink(data_dir)

# Fix missing column name in downloaded data file:
colnames(horus_orig)[ncol(horus_orig)] <- "Iside.n."

# Convert to ETAS.inlabru subset format ####

# ??? TODO: fill in this code

# Set as active data set in the package ####

# usethis::use_data(horus, overwrite = TRUE)
