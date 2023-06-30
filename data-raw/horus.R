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

# Convert time (Francesco)
# pad 0 in front of single numbers
month.string <- as.character(horus_orig$Mo)
day.string <- as.character(horus_orig$Da)
hour.string <- as.character(horus_orig$Ho)
minute.string <- as.character(horus_orig$Mi)
second.string <- as.character(horus_orig$Se)

month.string[horus_orig$Mo < 10] <- paste0('0', month.string[horus_orig$Mo < 10])
day.string[horus_orig$Da < 10] <- paste0('0', day.string[horus_orig$Da < 10])
hour.string[horus_orig$Ho < 10] <- paste0('0', hour.string[horus_orig$Ho < 10])
minute.string[horus_orig$Mi < 10] <- paste0('0', minute.string[horus_orig$Mi < 10])
second.string[horus_orig$Se < 10] <- paste0('0', second.string[horus_orig$Se < 10])

horus_orig$time_string <- paste0(horus_orig$Year,'-',
                                 month.string,'-',
                                 day.string,'T',
                                 hour.string,':',
                                 minute.string,':',
                                 second.string)

# Set as active data set in the package ####

# usethis::use_data(horus, overwrite = TRUE)
