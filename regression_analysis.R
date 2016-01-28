### INIT ####

# necessary libraries
library(GEOquery)

# path to data
DATADIR = "~/wrk/dippa-data/"

# Load data
mrna <- getGEO(filename = file.path(DATADIR,"GSE58212_series_matrix.txt.gz"))
mirna <- getGEO(filename = file.path(DATADIR,"GSE58210_series_matrix.txt.gz"))
prot <- read.delim(file.path(DATADIR,"Oslo2-RPPA_data.csv"), stringsAsFactors = F)

