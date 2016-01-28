### INIT ####

# necessary libraries
library(GEOquery)

# path to data
DATADIR <- "~/wrk/dippa-data/"

# Constants for analysis
SIGMA <- 0.1  # used for guaranteeing positivity in normalization

# Load data
eset_mrna <- getGEO(filename = file.path(DATADIR,"GSE58212_series_matrix.txt.gz"))
eset_mirna <- getGEO(filename = file.path(DATADIR,"GSE58210_series_matrix.txt.gz"))
df_prot <- read.delim(file.path(DATADIR,"Oslo2-RPPA_data.csv"), stringsAsFactors = F)

# Convert Expression sets to matrices
mrna <- as.data.frame(exprs(eset_mrna))
mirna <- as.data.frame(exprs(eset_mirna))
# Convert from GSM IDs to OSLO2 patient IDs (with dash to dot for col names)
smplnames <- pData(eset_mirna)[colnames(mirna),"title"]
colnames(mirna) <- sub("-", ".", sub("BreastTumor_", "", smplnames))
smplnames <- pData(eset_mrna)[colnames(mrna),"title"]
colnames(mrna) <- sub("-", ".", sub("BreastTumor_", "", smplnames))
rm(smplnames)


### PREPROCESS ####

# Reformat protein data
rownames(df_prot) <- df_prot$Gene
# Expand AKT1/2/3-genes in protein data
i <- match("AKT1/AKT2/AKT3", df_prot$Gene)
prot <- df_prot[c(1:i, i, i, (i+1):nrow(df_prot)), 3:ncol(df_prot)]
rownames(prot)[i:(i+2)] <- c("AKT1","AKT2","AKT3")

# Normalization for regression (row-wise)
#norm.func <- function(row) { ()/() + SIGMA }

