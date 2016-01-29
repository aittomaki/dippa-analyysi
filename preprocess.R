### INIT ####

# necessary libraries
library(GEOquery)
library(componentSkeleton) #Anduril library

# path to data
DATADIR <- "~/wrk/dippa-data/"

# Constants for analysis
SIGMA <- 0.1  # used for guaranteeing positivity in normalization

# Load data
eset_mrna <- getGEO(filename = file.path(DATADIR,"GSE58212_series_matrix.txt.gz"))
eset_mirna <- getGEO(filename = file.path(DATADIR,"GSE58210_series_matrix.txt.gz"))
df_prot <- CSV.read(file.path(DATADIR,"Oslo2-RPPA_data.csv"))
annot_mrna <- getGEO(filename = file.path(DATADIR, "platforms","GPL14550.soft"))




### PREPROCESS ####

# Convert ExpressionSets to matrices
mrna <- as.matrix(exprs(eset_mrna))
mirna <- as.matrix(exprs(eset_mirna))

# Convert from GSM IDs to OSLO2 patient IDs (with dash to dot for col names)
smplnames <- pData(eset_mirna)[colnames(mirna),"title"]
colnames(mirna) <- sub("-", ".", sub("BreastTumor_", "", smplnames))
smplnames <- pData(eset_mrna)[colnames(mrna),"title"]
colnames(mrna) <- sub("-", ".", sub("BreastTumor_", "", smplnames))
rm(smplnames)

# Convert probe accession IDs to gene names (use Anduril IDConvert for convenience)
CSV.write(Table(annot_mrna), file.path(DATADIR,"mrna_annot.csv"))
CSV.write(mrna, file.path(DATADIR,"mrna.csv"), first.cell="ID")
system("cd /home/viljami/wrk/dippa-data; anduril run-component IDConvert -I csv=mrna.csv -I conversionTable=mrna_annot.csv -P conversionColumn=GENE_SYMBOL -P keyColumn=ID -P collapseNumeric=mean -P unique=true -O csv=mrna_genes.csv --java-heap 4000")
mrna <- as.matrix(read.delim(file.path(DATADIR,"mrna_genes.csv"), row.names=1))

# Reformat protein data
rownames(df_prot) <- df_prot$Gene
# Expand AKT1/2/3 and GSK3A/B genes in protein data
i <- match("AKT1/AKT2/AKT3", df_prot$Gene)
prot <- df_prot[c(1:i, i, i, (i+1):nrow(df_prot)), 3:ncol(df_prot)]
rownames(prot)[i:(i+2)] <- c("AKT1","AKT2","AKT3")
i <- match("GSK3A/GSK3B", rownames(prot))
prot <- prot[c(1:i, i, (i+1):nrow(prot)), ]
rownames(prot)[i:(i+1)] <- c("GSK3A","GSK3B")
prot <- as.matrix(prot)

# Normalization for regression (row-wise)
#norm.func <- function(row) { ()/() + SIGMA }





### OUTPUT ####

# Write out datafiles
CSV.write(file.path(DATADIR,"mrna_genes.csv"), mrna, first.cell="GENE_SYMBOL")
CSV.write(file.path(DATADIR,"mirna.csv"), mirna, first.cell="ID")
CSV.write(file.path(DATADIR,"protein.csv"), prot, first.cell="GENE_SYMBOL")
