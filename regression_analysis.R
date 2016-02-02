### INIT ####

# necessary libraries
library(rstan)


### For running under Anduril, use this block
prot <- table1
gene <- table2
mirna <- table3
model_file <- param1
gene_names <- unlist(strsplit(param2, ","))
mirna_names <- unlist(strsplit(param3, ","))


### For running outside Anduril, uncomment the following
# # path to data
# DATADIR <- "~/wrk/dippa-data/"
# # Load data
# prot <- read.delim(file.path(DATADIR,"protein.csv"), row.names=1)
# gene <- read.delim(file.path(DATADIR,"mrna_genes.csv"), row.names=1)
# mirna <- read.delim(file.path(DATADIR,"mirna.csv"), row.names=1)
# model_file <- "simple_priors.stan"
# gene_names <- c("BRAF","SYK")
# mirna_names <- c("hsa-miR-638","hsa-miR-671-5p","hsa-miR-107")




### UNIVARIATE ####

samples <- colnames(prot)
fitlist <- list()

# Run through all (given) miRNA-protein pairs
for(g in gene_names) {
    for(m in mirna_names){
        datalist <- list(N=length(samples), J=1, P=as.numeric(prot[g,samples]), M=t(as.matrix(mirna[m,samples])), G=as.numeric(gene[g,samples]))
        fit <- stan(file=model_file, data=datalist)
        fitlist <- c(fitlist, list(fit))
        names(fitlist)[length(fitlist)] <- paste(g,m,sep="_")
    }
}

# Extract the miRNA coefficient and interval for first model
beta1 <- summary(fitlist[[1]])$summary["beta[1]",c("mean","2.5%","97.5%")]

