### INIT ####

# necessary libraries
library(rstan)

# path to data
DATADIR <- "~/wrk/dippa-data/"

# Constants for analysis
SIGMA <- 0.1  # used for guaranteeing positivity in normalization

# Load data
mrna <- read.delim(file.path(DATADIR,"mrna_genes.csv"), row.names=1)
mirna <- read.delim(file.path(DATADIR,"mirna.csv"), row.names=1)
prot <- read.delim(file.path(DATADIR,"protein.csv"), row.names=1)


### UNIVARIATE ####

samples <- colnames(prot)
gene_names <- c("BRAF","SYK")
mirna_names <- c("hsa-miR-638","hsa-miR-671-5p","hsa-miR-107")
fitlist <- list()

for(g in gene_names) {
    for(m in mirna_names){
        datalist <- list(N=length(samples), J=1, P=as.numeric(prot[g,samples]), M=t(as.matrix(mirna[m,samples])), G=as.numeric(mrna[g,samples]))
        fit <- stan(file="simple_priors.stan", data=datalist)
        fitlist <- c(fitlist, list(fit))
        names(fitlist)[length(fitlist)] <- paste(g,m,sep="_")
    }
}

# Extract the miRNA coefficient and interval for first model
beta1 <- summary(fitlist[[1]])$summary["beta[1]",c("mean","2.5%","97.5%")]

