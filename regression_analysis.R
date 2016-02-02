### INIT ####

# necessary libraries
library(rstan)


### For running under Anduril, use this block
prot <- as.matrix(table1[,2:ncol(table1)])
rownames(prot) <- table1[,1]
gene <- as.matrix(table2[,2:ncol(table2)])
rownames(gene) <- table2[,1]
mirna <- as.matrix(table3[,2:ncol(table3)])
rownames(mirna) <- table3[,1]
model_file <- param1
gene_names <- unlist(strsplit(param2, ","))
mirna_names <- unlist(strsplit(param3, ","))

print(paste("gene",class(gene),paste(dim(gene),collapse=",")))
print(paste("prot",class(prot),paste(dim(prot),collapse=",")))
print(paste("mirna",class(mirna),paste(dim(mirna),collapse=",")))

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
varnames <- c("beta[1]","gamma")
valuenames <- c("mean","2.5%","97.5%","n_eff","Rhat")
fitlist <- list()
posteriors <- matrix(nrow=0, ncol=length(valuenames))
names <- matrix(nrow=0, ncol=3)

# Run through all (given) miRNA-protein pairs
for(g in gene_names) {
    for(m in mirna_names){
        datalist <- list(N=length(samples), J=1, P=as.numeric((prot[g,samples])), M=as.matrix(mirna[m,samples]), G=as.numeric(gene[g,samples]))
        fit <- stan(file=model_file, data=datalist)
        fitlist <- c(fitlist, list(fit))
        names(fitlist)[length(fitlist)] <- paste(g,m,sep="_")
        post_new <- summary(fit)$summary[varnames, valuenames]
        posteriors <- rbind(posteriors, post_new)
        names_new <- cbind(rep(g,length(varnames)),rep(m,length(varnames)),varnames)
        names <- rbind(names, names_new)
    }
}
table.out <- as.data.frame(cbind(names, posteriors))
names(table.out) <- c("Gene","miRNA","coef",valuenames)

# Extract the miRNA coefficient and interval for first model
#beta1 <- summary(fitlist[[1]])$summary[c("beta[1]","gamma"),c("mean","2.5%","97.5%","n_eff","Rhat")]

