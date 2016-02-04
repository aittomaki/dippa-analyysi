
### INIT ####

# necessary libraries
library(rstan)


### For running under Anduril, use this block
# Data
prot <- as.matrix(table1[,2:ncol(table1)])
rownames(prot) <- table1[,1]
gene <- as.matrix(table2[,2:ncol(table2)])
rownames(gene) <- table2[,1]
mirna <- as.matrix(table3[,2:ncol(table3)])
rownames(mirna) <- table3[,1]
# Parameters
model_file <- param1
n_iter <- as.numeric(param2)
n_chains <- as.numeric(param3)
multicore <- as.boolean(param4)
gene_names <- unlist(strsplit(param5, ","))
mirna_names <- unlist(strsplit(param6, ","))
# Output files
fitted_models_file <- get.output(cf,'optOut1')
rm(optOut1)
table.out <- data.frame()


### For running outside Anduril, use this block
# # Load data
# DATADIR <- "~/wrk/dippa-data/"
# prot <- read.delim(file.path(DATADIR,"protein.csv"), row.names=1)
# gene <- read.delim(file.path(DATADIR,"mrna_genes.csv"), row.names=1)
# mirna <- read.delim(file.path(DATADIR,"mirna.csv"), row.names=1)
# # Parameters
# model_file <- "simple_priors.stan"
# gene_names <- c("BRAF","SYK")
# mirna_names <- c("hsa-miR-638","hsa-miR-671-5p","hsa-miR-107")
# n_iter <- 1000
# n_chains <- 4
# multicore <- true
# # Output files
# fitted_models_file <- "fitted_models_multivariate.rda"
# posteriors_file <- "posteriors_multivariate.csv"


if(length(gene_names) < 1)
    gene_names <- rownames(gene)
if(length(mirna_names) < 1)
    mirna_names <- rownames(mirna)
if(multicore) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
}



### MULTIVARIATE REGRESSION ####

samples <- colnames(prot)
fits <- list()
posteriors <- list()
fit <- NA

# Run through all (given) proteins
for(g in gene_names) {
    datalist <- list(N=length(samples), J=nrow(mirna), P=as.numeric((prot[g,samples])), M=t(as.matrix(mirna[,samples])), G=as.numeric(gene[g,samples]))
    fit <- stan(file=model_file, data=datalist, iter=n_iter, chains=n_chains, fit=fit, model_name=paste(g,"multi",sep="_"))

    fits <- c(fits, list(fit))
    names(fits)[length(fits)] <- paste(g,"multi",sep="_")

    post <- summary(fit)$summary
    post <- as.data.frame(cbind(Gene=rep(g,nrow(post)), miRNA=c(rownames(mirna),rep(NA,nrow(post)-nrow(mirna))), coef=rownames(post), post))
    posteriors <- c(posteriors, list(post))
    names(posteriors)[length(posteriors)] <- paste(g,"multi",sep="_")

    save(fits, file=fitted_models_file)
}
array.out <- posteriors
table.out <- do.call(rbind, posteriors)



### OUTPUT #####

# Save model lists
save(fits, posteriors, file=fitted_models_file)
# Output if not running under Anduril
if(exists("posteriors_file"))
    write.table(table.out, file=posteriors_file, sep="\t", row.names=FALSE, )
