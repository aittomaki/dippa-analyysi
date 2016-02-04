
### INIT ####

# necessary libraries
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


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
gene_names <- unlist(strsplit(param4, ","))
mirna_names <- unlist(strsplit(param5, ","))
# Output files
fitted_models_file <- get.output(cf,'optOut1')
rm(optOut1)


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
# # Output files
# fitted_models_file <- "fitted_models.rda"
# posteriors_uni_file <- "posteriors_univariate.csv"


if(length(gene_names) < 1)
    gene_names <- rownames(gene)
if(length(mirna_names) < 1)
    mirna_names <- rownames(mirna)
table.out <- data.frame()


### UNIVARIATE ####

samples <- colnames(prot)
fits.uni <- list()
posteriors <- matrix(nrow=0, ncol=length(valuenames))
names <- matrix(nrow=0, ncol=3)

# Run through all (given) miRNA-protein pairs
for(g in gene_names) {
    for(m in mirna_names){
        datalist <- list(N=length(samples), J=1, P=as.numeric((prot[g,samples])), M=as.matrix(mirna[m,samples]), G=as.numeric(gene[g,samples]))
        fit <- stan(file=model_file, data=datalist, model_name=paste(g,m,"uni",sep="_"), iter=n_iter, chains=n_chains)

        fits.uni <- c(fits.uni, list(fit))
        names(fits.uni)[length(fits.uni)] <- paste(g,m,"uni",sep="_")

        post_new <- summary(fit)$summary
        posteriors <- rbind(posteriors, post_new)
        names_new <- cbind(rep(g,nrow(post_new)),rep(m,nrow(post_new)),rownames(post_new))
        names <- rbind(names, names_new)
    }
}
table.out <- as.data.frame(cbind(names, posteriors))
names(table.out)[1:3] <- c("Gene","miRNA","coef")




### MULTIVARIATE ####

fits.multi <- list()
posteriors.multi <- list()

# Run through all (given) proteins
for(g in gene_names) {
    datalist <- list(N=length(samples), J=nrow(mirna), P=as.numeric((prot[g,samples])), M=t(as.matrix(mirna[,samples])), G=as.numeric(gene[g,samples]))
    fit <- stan(file=model_file, data=datalist, model_name=paste(g,"multi",sep="_"), iter=n_iter, chains=n_chains)

    fits.multi <- c(fits.multi, list(fit))
    names(fits.multi)[length(fits.multi)] <- paste(g,"multi",sep="_")

    post_new <- summary(fit)$summary
    post_new <- as.data.frame(cbind(coef=rownames(post_new), post_new))
    posteriors.multi <- c(posteriors.multi, list(post_new))
    names(posteriors.multi)[length(posteriors.multi)] <- paste(g,"multi",sep="_")

    save(fits.multi, file=fitted_models_file)
}
array.out <- posteriors.multi



### OUTPUT #####

# Save model lists
save(fits.uni, fits.multi, posteriors.multi, posteriors, file=fitted_models_file)
# Output if not running under Anduril
if(exists("posteriors_uni_file"))
    write.table(table.out, file=posteriors_uni_file, sep="\t", row.names=FALSE, )
