
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
multicore <- as.logical(param4)
gene_names <- unlist(strsplit(param5, ","))
# Output files
fitted_models_file <- get.output(cf,'optOut1')
rm(optOut1)
table.out <- data.frame()


### For running outside Anduril, use this block
# # Load data
# DATADIR <- "./data/"
# prot <- as.matrix(read.delim(file.path(DATADIR,"protein_normalized.csv"), row.names=1))
# gene <- as.matrix(read.delim(file.path(DATADIR,"gene_normalized.csv"), row.names=1))
# mirna <- as.matrix(read.delim(file.path(DATADIR,"mirna_normalized.csv"), row.names=1))
# # Parameters
# jobi <- commandArgs(TRUE)[1]
# model_file <- "shrinkage_prior.stan"
# gene_names <- colnames(prot)[jobi]
# n_iter <- 1000
# n_chains <- 4
# multicore <- FALSE
# # Output files
# fitted_models_file <- "fitted_models_multivariate.rda"
# posteriors_file <- "posteriors_multivariate.csv"


# Check that gene and mirna names are defined
if(!exists("gene_names") || length(gene_names) < 1)
    gene_names <- colnames(gene)


# Set rstan multicore options if wished
if(multicore) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = n_chains)
}



### MULTIVARIATE REGRESSION ####

samples <- rownames(prot)
fits <- list()
posteriors <- list()
fit <- NA

# Run through all (given) proteins
for(g in gene_names) {
    datalist <- list(nu=3, n=length(samples), d=ncol(mirna), P=prot[samples,g], M=mirna[samples,], G=gene[samples,g])
    fit <- stan(file=model_file, data=datalist, iter=n_iter, chains=n_chains, fit=fit, model_name=paste(g,"multi",sep="_"))

    fits <- c(fits, list(fit))
    names(fits)[length(fits)] <- paste(g,"multi",sep="_")

    post <- summary(fit)$summary
    post <- as.data.frame(cbind(Gene=rep(g,nrow(post)), miRNA=c(rownames(mirna),rep(NA,nrow(post)-nrow(mirna))), coef=rownames(post), post))
    posteriors <- c(posteriors, list(post))
    names(posteriors)[length(posteriors)] <- paste(g,"multi",sep="_")

    # Save results so far (in case of crash)
    save(fits, file=fitted_models_file)
    
    # Plot the simulated posteriors for parameters
    plot.file <- file.path(document.dir,sprintf('%s-%s-posteriors.png',instance.name,g))
    png(plot.file)
    print(plot(fit))
    dev.off()

    # Plot trace plot
    plot.file <- file.path(document.dir,sprintf('%s-%s-trace.png',instance.name,g))
    png(plot.file)
    print(traceplot(fit))
    dev.off()
}
array.out <- posteriors
table.out <- do.call(rbind, posteriors)



### OUTPUT #####

array.out <- posteriors
table.out <- do.call(rbind, posteriors)
document.out <- ""

# Save model lists
save(fits, posteriors, file=fitted_models_file)
# Output if not running under Anduril
if(exists("posteriors_file"))
    write.table(table.out, file=posteriors_file, sep="\t", row.names=FALSE)
