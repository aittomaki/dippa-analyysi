
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
mirna_names <- unlist(strsplit(param6, ","))
# Output files
fitted_models_file <- get.output(cf,'optOut1')
rm(optOut1)
table.out <- data.frame()


### For running outside Anduril, use this block
# # Load data
# DATADIR <- "~/wrk/dippa-data/"
# prot <- as.matrix(read.delim(file.path(DATADIR,"protein.csv"), row.names=1))
# gene <- as.matrix(read.delim(file.path(DATADIR,"mrna_genes.csv"), row.names=1))
# mirna <- as.matrix(read.delim(file.path(DATADIR,"mirna.csv"), row.names=1))
# # Parameters
# model_file <- "shrinkage_prior.stan"
# gene_names <- c("BRAF","SYK")
# mirna_names <- c("hsa.miR.638","hsa.miR.671.5p","hsa.miR.107")
# n_iter <- 1000
# n_chains <- 4
# multicore <- true
# # Output files
# fitted_models_file <- "fitted_models_univariate.rda"
# posteriors_file <- "posteriors_univariate.csv"


# Check that gene and mirna names are defined
if(!exists("gene_names") || length(gene_names) < 1)
    gene_names <- colnames(gene)
if(!exists("mirna_names") || length(mirna_names) < 1)
    mirna_names <- colnames(mirna)

# Set rstan multicore options if wished
if(multicore) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
}



### UNIVARIATE REGRESSION ####

samples <- rownames(prot)
fits <- list()
posteriors <- list()
fit <- NA

# Run through all (given) miRNA-protein pairs
for(g in gene_names) {
    for(m in mirna_names){
        datalist <- list(nu=3, n=length(samples), d=1, P=prot[samples,g], M=as.matrix(mirna[samples,m]), G=gene[samples,g])
        fit <- stan(file=model_file, data=datalist, iter=n_iter, chains=n_chains, fit=fit, model_name=paste(g,m,"uni",sep="_"))

        fits <- c(fits, list(fit))
        names(fits)[length(fits)] <- paste(g,m,"uni",sep="_")

        post <- summary(fit)$summary
        post <- as.data.frame(cbind(Gene=rep(g,nrow(post)), miRNA=rep(m,nrow(post)), coef=rownames(post), post))
        posteriors <- c(posteriors, list(post))
        names(posteriors)[length(posteriors)] <- paste(g,m,"uni",sep="_")

        # Save results so far (in case of crash)
        save(fits, file=fitted_models_file)
        
        # Plot the simulated posteriors for parameters
        plot.file <- file.path(document.dir,sprintf('%s-%s-%s-posteriors.png',instance.name,g,m))
        png(plot.file)
        print(plot(fit))
        dev.off()

        # Plot trace plot
        plot.file <- file.path(document.dir,sprintf('%s-%s-%s-trace.png',instance.name,g,m))
        png(plot.file)
        print(traceplot(fit))
        dev.off()
    }
}



### OUTPUT #####

array.out <- posteriors
table.out <- do.call(rbind, posteriors)
document.out <- ""

# Save model lists
save(fits, posteriors, file=fitted_models_file)
# Output if not running under Anduril
if(exists("posteriors_file"))
    write.table(table.out, file=posteriors_file, sep="\t", row.names=FALSE)

