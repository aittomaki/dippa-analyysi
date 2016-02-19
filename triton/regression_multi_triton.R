## Bayesian regression analysis script
## This version is meant to be run in Triton cluster!
## The script should be run as:
##   Rscript --vanilla regression_multi_triton.R $SLURM_ARRAY_TASK_ID

### INIT ####

# necessary libraries
require(rstan)
sessionInfo()

# which array ID are we at
jobi <- as.integer(commandArgs(TRUE)[1])

# Load data
print("Loading data...")
DATADIR <- "/triton/work/jaittoma/dippa-data"
prot <- as.matrix(read.delim(file.path(DATADIR,"protein_normalized.csv"), row.names=1))
gene <- as.matrix(read.delim(file.path(DATADIR,"gene_normalized.csv"), row.names=1))
mirna <- as.matrix(read.delim(file.path(DATADIR,"mirna_normalized.csv"), row.names=1))

# Parameters
g <- colnames(prot)[jobi] # name of the gene
model_file <- "shrinkage_prior.stan"
nu <- 3.0 # parameter for hyperpriors (student-t degrees of freedom)
n_iter <- 1000
n_chains <- 4
multicore <- FALSE
# Output files
OUTDIR <- "./output"
out_file <- file.path(OUTDIR,sprintf("fit-%d-%s.rda",jobi,g))

# Set rstan multicore options if wished
if(multicore) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = n_chains)
}



### REGRESSION ####

samples <- rownames(prot)
fit <- NA

print("Fitting model...")
datalist <- list(nu=nu, n=length(samples), d=ncol(mirna), P=prot[samples,g], M=mirna[samples,], G=gene[samples,g])
fit <- stan(file=model_file, data=datalist, iter=n_iter, chains=n_chains, fit=fit, model_name=g)

post <- summary(fit)$summary
post <- as.data.frame(cbind(Gene=rep(g,nrow(post)), miRNA=c(rownames(mirna),rep(NA,nrow(post)-nrow(mirna))), coef=rownames(post), post))

# Save fitted model
save(fit, post, file=out_file)

# # Plot the simulated posteriors for parameters
# plot.file <- file.path(document.dir,sprintf('%s-%s-posteriors.png',instance.name,g))
# png(plot.file)
# print(plot(fit))
# dev.off()

# # Plot trace plot
# plot.file <- file.path(document.dir,sprintf('%s-%s-trace.png',instance.name,g))
# png(plot.file)
# print(traceplot(fit))
# dev.off()
