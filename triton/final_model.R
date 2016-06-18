
## Simulate full model for chosen number of miRNA variables
#
# Fits the full model for given number of miRNA variables.
# The number of vars should be first computed with CV.
# This script then does projection var sel upto chosen number
# of vars, and then fits the full model for selected vars.


## INIT #####

WRKDIR  <- Sys.getenv("WRKDIR")
DATADIR <- file.path(WRKDIR,"dippa-data")
OUTDIR  <- file.path(WRKDIR,"dippa-analyysi","execute")

# Load necessary libraries
library(rstan)
library(rstanarm)
library(stats)
source(file.path(WRKDIR,"dippa-analyysi","triton","projection.R"))

# Check which job array ID we are at
jobi <- as.integer(commandArgs(TRUE)[1])

# Load data
print("Loading data...")
prot <- as.matrix(read.delim(file.path(DATADIR,"protein_normalized.csv"), row.names=1))
gene <- as.matrix(read.delim(file.path(DATADIR,"gene_normalized.csv"), row.names=1))
mirna <- as.matrix(read.delim(file.path(DATADIR,"mirna_normalized.csv"), row.names=1))
# Reformat data
d <- ncol(mirna)          #number of miRNA covariates
g <- colnames(prot)[jobi] #name of the gene
samples <- rownames(prot) #samplenames
n <- length(samples)      #num of samples
P <- as.numeric(prot[samples,g]) #protein expr vector
G <- as.numeric(gene[samples,g]) #gene expr vector
M <- mirna[samples,]             #miRNA expr matrix
# Get chosen number of miRNA variables
n_vars <- read.delim(file.path(DATADIR, "n_chosen_variables.csv"))
n_vars <- n_vars[match(g,n_vars[,1]),2]
if(is.na(n_vars)) n_vars <- 0
# Cleanup
rm(prot,gene,mirna,samples)

## Parameters
# Params for CV and simulation
model <- file.path(WRKDIR,"dippa-analyysi","stan","HS.stan")
nu <- 3.0 #parameter for hyperpriors (student-t degrees of freedom)
pn <- 13.75 #assumed number of meaningful covars, used for variance of tau prior, set small for more restrictive prior
n_iter <- 1000
n_chains <- 4
n_proj_samples <- 2000 #num of simulation samples to use for projection prediction
multicore <- FALSE
# Save params for reference later
params <- list(model=model, n=n, d=d, nu=nu, pn=pn, n_iter=n_iter, n_chains=n_chains, n_proj_samples=n_proj_samples, n_vars=n_vars)

# Output files
out_file <- file.path(OUTDIR,sprintf("finalmodel-%d-%s.rda",jobi,g))

# Use rstan multicore options if set
if(multicore) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = min(n_chains, parallel::detectCores()))
}


fit <- NA #last fit by rstan

## VARIABLE SELECTION #####

# Fit standard lm models for comparison
lm.fits <- list()

# Fit full model for gene variable only
print("Fitting gene only model...")
fit.gene <- stan_glm(P ~ G, prior=normal(0,5), prior_intercept=normal(0,5), iter=n_iter, chains=n_chains)
posterior.gene <- fit.gene$stan_summary
r2 <- 1 - var(residuals(fit.gene))/var(P)
r2.gene <- c( r2, 1-(1-r2)*(n-1)/(n-2) )
names(r2.gene) <- c("R2","R2.adjusted")

lm.fits[["gene_only"]] <- lm(P ~ G)

if(n_vars > 0) {

    # Fit full model for variable selection
    print("Fitting full model for projection variable selection...")
    datalist <- list(G=G, P=P, M=M, n=n, d=d, nu=nu, pn=pn)
    fit <- stan(file=model, data=datalist, iter=n_iter, chains=n_chains, fit=fit)

    ## Do the variable selection search
    # Get posterior samples of weights for full model
    e <- extract(fit)
    w <- rbind(e$w0, e$wg, t(e$w)) # stack the intercept and gene and miRNA weights
    sigma2 <- e$sigma^2
    # Take a sample of the weight posterior samples to improve projection speed
    isamp  <- sample.int(ncol(w), n_proj_samples)
    w.s      <- w[,isamp]
    sigma2.s <- sigma2[isamp]
    # Stack expression data and add column of 1s
    x <- cbind(rep(1,n), G, M)
    y <- P

    # Do the projection predictive variable selection!
    print("Doing projection variable selection...")
    spath <- lm_fprojsel(w.s, sigma2.s, x, n_vars+2)

    # Get the chosen miRNA vars
    chosen.mirnas <- colnames(x)[spath$chosen[3:length(spath$chosen)]]
    mirna.i <- match(chosen.mirnas, colnames(M))

    # Fit final full model for chosen miRNA vars
    #print("Fitting final full model...")
    #datalist <- list(G=G, P=P, M=M[,mirna.i], n=n, d=length(mirna.i), nu=nu, pn=pn)
    #fit <- stan(file=model, data=datalist, iter=n_iter, chains=n_chains, fit=fit)

    # Use the projected weights to do prediction
    w <- spath$w[,,n_vars+2]
    ypred <- x %*% w
    #r2 <- 1 - colSums((y-ypred)^2)/sum((y-mean(y))^2)
    resid <- y-ypred
    r2 <- 1-apply(resid,2,var)/var(y)
    r2.adj <- 1 - (1-r2)*(n-1)/(n-(n_vars+1)-1)

    # Save summary of marginal posterior distributions of full model
    sry <- summary(fit, probs=c(.025,.1,.25,.5,.75,.9,.975))$summary
    ikeepvars <- grep("w|tau|lambda", rownames(sry))
    posterior <- sry[ikeepvars,]

    # lm fit for comparison
    my.formula <- as.formula(paste("P ~ G + ", paste(chosen.mirnas, collapse=" + "),sep=""))
    lm.fits[["full_model"]] <- lm(my.formula, data=as.data.frame(x))

    # Save results
    save(chosen.mirnas, posterior, e, r2, r2.adj, spath, fit.gene, posterior.gene, r2.gene, lm.fits, params, file=out_file)
} else {
    # Save only gene model
    save(fit.gene, posterior.gene, r2.gene, lm.fits, params, file=out_file)
}
