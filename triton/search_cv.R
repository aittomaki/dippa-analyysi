## Cross-validation of the variable searching.
#
# Fits the full model K=10 times and does the variable selection separately
# for each CV fold, while measuring the performance of the found submodels
# on the validation data.
#
# Outputs performance metrics (MLPD, MSE) for each fold and
# posteriors for full models.

WRKDIR  <- Sys.getenv("WRKDIR")
DATADIR <- file.path(WRKDIR,"dippa-data")
OUTDIR  <- file.path(WRKDIR,"dippa-analyysi","execute")

# Load necessary libraries
library(rstan)
library(stats)
library(bayesboot)
library(Matrix)
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
# Cleanup
rm(prot,gene,mirna,samples)

## Parameters
# Params for CV and simulation
model <- file.path(WRKDIR,"dippa-analyysi","stan","HS.stan")
nu <- 3.0 #parameter for hyperpriors (student-t degrees of freedom)
pn <- 13.75 #assumed number of meaningful covars, used for variance of tau prior, set small for more restrictive prior
n_iter <- 1000
n_chains <- 4
n_proj_samples <- 1000 #num of simulation samples to use for projection prediction
MAX_VARS <- 200        #max num of covars to add into model in projection prediction
multicore <- FALSE
# Params for dLPD intervals
n_boot <- 5000
conf_level <- 0.95
U.factor <- 0.05
# Save params for reference later
params <- list(model=model, n=n, d=d, nu=nu, pn=pn, n_iter=n_iter, n_chains=n_chains, n_proj_samples=n_proj_samples, MAX_VARS=MAX_VARS, n_boot=n_boot, conf_level=conf_level, U.factor=U.factor)

# Output files
out_file <- file.path(OUTDIR,sprintf("CV-%d-%s.rda",jobi,g))

# Use rstan multicore options if set
if(multicore) {
	rstan_options(auto_write = TRUE)
	options(mc.cores = n_chains)
}





## Cross-validate the variable searching (to get good n of vars)
cvk <- 10 #number of folds in CV
# Variables for results
lpd <- matrix(0, n, MAX_VARS) #lpd for each validation sample in each CV fold for all submodels
se <- matrix(0, n, MAX_VARS)  #squared error for - '' -
lpd.full <- numeric(n) #lpd for each sample for full model
se.full <- numeric(n)  #squared error for - '' -
posterior <- list() #posterior distrs for weights of full model for each CV-fold
spath <- list() #variable forward selection paths for each CV-fold

fit <- NA #last fit by rstan

for (i in 1:cvk) {

	# Form training and validation set indices
	ival <- seq(i,n,cvk)
	itr <- setdiff(1:n,ival)
	nval <- length(ival)
	ntr <- length(itr)

	# Fit full model for this fold
	print(sprintf("Fitting full model for fold %d/%d...",i,cvk))
	datalist <- list(G=G[itr], P=P[itr], M=M[itr,], n=ntr, d=d, nu=nu, pn=pn)
	fit <- stan(file=model, data=datalist, iter=n_iter, chains=n_chains, fit=fit)

	# Save summary of marginal posterior distributions of full model
	sry <- summary(fit, probs=c(.025,.1,.25,.5,.75,.9,.975))$summary
	ikeepvars <- grep("w|tau|lambda", rownames(sry))
	posterior[[i]] <- sry[ikeepvars,]

	## Do the variable selection search
	# Get posterior samples of weights for full model
	e <- extract(fit)
	w <- rbind(e$w0, e$wg, t(e$w)) # stack the intercept and gene and miRNA weights
	sigma2 <- e$sigma^2

	# Form training and validation data (adding column of 1's for intercept)
	xtr <- cbind(rep(1,ntr), G[itr], M[itr,])
	xval <- cbind(rep(1,nval), G[ival], M[ival,])
	yval <- P[ival]

	# Calculate lpd and se for full model
	pd <- dnorm(yval, xval %*% w, sqrt(sigma2))
	lpd.full[ival] <- log(rowMeans(pd))
	ypred <- rowMeans(xval %*% w)
	se.full[ival] <- (yval-ypred)^2

	# Take a sample of the weight posterior samples to improve projection speed
	isamp  <- sample.int(ncol(w), n_proj_samples)
	w      <- w[,isamp]
	sigma2 <- sigma2[isamp]

	# Do the projection predictive variable selection!
	print(sprintf("Doing variable selection for fold %d/%d...",i,cvk))
	spath[[i]] <- lm_fprojsel(w, sigma2, xtr, MAX_VARS)

	# Calculate lpd and se for each submodel along the selection path
	# by making predictions for the observations in the validation set
	for (k in 1:MAX_VARS) {

		# projected parameters
		submodel <- lm_proj(w, sigma2, xtr, spath[[i]]$chosen[1:k])
		wp       <- submodel$w
		sigma2p  <- submodel$sigma2

		# squared error
		ypred <- rowMeans(xval %*% wp)
		se[ival,k] <- (yval-ypred)^2

		# log predictive density using the projected parameters
		pd <- dnorm(yval, xval %*% wp, sqrt(sigma2p))
		lpd[ival,k] <- log(rowMeans(pd))
	}
	# Drop the projected weight matrices from spath
	spath[[i]]$w <- NULL
}

## Compute utility differences for submodels

# Make a matrix of full model utility
lpd.full.m <- matrix(rep(lpd.full, MAX_VARS), ncol=MAX_VARS, byrow=F)

# Compute deltaLPD and its interval for each submodel
# Use Bayesian bootstrap to compute posterior interval
bb <- function(x, nboot=n_boot, credMass=conf_level) {
    babo <- bayesboot(x, weighted.mean, R=nboot, use.weights=T)
    cred.int <- quantile(babo[,1], probs=c(0.025, 0.25, 0.50, 0.75, 0.975))
}
deltaLPD <- lpd - lpd.full.m
deltaMLPD <- colMeans(deltaLPD)
cred.int <- apply(deltaLPD, 2, bb)
# Also compute gaussian interval for comparison
a <- (1-conf_level)/2
Q <- qnorm(1-a)
se.int <- apply(deltaLPD, 2, function(x) sqrt(var(x)/length(x)))
se.int <- rbind(deltaMLPD-Q*se.int, deltaMLPD+Q*se.int)
# Make a data frame of dMLPD's and intervals
util <- t(rbind(1:MAX_VARS, deltaMLPD, cred.int, se.int))
util <- as.data.frame(util)
names(util) <- c("n","dMLPD","bb0.025","bb0.25","bb0.50","bb0.75","bb0.975","se.lower","se.upper")


# Compute decision limit for num of covars to choose (just for convenience)
U <- U.factor*mean(lpd[,1]-lpd.full)

# Save results
save(util, lpd, lpd.full, se, se.full, U, spath, posterior, params, file=out_file)
