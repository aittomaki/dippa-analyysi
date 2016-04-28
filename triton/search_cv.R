## Cross-validation of the variable searching.
#
# Fits the full model K=10 times and does the variable selection separately
# for each CV fold, while measuring the performance of the found submodels
# on the validation data.
#
# Outputs performance metrics (MLPD, MSE) for each fold and
# posteriors for full models.

FILEDIR <- "/triton/work/jaittoma"

# Load necessary libraries
require(rstan)
source(file.path(FILEDIR,"dippa-analyysi","triton","projection.R"))
sessionInfo()

# Check which job array ID we are at
jobi <- as.integer(commandArgs(TRUE)[1])

# Load data
print("Loading data...")
DATADIR <- file.path(FILEDIR,"dippa-data")
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

# Parameters
modelfile <- file.path(FILEDIR,"dippa-analyysi","stan","shrinkage_prior.stan")
nu <- 3.0 # parameter for hyperpriors (student-t degrees of freedom)
n_iter <- 1000
n_chains <- 4
n_proj_samples <- 200 #num of simulation samples to use for projection prediction
MAX_VARS <- 50        #max num of covars to add into model in projection prediction
multicore <- FALSE

# Output files
OUTDIR <- file.path(FILEDIR,"dippa-analyysi","execute")
out_file <- file.path(OUTDIR,sprintf("CV-%d-%s.rda",jobi,g))

# Use rstan multicore options if set
if(multicore) {
	rstan_options(auto_write = TRUE)
	options(mc.cores = n_chains)
}





# Cross-validate the variable searching (to get good n of vars)
cvk <- 10
lpd <- matrix(0, n, MAX_VARS) # lpd for each validation sample in each CV fold for each submodel
se <- matrix(0, n, MAX_VARS)  # square error for - '' -
lpd.full <- numeric(n) # lpd for each sample for full model
se.full <- numeric(n)  # square error for - '' -
mlpd <- matrix(0, cvk, MAX_VARS) # lpd for each validation set
mse <- matrix(0, cvk, MAX_VARS) # mse for - '' -
mlpd.full <- numeric(cvk) # lpd for each sample for full model
mse.full <- numeric(cvk)  # square error for - '' -


fit <- NA
posterior <- list()
spath <- list()

for (i in 1:cvk) {

	## Form the training and validation sets
	ival <- seq(i,n,cvk)
	itr <- setdiff(1:n,ival)
	nval <- length(ival)
	ntr <- length(itr)
	
	## Fit the full model for this CV-fold
	print(sprintf("Fitting full model for fold %d/%d...",i,cvk))
	datalist <- list(G=G[itr], P=P[itr], M=M[itr,], n=ntr, d=d, nu=nu)
	fit <- stan(file=modelfile, data=datalist, iter=n_iter, chains=n_chains, fit=fit)
	# save summary of posterior distributions
	sry <- summary(fit, probs=c(.025,.1,.25,.5,.75,.9,.975))$summary
	ikeepvars <- grep("w|tau|lambda", rownames(sry))
	posterior[[i]] <- sry[ikeepvars,]
	
	## Do the variable selection search
	# get weights for full model
	e <- extract(fit)
	w <- rbind(e$w0, e$wg, t(e$w)) # stack the intercept and gene and miRNA weights
	sigma2 <- e$sigma^2
	# form training and validation data
	xtr <- cbind(rep(1,ntr), G[itr], M[itr,])
	xval <- cbind(rep(1,nval), G[ival], M[ival,])
	yval <- P[ival]
	# calculate lpd and se for full model
	pd <- dnorm(yval, xval %*% w, sqrt(sigma2))
	lpd.full[ival] <- log(rowMeans(pd))
	mlpd.full[i] <- mean(log(rowMeans(pd)))
	ypred <- rowMeans(xval %*% w)
	se.full[ival] <- (yval-ypred)^2
	mse.full[i] <- mean((yval-ypred)^2)

	# take a smaller sample of the weight samples to improve projection speed
	isamp  <- sample.int(ncol(w), n_proj_samples)
	w      <- w[,isamp]
	sigma2 <- sigma2[isamp]

	# Do the projection predictive variable selection!
	print(sprintf("Doing variable selection for fold %d/%d...",i,cvk))
	spath[[i]] <- lm_fprojsel(w, sigma2, xtr, MAX_VARS)

	## Calculate lpd and se for each submodel along the selection path
	## by making predictions for the observations in the validation set
	for (k in 1:MAX_VARS) {

		# projected parameters
		submodel <- lm_proj(w, sigma2, xtr, spath[[i]]$chosen[1:k])
		wp       <- submodel$w
		sigma2p  <- submodel$sigma2
		
		# squared error
		ypred <- rowMeans(xval %*% wp)
		se[ival,k] <- (yval-ypred)^2
		mse[i,k] <- mean((yval-ypred)^2)
		
		# log predictive density using the projected parameters
		pd <- dnorm(yval, xval %*% wp, sqrt(sigma2p))
		lpd[ival,k] <- log(rowMeans(pd))
		mlpd[i,k] <- mean(log(rowMeans(pd)))
	}
}

# Save results
save(lpd, lpd.full, mlpd, mlpd.full, se, se.full, mse, mse.full, spath, posterior, file=out_file)
