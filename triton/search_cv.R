#
# cross-validation of the variable searching. fits the full model K=10 times
# and does the variable selection each time separately, while measuring the performance
# of the found models on the validation data.
#

FILEDIR <- "/triton/work/jaittoma/dippa-analyysi"

# necessary libraries
require(rstan)
sessionInfo()
source(file.path(FILEDIR,"triton","projection.R"))

# which array ID are we at
jobi <- as.integer(commandArgs(TRUE)[1])

# Load data
print("Loading data...")
DATADIR <- "/triton/work/jaittoma/dippa-data"
prot <- as.matrix(read.delim(file.path(DATADIR,"protein_normalized.csv"), row.names=1))
gene <- as.matrix(read.delim(file.path(DATADIR,"gene_normalized.csv"), row.names=1))
mirna <- as.matrix(read.delim(file.path(DATADIR,"mirna_normalized.csv"), row.names=1))
d <- ncol(mirna)
g <- colnames(prot)[jobi] # name of the gene
samples <- rownames(prot)
n <- length(samples)
P <- as.numeric(prot[samples,g])
G <- as.numeric(gene[samples,g])
M <- mirna[samples,]

# Parameters
model_file <- file.path(FILEDIR,"stan","shrinkage_prior.stan")
nu <- 2.0 # parameter for hyperpriors (student-t degrees of freedom)
n_iter <- 1000
n_chains <- 4
multicore <- FALSE
# Output files
OUTDIR <- file.path(FILEDIR,"execute")
out_file <- file.path(OUTDIR,sprintf("CV-%d-%s.rda",jobi,g))

# Set rstan multicore options if wished
if(multicore) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = n_chains)
}


# cross-validate the variable searching
cvk <- 5
lpd <- matrix(0, cvk, 100) # lpd for each validation set
mse <- matrix(0, cvk, 100) # mse for - '' -
lpdfull <- rep(0, cvk)
msefull <- rep(0, cvk)

fit <- list(NA)
spath <- list()

for (i in 1:cvk) {

	# form the training and validation sets
	ival <- seq(i,n,cvk)
	itr <- setdiff(1:n,ival)
	nval <- length(ival)
	ntr <- length(itr)
	
	# fit the full model
	print(sprintf("Fitting full model for fold %d/%d...",i,cvk))
	datalist <- list(G=G[itr], P=P[itr], M=M[itr,], n=ntr, d=d, nu=nu)
	fit[[i]] <- stan(file=model_file, data=datalist, iter=n_iter, chains=n_chains, fit=fit[[1]])
	e <- extract(fit[[i]])
	
	# perform the variable selection
	w <- rbind(e$w0, e$wg, t(e$w)) # stack the intercept and gene and miRNA weights
	sigma2 <- e$sigma^2
	# combine vector of ones, gene vector and miRNA matrix as input matrix
	xtr <- cbind(rep(1,ntr), G[itr], M[itr,])
	spath[[i]] <- lm_fprojsel(w, sigma2, xtr)
	
	# make predictions for the observations in the validation set
	xval <- cbind(rep(1,nval), G[ival], M[ival,])
	yval <- P[ival]
	for (k in 1:100) {

		# projected parameters
		submodel <- lm_proj(w, sigma2, xtr, spath[[i]]$chosen[1:k])
		wp <- submodel$w
		sigma2p <- submodel$sigma2
		
		# mean squared error
		ypred <- rowMeans(xval %*% wp)
		mse[i,k] <- mean((yval-ypred)^2)
		
		# mean log predictive density using the projected parameters
		pd <- dnorm(yval, xval %*% wp, sqrt(sigma2p))
		lpd[i,k] <- mean(log(rowMeans(pd)))
	}
	# calculate mse and lpd for full model with validation set
	ypredfull <- rowMeans(xval %*% w)
	msefull[i] <-  mean((yval-ypredfull)^2)
	pdfull <- dnorm(yval, xval %*% w, sqrt(sigma2))
	lpdfull[i] <- mean(log(rowMeans(pdfull)))
}

# Save results
save(lpd, mse, fit, spath, file=out_file)
