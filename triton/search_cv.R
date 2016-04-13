#
# cross-validation of the variable searching. fits the full model K=10 times
# and does the variable selection each time separately, while measuring the performance
# of the found models on the validation data.
#

FILEDIR <- "/triton/work/jaittoma"

# necessary libraries
require(rstan)
source(file.path(FILEDIR,"dippa-analyysi","triton","projection.R"))
sessionInfo()

# which array ID are we at
jobi <- as.integer(commandArgs(TRUE)[1])

# Load data
print("Loading data...")
DATADIR <- file.path(FILEDIR,"dippa-data")
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
model_file <- file.path(FILEDIR,"dippa-analyysi","stan","shrinkage_prior.stan")
nu <- 3.0 # parameter for hyperpriors (student-t degrees of freedom)
n_iter <- 1000
n_chains <- 4
n_proj_samples <- 200
multicore <- FALSE
do.plots <- FALSE
MAX_VARS <- 50
# Output files
OUTDIR <- file.path(FILEDIR,"dippa-analyysi","execute")
out_file <- file.path(OUTDIR,sprintf("CV-%d-%s.rda",jobi,g))

# Set rstan multicore options if wished
if(multicore) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = n_chains)
}


# cross-validate the variable searching
cvk <- 10
lpd <- matrix(0, cvk, MAX_VARS) # lpd for each validation set
mse <- matrix(0, cvk, MAX_VARS) # mse for - '' -
lpdfull <- rep(0, cvk)
msefull <- rep(0, cvk)

fit <- NA
posterior <- list()
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
	fit <- stan(file=model_file, data=datalist, iter=n_iter, chains=n_chains, fit=fit)
	e <- extract(fit)
	# save posterior summary
	sry <- summary(fit, probs=c(.025,.1,.25,.5,.75,.9,.975))$summary
	ikeepvars <- grep("w|tau|lambda", rownames(sry))
	posterior[[i]] <- sry[ikeepvars,]
	
	# perform the variable selection
	w <- rbind(e$w0, e$wg, t(e$w)) # stack the intercept and gene and miRNA weights
	sigma2 <- e$sigma^2
	# take a sample of the weight samples to improve projection speed
	isamp <- sample.int(ncol(w),n_proj_samples)
	w <- w[,isamp]
	sigma2 <- sigma2[isamp]
	# combine vector of ones, gene vector and miRNA matrix as input matrix
	xtr <- cbind(rep(1,ntr), G[itr], M[itr,])
	spath[[i]] <- lm_fprojsel(w, sigma2, xtr, MAX_VARS)
	
	# make predictions for the observations in the validation set
	xval <- cbind(rep(1,nval), G[ival], M[ival,])
	yval <- P[ival]
	# calculate mse and lpd for full model with validation set
	ypredfull <- rowMeans(xval %*% w)
	msefull[i] <-  mean((yval-ypredfull)^2)
	pdfull <- dnorm(yval, xval %*% w, sqrt(sigma2))
	lpdfull[i] <- mean(log(rowMeans(pdfull)))
	# then for each submodel along the selection path
	for (k in 1:MAX_VARS) {

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
}

# Save results
save(lpd, mse, lpdfull, msefull, spath, posterior, file=out_file)

if(do.plots) {
    # Plot a bit
    #png(file=file.path(OUTDIR,sprintf("CV-%d-%s.png",jobi,g)))
    #layout(matrix(c(1,2)))
    #plot(colMeans(lpd-matrix(rep(lpdfull,ncol(lpd)),nrow=nrow(lpd))), xlab="nvar", ylab="dLPD")
    #title("dLPD")
    #abline(h=0);
    #plot(colMeans(mse-matrix(rep(msefull,ncol(mse)),nrow=nrow(mse))), xlab="nvar", ylab="dMSE")
    #title("dMSE")
    #abline(h=0)
    #dev.off()

    # Better plots
    library(ggplot2)
    library(reshape)
    source("multiplot.R")
    summaryfun <- "mean_se"
    theme_set(theme_bw())
    dlpd <- lpd-matrix(rep(lpdfull,ncol(lpd)),nrow=nrow(lpd))
    dlpd <- melt(dlpd, varnames=c("CV","nvar"))
    p1 <- ggplot(dlpd, aes(nvar, value)) + geom_point(size=0.3)
    p1 <- p1 + stat_summary(fun.data=summaryfun, color="red") + geom_abline(slope=0)
    p1 <- p1 + labs(title=expression(Delta~LPD), x="variables", y=expression(Delta~LPD))
    dmse <- mse-matrix(rep(msefull,ncol(mse)),nrow=nrow(mse))
    dmse <- melt(dmse, varnames=c("CV","nvar"))
    p2 <- ggplot(dmse, aes(nvar, value)) + geom_point(size=0.3)
    p2 <- p2 + stat_summary(fun.data=summaryfun, color="red") + geom_abline(slope=0)
    p2 <- p2 + labs(title=expression(Delta~MSE), x="variables", y=expression(Delta~MSE))
    png(file=file.path(OUTDIR,sprintf("CV-%d-%s.png",jobi,g)),height=600,width=800)
    print(multiplot(p1,p2))
    dev.off()
}
