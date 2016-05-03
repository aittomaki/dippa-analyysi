# Script for processing utility metrics from projection predictive CV

WRKDIR <- Sys.getenv("WRKDIR")
DATADIR <- file.path(WRKDIR,"dippa-analyysi","execute")
RESULTDIR <- file.path(DATADIR,"processed")

# Load necessary libraries
library(bayesboot)
library(HDInterval)
library(ggplot2)

# Parameters
n_boot <- 5000
credible.level <- 0.95
a <- (1-credible.level)/2
Q <- qnorm(1-a)
U.factor <- 0.05

# List result files
files <- list.files(DATADIR, pattern = "CV-.*.rda")
str(files)

# Go through result files in DATADIR
for(f in files) {

    # Load data
    load(file.path(DATADIR,f))
    gene <- sub("CV-\\d+-(\\w+).rda","\\1", f)

    # Compute decision limit for num of covars to choose
    U <- U.factor*mean(lpd[,1]-lpd.full)

    # Make a matrix version of full model LPD
    MAX_VARS <- ncol(lpd)
    lpd.full <- matrix(rep(lpd.full, MAX_VARS), ncol=MAX_VARS, byrow=F)

    # Compute dMLPD and its interval for each submodel
    # Use Bayesian bootstrap to compute posterior interval
    bb <- function(x, nboot=n_boot, credMass=credible.level) {
        babo <- bayesboot(x, weighted.mean, R=nboot, use.weights=T)
        highdi <- hdi(babo, credMass=credMass)
    }
    dLPD <- lpd - lpd.full
    dMLPD <- colMeans(dLPD)
    cred.int <- apply(dLPD, 2, bb)
    # Also compute gaussian interval for comparison
    se.int <- apply(dLPD, 2, function(x) sqrt(var(x)/length(x)))
    se.int <- rbind(dMLPD-Q*se.int, dMLPD+Q*se.int)
    # Make a data frame of dMLPD's and intervals
    util <- t(rbind(1:MAX_VARS, dMLPD, cred.int, se.int))
    util <- as.data.frame(util)
    names(util) <- c("n","dMLPD","hdi.lower","hdi.upper","se.lower","se.upper")

    # Get chosen number of variables to include in final model
    n.chosen <- min(which(util$hdi.lower > U))
    print(n.chosen)


    # Make a plot of dMLPD vs num of vars in submodel
    theme_set(theme_bw())
    g <- ggplot(util, aes(n,dMLPD))
    g <- g + geom_hline(yintercept = 0, linetype="dashed")
    g <- g + geom_hline(yintercept = U, linetype="dashed", color="red")
    g <- g + geom_errorbar(aes(ymin=hdi.lower, ymax=hdi.upper), color="gray", width=0.4)
    g <- g + geom_errorbar(aes(ymin=se.lower, ymax=se.upper), color="red", width=0, alpha=0.4)
    g <- g + geom_line()
    g <- g + ggtitle(bquote(paste(.(gene), " half-Cauchy ", tau, " prior")))
    g <- g + ylab(bquote(Delta~MLPD))
    ggsave(file.path(RESULTDIR,paste(f,".png",sep="")), g)


    # Save output
    save(util,file=file.path(RESULTDIR,paste(f,"-mlpd.rda",sep="")))


    # Plot an example of the BB distribution
    babo <- bayesboot(dLPD[,10], weighted.mean, R=n_boot, use.weights=T)
    png(file.path(RESULTDIR,paste(f,"-BBdistr_m10.png",sep="")))
    print(plot(babo))
    dev.off()

}
