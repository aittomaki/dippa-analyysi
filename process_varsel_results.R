WRKDIR <- Sys.getenv("WRKDIR")
DATADIR <- file.path(WRKDIR,"dippa-analyysi","execute")
RESULTDIR <- file.path(DATADIR,"processed")

# Load necessary libraries
library(bayesboot)
library(HDInterval)
library(ggplot2)

# Parameters
credible.level <- 0.95
a <- (1-credible.level)/2
Q <- qnorm(1-a)
n_boot <- 5000

# List result files
files <- list.files(DATADIR, pattern = "CV-.*.rda")
print(str(files))
# Go through results
for(f in files) {
# Load data
load(file.path(DATADIR,f))


# Compute reference MLPD's
mlpd.0 <- mean(lpd[,1])
mlpd.full <- mean(lpd.full)
U <- 0.05*(mlpd.0-mlpd.full) # decision limit for variable number to choose

# Compute MLPD and credible intervals for submodels
# Use Bayesian bootstrap to compute posterior interval
bb <- function(x, nboot=n_boot, credMass=credible.level) {
	babo <- bayesboot(x, weighted.mean, R=nboot, use.weights=T)
	highdi <- hdi(babo, credMass=credMass)
}
d.mlpd <- apply(lpd, 2, mean)
credint <- apply(lpd, 2, bb)
se <- apply(lpd, 2, function(x) sqrt(var(x)/length(x)))
se <- rbind(d.mlpd-Q*se, d.mlpd+Q*se)
d.mlpd <- t(rbind(d.mlpd, credint, se))
# Compute delta MLPD, change into a data frame
d.mlpd <- d.mlpd-mlpd.full
d.mlpd <- as.data.frame(d.mlpd)
names(d.mlpd) <- c("dMLPD","hdi.lower","hdi.upper","se.lower","se.upper")
d.mlpd$n <- 1:nrow(d.mlpd)

# Get chosen number of variables to include in final model
n.chosen <- min(which(d.mlpd$lower > U))
print(n.chosen)


# Make a plot of dMLPD vs num of vars in submodel
theme_set(theme_bw())
g <- ggplot(d.mlpd, aes(n,dMLPD))
g <- g + geom_hline(yintercept = 0, linetype="dashed")
g <- g + geom_hline(yintercept = U, linetype="dashed", color="red")
g <- g + geom_errorbar(aes(ymin=hdi.lower, ymax=hdi.upper), color="gray", width=0.4)
g <- g + geom_errorbar(aes(ymin=se.lower, ymax=se.upper), color="red", width=0, alpha=0.4)
g <- g + geom_line()
g <- g + ggtitle(f)
ggsave(file.path(RESULTDIR,paste(f,".png",sep="")), g)






## COMPUTE THE FOLD-WISE MLPD VERSION TO COMPARE!
mlpd.cv.0 <- mean(lpd.cv[,1]) #this is for foldwise mlpd's
mlpd.cv.full <- mean(lpd.cv.full) #--"--
U.cv <- 0.05*(mlpd.cv.0-mlpd.cv.full)
d.mlpd.cv <- apply(lpd.cv, 2, mean)
credint.cv <- apply(lpd.cv, 2, bb)
se.cv <- apply(lpd.cv, 2, function(x) sqrt(var(x)/length(x)))
se.cv <- rbind(d.mlpd.cv-Q*se.cv, d.mlpd.cv+Q*se.cv)
d.mlpd.cv <- t(rbind(d.mlpd.cv, credint.cv, se.cv))
# Compute delta MLPD, change into a data frame
d.mlpd.cv <- d.mlpd.cv-mlpd.cv.full
d.mlpd.cv <- as.data.frame(d.mlpd.cv)
names(d.mlpd.cv) <- c("dMLPD","hdi.lower","hdi.upper","se.lower","se.upper")
d.mlpd.cv$n <- 1:nrow(d.mlpd.cv)

# Get chosen number of variables to include in final model
n.chosen.cv <- min(which(d.mlpd.cv$lower > U))
print(n.chosen.cv)

# Make a plot of dMLPD vs num of vars in submodel
g <- ggplot(d.mlpd.cv, aes(n,dMLPD))
g <- g + geom_hline(yintercept = 0, linetype="dashed")
g <- g + geom_hline(yintercept = U.cv, linetype="dashed", color="red")
g <- g + geom_errorbar(aes(ymin=hdi.lower, ymax=hdi.upper), color="gray", width=0.4)
g <- g + geom_errorbar(aes(ymin=se.lower, ymax=se.upper), color="red", width=0, alpha=0.4)
g <- g + geom_line()
g <- g + ggtitle(paste(f,"fold-wise MLPDs"))
ggsave(file.path(RESULTDIR,paste(f,"-cv.png",sep="")), g)







# Save output
save(d.mlpd,d.mlpd.cv,file=file.path(RESULTDIR,paste(f,"-mlpd.rda",sep="")))



# Plot an example of the BB distribution
babo <- bayesboot(lpd[,10], weighted.mean, R=n_boot, use.weights=T)
png(file.path(RESULTDIR,paste(f,"-BBdistr_m10.png",sep="")))
print(plot(babo))
dev.off()

}
