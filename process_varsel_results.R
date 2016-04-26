FILEDIR <- "/home/viljami/wrk/tritonfake/"

# Load necessary libraries
library(bayesboot)
libray(HDInterval)
library(ggplot2)

# Load data
nom <- "CV-5-ANXA1"
load(file.path(FILEDIR,"dippa-analyysi","execute",paste(nom,".rda",sep="")))

# Parameters
credible.level <- 0.95
a <- (1-credible.level)/2
Q <- qnorm(1-a)
n_boot <- 10000

# Compute reference MLPD's
mlpd.0 <- mean(lpd[,1])
mlpd.full <- mean(lpd.full)
U <- 0.05*(mlpd.0-mlpd.full) # decision limit for variable number to choose

# Compute MLPD and credible intervals for submodels
# Use Bayesian bootstrap to compute posterior interval
bb <- function(x, n_boot=n_boot, credMass=credible.level) {
	babo <- bayesboot(x, weighted.mean, R=n_boot, use.weights=T)
	highdi <- hdi(babo, credMass=credMass)
}
d.mlpd <- apply(lpd, 2, mean)
credint <- apply(lpd, 2, bb)
se <- sapply(lpd, 2, function(x) sqrt(var(x)/length(x)))
se <-  rbind(d.mlpd-Q*se, d.mlpd+Q*se)
d.mlpd <- t(rbind(d.mlpd, credint, se))
# Compute delta MLPD, change into a data frame
d.mlpd <- d.mlpd-mlpd.full
d.mlpd <- as.data.frame(d.mlpd)
names(d.mlpd) <- c("dMLPD","hdi.lower","hdi.upper","se.lower","se.upper")
mlpd$n <- 1:nrow(mlpd)

# Get chosen number of variables to include in final model
n.chosen <- min(which(mlpd$lower > U))
print(n.chosen)

# Save output
save(mlpd,file=file.path(FILEDIR,"dippa-analyysi","execute",paste(nom,"-mlpd",".rda",sep="")))

# Make a plot of dMLPD vs num of vars in submodel
theme_set(theme_bw())
g <- ggplot(mlpd, aes(n,dMLPD))
g <- g + geom_hline(yintercept = 0, linetype="dashed")
g <- g + geom_hline(yintercept = U, color="red")
g <- g + geom_errorbar(aes(ymin=hdi.lower, ymax=hdi.upper), color="gray", width=0)
g <- g + geom_errorbar(aes(ymin=se.lower, ymax=se.upper), color="red", width=0)
g <- g + geom_line()
g <- g + ggtitle(nom)
ggsave(file.path(FILEDIR,"dippa-analyysi","execute",paste(nom,".png",sep="")), g)

# Plot an example of the BB distribution
babo <- bayesboot(lpd[,4], weighted.mean, R=n_boot, use.weights=T))
png(file.path(FILEDIR,"dippa-analyysi","execute",paste(nom,"-BBdistr_m4",".png",sep=""))))
print(plot(babo))
dev.off()
