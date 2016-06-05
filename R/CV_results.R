
# Script for processing utility metrics from projection predictive CV
#
# Significance level (alpha) is fixed to 0.5 or 0.975 since the 0.95
# HDI is computed in CV. Other levels need recomputing of HDI
# (or gaussian approximation).

# Input, params and output
if(exists("param1")) { # Anduril
    RESULTDIR <- param1
    PLOTDIR <- document.dir
    U.factor <- as.numeric(param2)
} else { # non-Anduril
    RESULTDIR <- "/home/viljami/wrk/cvresults"
    PLOTDIR <- "/home/viljami/wrk/cvresults/plots"
    OUTFILE <- "/home/viljami/wrk/cvresults/selected_varnums.csv"
    if(!dir.exists(PLOTDIR)) dir.create(PLOTDIR)
    U.factor <- 0.05
}
# Params for recomputing credible interval for dMLPD
redo.cred.int <- TRUE
n_boot <- 5000


library(bayesboot)
library(xtable)
library(reshape2)
library(ggplot2)
theme_set(theme_bw())

resfiles <- list.files(RESULTDIR, pattern = "*.rda")
# Get gene names from filenames
genes <- sub("CV-\\d+-(\\w+).rda", "\\1", resfiles)
# Result matrix, chosen number of vars to include
varnums <- matrix(0, nrow=length(genes), ncol=2)

# Helper function for getting lowest index higher than threshold
# Gives index-2 i.e. number of miRNAs!
get_n_miRNA <- function(x, thresh) {
    lowest <- NA
    lt <- which(x > thresh)
    if(length(lt) > 0)
        lowest <- max(min(lt)-2,0)
    lowest
}

for (i in 1:length(resfiles)) {
    f <- resfiles[i]
    gene <- genes[i]
    # Load CV results for current gene
    load(file.path(RESULTDIR,f))

    # Convert varnum from 1:n to 0:n-1 (intercept not a var)
    if(util$n[1] == 1)
        util$n <- 0:(nrow(util)-1)

    # Compute U (=decision threshold for varnum)
    U <- U.factor*mean(lpd[,1]-lpd.full)

    # Recompute credible interval
    if(redo.cred.int) {
        MAX_VARS <- ncol(lpd)
        lpd.full.m <- matrix(rep(lpd.full, MAX_VARS), ncol=MAX_VARS, byrow=F)
        bb <- function(x, nboot=n_boot) {
            babo <- bayesboot(x, weighted.mean, R=nboot, use.weights=T)
            cred.int <- quantile(babo[,1], probs=c(0.025, 0.50, 0.975))
        }
        deltaLPD <- lpd - lpd.full.m
        cred.int <- t(apply(deltaLPD, 2, bb))
        util$ci.lower <- cred.int[,"2.5%"]
        util$ci.upper <- cred.int[,"97.5%"]
        util$median <- cred.int[,"50%"]
    } else {
        util$ci.lower <- util$hdi.lower
        util$ci.upper <- util$hdi.upper
        util$median <- util$dMLPD
    }
    # Compute selected number of miRNAs!
    varnum.50 <- get_n_miRNA(util$median, U)
    varnum.975 <- get_n_miRNA(util$ci.lower, U)
    varnums[i,] <- c(varnum.50, varnum.975)

    g <- ggplot(data=util, aes(x=n, y=dMLPD))
    g <- g + geom_hline(yintercept=0)
    g <- g + geom_hline(yintercept=U, color="red")
    g <- g + geom_line()
    g <- g + geom_errorbar(aes(ymin=ci.lower,ymax=ci.upper), width=0, alpha=0.4)
    g <- g + labs(title=gene, x="N variables", y=bquote(Delta~MLPD))

    plot.file <- file.path(PLOTDIR, sprintf("%s_CV_path.png",gene))
    ggsave(plot.file, g, height=3, width=4, dpi=600)
}

# Convert result matrix to df
varnums <- data.frame(genes, varnums)
names(varnums) <- c("gene","a_0.50","a_0.975")
varnums <- varnums[order(varnums$gene),]

# Make a histogram of num of vars selected
d <- melt(varnums, id.vars="gene", value.name="N_variables", variable.name="confidence")
g <- ggplot(d, aes(x=N_variables, group=confidence))
g <- g + geom_histogram(aes(fill=confidence), alpha=0.4)
plot.file <- file.path(PLOTDIR, "ZZ_variable_number_hist.pdf")
ggsave(plot.file, g, height=3, width=4, dpi=600)

# Output
table.out <- varnums
document.out <- print(xtable(varnums))
if(!exists("param1")) { # non-Anduril
    write.table(table.out, file=OUTFILE, sep="\t", row.names=F)
}
