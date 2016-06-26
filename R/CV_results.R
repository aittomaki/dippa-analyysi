
# Script for processing utility metrics from projection predictive CV
#
# Significance level (alpha) is fixed to 0.5 or 0.975 since the 0.95
# HDI is computed in CV. Other levels need recomputing of HDI
# (or gaussian approximation).

# Input, params and output
if(exists("param1")) { # Anduril
    RESULTDIR <- param1
    PLOTDIR <- document.dir
    U.factor <- as.numeric(unlist(strsplit(param2, ",")))
    redo.cred.int <- as.logical(param3)
    if(redo.cred.int == "") redo.cred.int <- FALSE
} else { # non-Anduril
    RESULTDIR <- "/home/viljami/wrk/cvresults"
    PLOTDIR <- "/home/viljami/wrk/cvresults/plots"
    OUTFILE <- "/home/viljami/wrk/cvresults/selected_varnums.csv"
    if(!dir.exists(PLOTDIR)) dir.create(PLOTDIR)
    U.factor <- c(0.05, 0.1, 0.2)
    redo.cred.int <- FALSE
}
# Params for recomputing credible interval for dMLPD
n_boot <- 5000


library(bayesboot)
library(xtable)
library(ggplot2)
theme_set(theme_bw())

# Get result file names
resfiles <- list.files(RESULTDIR, pattern = "*.rda")
# Get gene names from filenames
genes <- sub("CV-\\d+-(\\w+).rda", "\\1", resfiles)
# Sort into alphabetical order by gene
isort <- order(genes)
resfiles <- resfiles[isort]
genes <- genes[isort]

# Result matrix, chosen number of vars to include
varnums <- matrix(0, nrow=length(genes), ncol=4*length(U.factor))

# Helper function for getting lowest index higher than threshold
# Gives index-2 i.e. number of miRNAs!
get_n_miRNA <- function(x, thresh) {
    lowest <- NA
    lt <- which(x > thresh)
    if(length(lt) > 0)
        lowest <- max(min(lt)-2, 0)
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
            cred.int <- quantile(babo[,1], probs=c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975))
        }
        deltaLPD <- lpd - lpd.full.m
        cred.int <- t(apply(deltaLPD, 2, bb))
        util$bb0.025 <- cred.int[,"2.5%"]
        util$bb0.10 <- cred.int[,"10%"]
        util$bb0.25 <- cred.int[,"25%"]
        util$bb0.50 <- cred.int[,"50%"]
        util$bb0.75 <- cred.int[,"75%"]
        util$bb0.90 <- cred.int[,"90%"]
        util$bb0.975 <- cred.int[,"97.5%"]
        params$n_boot=n_boot
        save(util, lpd, lpd.full, se, se.full, U, spath, posterior, params, file=file.path(RESULTDIR,f))
    }

    # Compute selected number of miRNAs!
    for (j in 1:length(U)) {
        varnums[i,(j*4-3)] <- get_n_miRNA(util$bb0.50, U[j])
        varnums[i,(j*4-2)] <- get_n_miRNA(util$bb0.25, U[j])
        varnums[i,(j*4-1)] <- get_n_miRNA(util$bb0.10, U[j])
        varnums[i,(j*4)] <- get_n_miRNA(util$bb0.025, U[j])
    }

    # Plot var selection path and utility for current gene
    g <- ggplot(data=util, aes(x=n, y=dMLPD))
    g <- g + geom_hline(yintercept=0)
    g <- g + geom_hline(yintercept=U, color="red")
    g <- g + geom_line()
    g <- g + geom_errorbar(aes(ymin=bb0.025, ymax=bb0.975), width=0, color="grey60")
    g <- g + geom_errorbar(aes(ymin=bb0.25, ymax=bb0.75), width=0, color="grey40")
    g <- g + labs(title=gene, x="N variables", y=bquote(Delta~MLPD))

    plot.file <- file.path(PLOTDIR, sprintf("%s_CV_path.png",gene))
    ggsave(plot.file, g, height=5, width=7, dpi=600)
}

# Convert result matrix to df
varnums <- data.frame(genes, varnums)
clnms <- paste(rep(paste("U", U.factor, sep=""), each=length(U)), c("a0.50","a0.75","a0.90","a0.975"), sep="_")
names(varnums) <- c("gene", clnms)
varnums <- varnums[order(varnums$gene),]

# Output
table.out <- varnums
document.out <- print(xtable(varnums), print.results=F)
if(!exists("param1")) { # non-Anduril
    write.table(table.out, file=OUTFILE, sep="\t", row.names=F)
}
