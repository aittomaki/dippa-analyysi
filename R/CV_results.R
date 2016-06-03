
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

library(xtable)
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

    # Compute U (=decision threshold for varnum)
    U <- U.factor*mean(lpd[,1]-lpd.full)

    # Compute selected number of miRNAs!
    varnum.50 <- get_n_miRNA(util$dMLPD, U)
    varnum.975 <- get_n_miRNA(util$hdi.lower, U)
    varnums[i,] <- c(varnum.50, varnum.975)

    g <- ggplot(data=util, aes(x=n, y=dMLPD))
    g <- g + geom_hline(yintercept=0)
    g <- g + geom_hline(yintercept=U, color="red")
    g <- g + geom_line()
    g <- g + geom_errorbar(aes(ymin=hdi.lower,ymax=hdi.upper), width=0, alpha=0.4)
    g <- g + ggtitle(gene)

    plot.file <- file.path(PLOTDIR, sprintf("%s_CV_path.png",gene))
    ggsave(plot.file, g, height=3, width=4, dpi=600)
}

# Convert result matrix to df
varnums <- data.frame(genes, varnums)
names(varnums) <- c("gene","a_0.50","a_0.975")
varnums <- varnums[order(varnums$gene),]

# Output
table.out <- varnums
document.out <- print(xtable(varnums))
if(!exists("param1")) { # non-Anduril
    write.table(table.out, file=OUTFILE, sep="\t", row.names=F)
}
