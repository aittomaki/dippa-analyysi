### Script for reading simulation results for full models

library(rstan)
library(reshape)

if(exists("param1")) { # Anduril
    RESULTDIR <- param1
    PLOTDIR <- document.dir
    OUTFILE <- get.output(cf,'optOut1')
    rm(optOut1)
} else { # non-Anduril
    RESULTDIR <- "/home/viljami/output"
    PLOTDIR <- "/home/viljami/output/plots"
    OUTFILE <- "/home/viljami/output/sigs.rda"
}

sigs95 <- list()
sigs80 <- list()
sigs50 <- list()

files <- list.files(RESULTDIR, pattern = "fit-.*.rda")
genes <- sub("fit-\\d+-(\\w+).rda","\\1", files)
# Sort into gene alphabetical order
i <- order(genes)
files <- files[i]
genes <- genes[i]

for(i in 1:length(files)) {
    f <- files[i]
    g <- genes[i]

    load(file.path(RESULTDIR,f))
    sry <- summary(fit, probs=c(.025,.1,.25,.75,.9,.975))$summary

    sig <- which(sign(sry[,"2.5%"]) == sign(sry[,"97.5%"]))
    sig <- sig[grep("w",names(sig))]
    sigs95[[g]] <- sig

    sig <- which(sign(sry[,"10%"]) == sign(sry[,"90%"]))
    sig <- sig[grep("w",names(sig))]
    sigs80[[g]] <- sig

    sig <- which(sign(sry[,"25%"]) == sign(sry[,"75%"]))
    sig <- sig[grep("w",names(sig))]
    sigs50[[g]] <- sig

    # Plot posteriors and trace for sig50 vars
    pars <- names(sig)
    # Always include intercept
    if(!any(pars == "w0"))
        pars <- c("w0",pars)
    png(file=file.path(PLOTDIR, sprintf("%s-posterior.png",g)))
    print(plot(fit, pars=pars))
    dev.off()
    png(file=file.path(PLOTDIR, sprintf("%s-trace.png",g)))
    print(traceplot(fit, pars=pars))
    dev.off()
}

# Make summary table
num.sigs <- cbind(sapply(sigs95,length), sapply(sigs80,length), sapply(sigs50,length))
colnames(num.sigs) <- c("95%","80%","50%")
# Plot count polygons for different posterior intervals
counts <- melt(num.sigs, varnames=c("gene","postInterval"))
colnames(counts)[3] <- "n_significants"
png(file=file.path(PLOTDIR, "signif_counts.png"))
g <- ggplot(counts, aes(n_significants)) + geom_bar(aes(fill=postInterval)) + facet_grid(postInterval ~ .)
print(g)
dev.off()
num.sigs <- cbind(postInterval=colnames(num.sigs), as.data.frame(t(apply(num.sigs, 2, summary))))
#num.sigs <- data.frame(sigLevel=c("95%","80%","50%"),rbind(summary(sapply(sigs95,length)),summary(sapply(sigs80,length)),summary(sapply(sigs50,length))))
# Save output
save(num.sigs, sigs95, sigs80, sigs50, file=OUTFILE)
if(exists("param1"))
    table.out <- num.sigs

# Cleanup
rm(OUTFILE,RESULTDIR,PLOTDIR,sig,g,f,i,pars,counts,files,genes)
