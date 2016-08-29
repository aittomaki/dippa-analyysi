# Script for processing results from final model simulations

library(rstanarm)
library(xtable)
library(ggplot2)
library(reshape2)
library(plyr)
theme_set(theme_bw())

# Helper function for conversion from df to matrix
matrix.plz <- function(x) {
    m <- as.matrix(x[,-1])
    rownames(m) <- as.character(x[,1])
    m
}

# Input, params and output
if(exists("param1")) { # Anduril
    RESULTDIR <- param1
    thresh    <- param2 # e.g. U0.2_a0.9
    PLOTDIR   <- document.dir
    LATEX.CSV <- get.output(cf,'optOut1')
    COEFS.CSV <- get.output(cf,'optOut2')
    cv.chosen.vars <- table1
    prot           <- table2
    gene           <- table3
    mirna          <- table4
    rm(optOut1)
    rm(optOut2)
} else { # non-Anduril
    RESULTDIR <- "/home/viljami/wrk/finalmodel_results"
    PLOTDIR   <- "/home/viljami/wrk/finalmodel_results/plots"
    OUTFILE   <- "/home/viljami/wrk/finalmodel_results/foo.csv"
    LATEX.CSV <- "/home/viljami/wrk/finalmodel_results/ltable.csv"
    COEFS.CSV <- "/home/viljami/wrk/finalmodel_results/coefs.csv"
    WRKDIR    <- Sys.getenv("WRKDIR")
    cv.chosen.vars <- read.delim(file.path(WRKDIR,"dippa-data","n_chosen_variables.csv"))
    prot    <- read.delim(file.path(WRKDIR,"dippa-data","protein_normalized.csv"))
    gene    <- read.delim(file.path(WRKDIR,"dippa-data","gene_normalized.csv"))
    mirna   <- read.delim(file.path(WRKDIR,"dippa-data","mirna_normalized.csv"))
    if(!dir.exists(PLOTDIR)) dir.create(PLOTDIR)
    thresh <- "U0.2_a0.5"
}
SIGDIGS <- 3
# Sanity checking
if(!(thresh %in% colnames(cv.chosen.vars)))
    stop(sprintf("Threshold '%s' not available!", thresh))

# Convert experssion data into matrices
samples <- as.character(prot[,1])
prot <- matrix.plz(prot)[samples,]
gene <- matrix.plz(gene)[samples,]
mirna <- matrix.plz(mirna)[samples,]

# Get result file names
resfiles <- list.files(RESULTDIR, pattern = "*.rda")
# Get gene names from filenames
genes <- sub("finalmodel-\\d+-(\\w+).rda", "\\1", resfiles)
# Sort into alphabetical order by gene
gord <- order(genes)
genes <- genes[gord]
resfiles <- resfiles[gord]
rm(gord)

# Result table data frame
d <- data.frame(gene=genes, R2_gene=0, R2adj_gene=0,
                R2_full=NA, R2adj_full=NA, full_model_found="",
                gene_only_significant="", gene_full_significant="",
                n_miRNAs=0, chosen_miRNAs="",
                n_significant_miRNAs=0, significant_miRNAs="", significant_miRNA_inds="", stringsAsFactors=F)
# Data frame for coefs for all models (compile as a list)
d.coefs <- list()

for (i in 1:length(resfiles)) {
#for (i in 1:2) {
    f <- resfiles[i]
    g <- genes[i]

    # Load CV results for current gene
    load(file.path(RESULTDIR,f))

    # Get number of covariates chosen in CV (includes gene!)
    n_vars <- cv.chosen.vars[match(g, cv.chosen.vars$gene), thresh]
    if(!is.na(n_vars) && n_vars > params$n_vars) {
        warning(sprintf("%s : n_vars too big (%d > %d)!",g, n_vars, params$n_vars))
        n_vars <- NA
    }
    d[i,"n_miRNAs"] <- max(n_vars - 1, 0)
    d[i,"full_model_found"] <- ifelse(!is.na(n_vars), "yes", "no")

    # Gene only model results
    d[i,"R2_gene"] <- r2.gene["R2"]
    d[i,"R2adj_gene"] <- r2.gene["R2.adjusted"]
    # Check if gene significant for gene only model
    postsum.gene <- posterior_interval(fit.gene, prob = 0.95)
    signifs <- sign(postsum.gene["G",1]) == sign(postsum.gene["G",2])
    d[i,"gene_only_significant"] <- ifelse(signifs, "yes", "no")

    # Full model results
    if(!is.na(n_vars)) {
        if(n_vars > 1) {
            # Get names of chosen mirnas
            mirna_vars <- spath$chosen[3:(n_vars+1)]-2 #1 is constant, 2 is gene
            mirna_names <- colnames(mirna)[mirna_vars]
            mirna_names <- gsub("\\.", "-", sub("\\.$", "\\*", mirna_names))
            d[i,"chosen_miRNAs"] <- paste(mirna_names, collapse=",")

            x <- cbind(1, gene[,g], mirna[,mirna_vars])
        } else {
            # no miRNAs chosen
            mirna_names <- c()
            x <- cbind(1, gene[,g])
        }

        # Projected w's weren't computed for w0+g-models!
        # So R2's can't be computed (instead NA).
        # If that is corrected, then remove this if (and else)!
        if(n_vars > 1) {
            # Create regression matrices
            x <- x[,1:(n_vars+1)] # drops gene if only constant selected
            y <- prot[,g]
            n <- length(y)
            # Get projected weights for chosen mirnas
            w <- as.matrix(spath$w[[n_vars+1]][spath$chosen[1:(n_vars+1)],])
            # Compute predicted prot expr
            ypred <- x %*% w
            # Compute R2 and adjusted R2
            R2 <- 1 - colSums((y-ypred)^2)/sum((y-mean(y))^2)
            R2.adj <- 1 - (1-R2)*(n-1)/(n-n_vars+1)
            resid <- y-ypred
            R2.var <- 1-apply(resid,2,var)/var(y)
            R2.var.adj <- 1 - (1-R2.var)*(n-1)/(n-n_vars+1)

            d[i,"R2_full"] <- median(R2)
            d[i,"R2adj_full"] <- median(R2.adj)

            # Posterior summary for projected weights
            postsum <- t(apply(w, 1, quantile, probs=c(0.025,0.1,0.25,0.5,0.75,0.9,0.975)))
            rownames(postsum) <- c("w0", "gene", mirna_names)[1:(n_vars+1)]
            sdev <- as.vector(apply(w, 1, sd))
            wgt <- as.vector(apply(w, 1, function(x) {ifelse(median(x)>0, sum(x>0)/length(x), sum(x<0)/length(x))}))
            signifs <- sign(postsum[,"2.5%"]) == sign(postsum[,"97.5%"])
            cfs <- data.frame(gene=g, variable=rownames(postsum), median=postsum[,"50%"], IQR=postsum[,"75%"]-postsum[,"25%"], sd=sdev, significant=ifelse(signifs, "yes", "no"), weight=wgt)
            d.coefs <- c(d.coefs, list(cfs))
            if(n_vars > 0) {
                d[i,"gene_full_significant"] <- ifelse(signifs["gene"], "yes", "no")
                if(n_vars > 1) {
                    d[i,"n_significant_miRNAs"] <- sum(signifs[3:length(signifs)])
                    d[i,"significant_miRNAs"] <- paste(grep("miR",rownames(postsum)[signifs], value=T), collapse=",")
                    d[i,"significant_miRNA_inds"] <- paste(which(signifs[3:length(signifs)]), collapse=",")
                }
            }
        } else {
            d[i,"R2_full"] <- NA
            d[i,"R2adj_full"] <- NA
        }
    }
}

# Output
table.out <- d
if(!exists("param1")) { # non-Anduril
    write.table(table.out, file=OUTFILE, sep="\t", row.names=F)
}
d.coefs <- do.call(rbind, d.coefs)
write.table(d.coefs, file=COEFS.CSV, sep="\t", row.names=F)


# Make a version of the table for latex in Anduril
numf <- paste0("%.",SIGDIGS,"f")
dlat <- d[,c("gene","R2_gene","R2_full","n_miRNAs","significant_miRNAs")]
# Combine n_miRNAs and n_significant_miRNAs
nonNA <- !is.na(dlat$n_miRNAs) & dlat$n_miRNAs > 0
dlat$n_miRNAs[nonNA] <- paste0(d$n_miRNAs, " (", d$n_significant_miRNAs, ")")[nonNA]
# Add gene significants as asterisks
gene_only_sigs <- revalue(d$gene_only_significant, c(yes="$^{\\ast}$", no=""))
gene_full_sigs <- revalue(d$gene_full_significant, c(yes="$^{\\ast}$", no=""))
dlat$R2_gene <- paste0(sprintf(numf, d$R2_gene), " (", sprintf(numf, d$R2adj_gene), ")", gene_only_sigs)
dlat$R2_full[nonNA] <- paste0(sprintf(numf, d$R2_full), " (", sprintf(numf, d$R2adj_full), ")", gene_full_sigs)[nonNA]
# Round R2_adjs
write.table(dlat, file=LATEX.CSV, sep="\t", row.names=F)

# Make a plot of var num vs R2 and num signif
# #dm <- subset(d, select=c(gene, n_miRNAs, n_significant_miRNAs, R2_full, R2_full_adj))
# dm <- subset(d, select=c(gene, n_miRNAs, n_significant_miRNAs, R2_full))
# dm$delta_R2adj <- d$R2adj_full - d$R2adj_gene
# dm <- melt(dm, id.vars=c("gene","n_miRNAs"))
# #dm$adjusted[dm$variable == "R2_full"] <- "no"
# #dm$adjusted[dm$variable == "R2_full_adj"] <- "yes"
# #dm$adjusted[dm$variable == "n_significant_miRNAs"] <- "NA"
# #dm$variable[dm$variable == "R2_full_adj"] <- "R2_full"
# # Pretty the variable names
# dm$variable <- revalue(dm$variable, c(
#     n_significant_miRNAs = "N~significant~miRNAs",
#     R2_full = "R[full]^2",
#     delta_R2adj = "Delta~bar(R)^2"
# ))
# g <- ggplot(dm, aes(x = n_miRNAs, y = value))#, color=adjusted, size=adjusted))
# g <- g + geom_point(alpha=0.6) + geom_line(stat="smooth", method="loess", alpha=0.3)
# g <- g + geom_ribbon(stat="smooth", method="loess", alpha=0.05)##, aes(color=NULL, group=adjusted))
# g <- g + facet_grid(variable ~ ., scales = "free", switch="y", labeller=label_parsed)
# #g <- g + scale_size_manual(values=c(1,1,0.5), guide=F)
# #g <- g + scale_color_discrete(name="R2 adjusted", breaks=c("no","yes"))
# #g <- g + ggtitle(parse(text=sub("U(.*)_a(.*)", "alpha:\\2~~~gamma:\\1", thresh)))
# g <- g + labs(x="N miRNAs", y=NULL)
# g <- g + theme(strip.text.y = element_text(size = 12))
# plot.file <- file.path(PLOTDIR, sprintf("n_miRNAs_R2s_%s.pdf", thresh))
# ggsave(plot.file, g, height=7, width=9, dpi=600)

# Make a latex table version in R (NOT USED)
dltx <- dlat
#colnames(dltx) <- c("Gene","$R^2_{gene}$","$R^2_{full}$","$\\bar{R}^2_{gene}$","$\\bar{R}^2_{full}$","$N_{miRNA}$","Chosen miRNAs")
#ltable <- xtable(dltx, align="llrrrrrp{5cm}", digits=SIGDIGS, caption="Properties of final regression models.")
colnames(dltx) <- c("Gene","$R^2_{gene}$","$R^2_{full}$","$N_{miRNA}$","Chosen miRNAs")
ltable <- xtable(dltx, align="llllrp{8cm}", digits=SIGDIGS, caption="Properties of final regression models.")
janitor <- function(x){paste0('{\\textbf{', x,'}}')}
document.out <- print(ltable, sanitize.colnames.function=janitor, print.results=F, include.rownames=F)
