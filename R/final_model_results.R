# Script for processing results from final model simulations

library(rstanarm)
library(xtable)
library(ggplot2)
library(reshape2)
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
    cv.chosen.vars <- table1
    prot           <- table2
    gene           <- table3
    mirna          <- table4
} else { # non-Anduril
    RESULTDIR <- "/home/viljami/wrk/finalmodel_results"
    PLOTDIR   <- "/home/viljami/wrk/finalmodel_results/plots"
    OUTFILE   <- "/home/viljami/wrk/finalmodel_results/foo.csv"
    WRKDIR    <- Sys.getenv("WRKDIR")
    cv.chosen.vars <- read.delim(file.path(WRKDIR,"dippa-data","n_chosen_variables.csv"))
    prot    <- read.delim(file.path(WRKDIR,"dippa-data","protein_normalized.csv"))
    gene    <- read.delim(file.path(WRKDIR,"dippa-data","gene_normalized.csv"))
    mirna   <- read.delim(file.path(WRKDIR,"dippa-data","mirna_normalized.csv"))
    if(!dir.exists(PLOTDIR)) dir.create(PLOTDIR)
    thresh <- "U0.2_a0.9"
}
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

# Result data frame
d <- data.frame(gene=genes, R2_gene=0, R2_gene_adj=0, R2_full=NA, R2_full_adj=NA, full_model_found="", gene_significant="", n_miRNAs=0, chosen_miRNAs="", n_significant_miRNAs=0, significant_miRNAs="", significant_miRNA_inds="", stringsAsFactors=F)

for (i in 1:length(resfiles)) {
#for (i in 1:2) {
    f <- resfiles[i]
    g <- genes[i]

    # Load CV results for current gene
    load(file.path(RESULTDIR,f))

    # Get number of covariates chosen in CV (includes gene!)
    n_vars <- cv.chosen.vars[match(g, cv.chosen.vars$gene), thresh]

    # Gene only model results
    d[i,"full_model_found"] <- ifelse(!is.na(n_vars), "yes", "no")
    d[i,"R2_gene"] <- r2.gene["R2"]
    d[i,"R2_gene_adj"] <- r2.gene["R2.adjusted"]

    if(!is.na(n_vars) && n_vars > params$n_vars) {
        warning(sprintf("%s : n_vars too big (%d > %d)!",g, n_vars, params$n_vars))
        n_vars <- 0
    }

    # Full model results
    if(!is.na(n_vars) && n_vars > 1) {
        # Get names of chosen mirnas
        mirna_vars <- spath$chosen[3:(n_vars+1)]-2 #1 is constant, 2 is gene
        mirna_names <- colnames(mirna)[mirna_vars]
        mirna_names <- gsub("\\.", "-", sub("\\.$", "\\*", mirna_names))

        # Create regression matrices
        y <- prot[,g]
        x <- cbind(1, gene[,g], mirna[,mirna_vars])
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

        # Rename chosen miRNAs to include dashes and asterisks
        chosen.mirnas <- gsub("\\.", "-", sub("\\.$", "\\*", chosen.mirnas))
        # Correct the adjusted R2s, they were computed wrong
        r2.adj <- 1 - (1-r2)*(params$n-1)/(params$n-params$n_vars+1)

        d[i,"R2_full"] <- median(R2)
        d[i,"R2_full_adj"] <- median(R2.adj)
        d[i,"n_miRNAs"] <- n_vars - 1
        d[i,"chosen_miRNAs"] <- paste(mirna_names, collapse=", ")

        # posterior summary for projected weights
        postsum <- t(apply(w, 1, quantile, probs=c(0.025,0.1,0.5,0.9,0.975)))
        rownames(postsum) <- c("w0", "gene", mirna_names)
        signifs <- sign(postsum[,"2.5%"]) == sign(postsum[,"97.5%"])
        d[i,"gene_significant"] <- ifelse(signifs["gene"], "yes", "no")
        d[i,"n_significant_miRNAs"] <- sum(signifs[3:length(signifs)])
        d[i,"significant_miRNAs"] <- paste(grep("miR",rownames(postsum)[signifs]), collapse=", ")
        d[i,"significant_miRNA_inds"] <- paste(which(signifs[3:length(signifs)]), collapse=",")
    } else {
        postsum <- posterior_interval(fit.gene, prob = 0.95)
        signifs <- sign(postsum["G",1]) == sign(postsum["G",2])
        d[i,"gene_significant"] <- ifelse(signifs, "yes", "no")
    }
}

# Output
table.out <- d
if(!exists("param1")) { # non-Anduril
    write.table(table.out, file=OUTFILE, sep="\t", row.names=F)
}

# Make a plot of var num vs R2 and num signif
#dm <- subset(d, select=c(gene, n_miRNAs, n_significant_miRNAs, R2_full, R2_full_adj))
dm <- subset(d, select=c(gene, n_miRNAs, n_significant_miRNAs, R2_full))
dm$delta_R2_adj <- d$R2_full_adj - d$R2_gene_adj
dm <- melt(dm, id.vars=c("gene","n_miRNAs"))
#dm$adjusted[dm$variable == "R2_full"] <- "no"
#dm$adjusted[dm$variable == "R2_full_adj"] <- "yes"
#dm$adjusted[dm$variable == "n_significant_miRNAs"] <- "NA"
#dm$variable[dm$variable == "R2_full_adj"] <- "R2_full"
g <- ggplot(dm, aes(x = n_miRNAs, y = value))#, color=adjusted, size=adjusted))
g <- g + geom_point(alpha=0.6) + geom_line(stat="smooth", method="loess", alpha=0.3)
g <- g + geom_ribbon(stat="smooth", method="loess", alpha=0.05)##, aes(color=NULL, group=adjusted))
g <- g + facet_grid(variable ~ ., scales = "free", switch="y")
#g <- g + scale_size_manual(values=c(1,1,0.5), guide=F)
#g <- g + scale_color_discrete(name="R2 adjusted", breaks=c("no","yes"))
g <- g + ggtitle(parse(text=sub("U(.*)_a(.*)", "alpha:\\2~~~gamma:\\1", thresh)))
g <- g + labs(x="N miRNAs", y=NULL)
plot.file <- file.path(PLOTDIR, sprintf("n_miRNAs_R2s_%s.png", thresh))
ggsave(plot.file, g, height=7, width=9, dpi=600)

# Make a latex table version
dltx <- d[,c("gene","R2_gene","R2_full","R2_gene_adj","R2_full_adj","full_model_found","gene_significant","n_miRNAs","n_significant_miRNAs","chosen_miRNAs")]
colnames(dltx) <- c("Gene","$R^2_{gene}$","$R^2_{full}$","$\\bar{R}^2_{gene}$","$\\bar{R}^2_{full}$","Model found","Gene signif","$N_{miRNA}$","$N_{signif miRNA}$","Chosen miRNAs")
ltable <- xtable(dltx, align="llrrrrccrrp{5cm}", digits=3, caption="Properties of final regression models.")
janitor <- function(x){paste0('{\\textbf{', x,'}}')}
document.out <- print(ltable, sanitize.colnames.function=janitor, print.results=F, include.rownames=F)
