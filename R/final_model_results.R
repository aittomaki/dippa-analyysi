
# Script for processing results from final models

# Input, params and output
if(exists("param1")) { # Anduril
    RESULTDIR <- param1
    PLOTDIR <- document.dir
    thresh <- param2
} else { # non-Anduril
    RESULTDIR <- "/home/viljami/wrk/finalmodel_results"
    PLOTDIR <- "/home/viljami/wrk/finalmodel_results/plots"
    OUTFILE <- "/home/viljami/wrk/finalmodel_results/foo.csv"
    if(!dir.exists(PLOTDIR)) dir.create(PLOTDIR)
    thresh <- "U0.2_a0.9"
    table1 <- read.delim("/home/viljami/wrk/tritonfake/dippa-data/n_chosen_variables.csv")
}

library(rstanarm)
library(xtable)
library(ggplot2)
library(reshape2)
theme_set(theme_bw())

# Get result file names
resfiles <- list.files(RESULTDIR, pattern = "*.rda")
# Get gene names from filenames
genes <- sub("finalmodel-\\d+-(\\w+).rda", "\\1", resfiles)
# Sort into alphabetical order by gene
isort <- order(genes)
resfiles <- resfiles[isort]
genes <- genes[isort]

# Result data frame
d <- data.frame(gene=genes, R2_gene=0, R2_gene_adjusted=0, R2_fullmodel=NA, R2_fullmodel_adjusted=NA, full_model_found="", gene_significant="", n_miRNAs=0, chosen_miRNAs="", n_significant_miRNAs=0, significant_miRNAs="", significant_miRNA_inds="", stringsAsFactors=F)

for (i in 1:length(resfiles)) {
#for (i in 1:2) {
    f <- resfiles[i]
    gene <- genes[i]
    # Load CV results for current gene
    load(file.path(RESULTDIR,f))

    n_vars <- table1[match(gene, table1$gene), thresh]

    # Rename chosen miRNAs to include dashes and asterisks
    chosen.mirnas <- sub("\\.$", "\\*", chosen.mirnas)
    chosen.mirnas <- gsub("\\.", "-", chosen.mirnas)

    # Gene only model results
    d[i,"full_model_found"] <- ifelse(params$n_vars > 0, "yes", "no")
    d[i,"R2_gene"] <- r2.gene["R2"]
    d[i,"R2_gene_adjusted"] <- r2.gene["R2.adjusted"]

    # Full model results
    if(params$n_vars > 1) {

        # Correct the adjusted R2s, they were computed wrong
        r2.adj <- 1 - (1-r2)*(params$n-1)/(params$n-params$n_vars+1)

        d[i,"R2_fullmodel"] <- median(r2)
        d[i,"R2_fullmodel_adjusted"] <- median(r2.adj)
        d[i,"n_miRNAs"] <- params$n_vars - 1
        d[i,"chosen_miRNAs"] <- paste(chosen.mirnas, collapse=", ")

        # posterior summary for projected weights
        postsum <- t(apply(as.matrix(spath$w[[length(spath$chosen)]][spath$chosen,]), 1, quantile, probs=c(0.025,0.1,0.5,0.9,0.975)))
        rownames(postsum) <- c("w0","gene",chosen.mirnas)
        signifs <- sign(postsum[,"2.5%"]) == sign(postsum[,"97.5%"])
        d[i,"gene_significant"] <- ifelse(signifs["gene"], "yes", "no")
        d[i,"n_significant_miRNAs"] <- sum(signifs[3:length(signifs)])
        d[i,"significant_miRNAs"] <- paste(chosen.mirnas[which(signifs[3:length(signifs)])], collapse=", ")
        d[i,"significant_miRNA_inds"] <- paste(which(signifs[3:length(signifs)]), collapse=",")
    } else {
        postsum <- posterior_interval(fit.gene, prob = 0.95, pars = "G")
        signifs <- sign(postsum[1,1]) == sign(postsum[1,2])
        d[i,"gene_significant"] <- ifelse(signifs, "yes", "no")
    }
}

# Output
table.out <- d
if(!exists("param1")) { # non-Anduril
    write.table(table.out, file=OUTFILE, sep="\t", row.names=F)
}

# Make a plot of var num vs R2 and num signif
dm <- subset(d, select=c(gene, n_miRNAs, n_significant_miRNAs, R2_fullmodel, R2_fullmodel_adjusted))
dm <- melt(dm, id.vars=c("gene","n_miRNAs"))
dm$adjusted[dm$variable == "R2_fullmodel"] <- "no"
dm$adjusted[dm$variable == "R2_fullmodel_adjusted"] <- "yes"
dm$adjusted[dm$variable == "n_significant_miRNAs"] <- "NA"
dm$variable[dm$variable == "R2_fullmodel_adjusted"] <- "R2_fullmodel"
g <- ggplot(dm, aes(x = n_miRNAs, y = value, color=adjusted))
g <- g + geom_point() + geom_smooth(method="lm")
g <- g + facet_grid(variable ~ ., scales = "free")
plot.file <- file.path(PLOTDIR, "n_miRNAs_R2s.png")
ggsave(plot.file, g, height=7, width=9, dpi=600)

# Make a latex table version
dltx <- d[,c("gene","R2_gene","R2_fullmodel","R2_gene_adjusted","R2_fullmodel_adjusted","full_model_found","gene_significant","n_miRNAs","n_significant_miRNAs","chosen_miRNAs")]
colnames(dltx) <- c("Gene","$R^2_{gene}$","$R^2_{full}$","$\\bar{R}^2_{gene}$","$\\bar{R}^2_{full}$","Model found","Gene signif","$N_{miRNA}$","$N_{signif miRNA}$","Chosen miRNAs")
ltable <- xtable(dltx, align="llrrrrccrrp{5cm}", digits=3, caption="Properties of final regression models.")
janitor <- function(x){paste0('{\\textbf{', x,'}}')}
document.out <- print(ltable, sanitize.colnames.function=janitor, print.results=F, include.rownames=F)
