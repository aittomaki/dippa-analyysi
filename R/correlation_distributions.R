# Script for computing correlations between epxression vars from different data

## INIT ####

library(reshape2)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
theme_set(theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()))

# Helper function for conversion from df to matrix
matrix.plz <- function(x) {
    m <- as.matrix(x[,-1])
    rownames(m) <- as.character(x[,1])
    m
}

# Input, params and output
if(exists("param1")) { # Anduril
    PLOTDIR    <- document.dir
    cormethod <- param1
    prot    <- table1
    gene    <- table2
    mirna   <- table3
    tarbase <- table4
} else { # non-Anduril
    PLOTDIR    <- "/home/viljami/tmp/plots"
    OUTFILE    <- "/home/viljami/wrk/tmp/plots/foo.csv"
    cormethod <- "pearson"
    WRKDIR <- Sys.getenv("WRKDIR")
    prot   <- read.delim(file.path(WRKDIR,"dippa-data","protein_normalized.csv"))
    gene   <- read.delim(file.path(WRKDIR,"dippa-data","gene_normalized.csv"))
    mirna  <- read.delim(file.path(WRKDIR,"dippa-data","mirna_normalized.csv"))
    tabase <- read.delim(file.path(WRKDIR,"dippa-data","tarbase_renamed.csv"))
    if(!dir.exists(PLOTDIR)) dir.create(PLOTDIR)
}

# Convert experssion data into matrices
samples <- as.character(prot[,1])
prot <- matrix.plz(prot)[samples,]
gene <- matrix.plz(gene)[samples,colnames(prot)]
mirna <- matrix.plz(mirna)[samples,]

# Form a list of validated pairs present in the data
# First clean mirna names (- and * to .)
tarbase$mirna <- gsub("-", ".", gsub("*", "-", tarbase$mirna, fixed=TRUE), fixed=TRUE)
mirna.names <- colnames(mirna)
tarbase <- filter(tarbase, geneName %in% colnames(gene) & mirna %in% mirna.names & positive_negative == "POSITIVE")
# Add variable indicating pair
tarbase <- mutate(tarbase, pair = paste(geneName, mirna, sep="-"))

# Plot params
h <- 4
w <- 5
dpi <- 300



## PROT VS GENE ####

# Correlation between proteins and genes
pg.cor <- melt(cor(prot, gene, method=cormethod))
names(pg.cor) <- c("protein","mrna","correlation")
# Add variable indicating of gene and prot match
pg.cor$group <- ifelse(pg.cor$protein == pg.cor$mrna, "matched", "unmatched")

# Plot
g <- ggplot(pg.cor, aes(x = correlation, group=group))
g <- g + geom_vline(xintercept = 0, color = "gray")
g <- g + geom_density(aes(color = group), size=1, trim=T)
g <- g + labs(x=sprintf("%s correlation", cormethod))
g <- g + scale_x_continuous(limit = c(-1,1))
g <- g + scale_color_brewer(name="protein-mRNA\npair", palette="Set1")
g <- g + guides(color = guide_legend(override.aes = list(fill = brewer.pal(3,"Set1")[1:2])))

plot.file <- file.path(PLOTDIR, "protein-gene-correlation.pdf")
ggsave(plot.file, g, height=h, width=w, dpi=dpi)






## MIRNA VS GENE ####
# Compute correlations
gm.cor <- melt(cor(gene, mirna, method=cormethod))
names(gm.cor) <- c("mrna","mirna","correlation")
gm.cor <- mutate(gm.cor, pair = paste(mrna, mirna, sep="-"))

# Make new df with different groups of mRNA-miRNA pairs
gm.data <- filter(gm.cor, pair %in% tarbase$pair)
gm.data$group <- "validated"
gm.smp <- sample_n(gm.cor, 3000, replace=T)
gm.smp$group <- "random"
gm.data <- rbind(gm.data, gm.smp)

# Plot
g <- ggplot(gm.data, aes(x = correlation, group = group))
g <- g + geom_vline(xintercept = 0, color = "gray", alpha=0.5)
g <- g + geom_density(aes(color = group), size=1, trim=T)
g <- g + labs(x=sprintf("%s correlation", cormethod))
g <- g + scale_x_continuous(limit = c(-1,1))
g <- g + scale_color_brewer(name="mRNA-miRNA\npair", palette="Set1")
g <- g + guides(color = guide_legend(override.aes = list(fill = brewer.pal(3,"Set1")[1:2])))

plot.file <- file.path(PLOTDIR, "gene-mirna-correlation.pdf")
ggsave(plot.file, g, height=h, width=w, dpi=dpi)





## MIRNA VS PROTEIN ####
# Compute correlations
pm.cor <- melt(cor(prot, mirna, method=cormethod))
names(pm.cor) <- c("protein","mirna","correlation")
pm.cor <- mutate(pm.cor, pair = paste(protein, mirna, sep="-"))

# Make new df with different groups of mRNA-miRNA pairs
pm.data <- filter(pm.cor, pair %in% tarbase$pair)
pm.data$group <- "validated"
pm.smp <- sample_n(pm.cor, 3000, replace=T)
pm.smp$group <- "random"
pm.data <- rbind(pm.data, pm.smp)

# Plot
g <- ggplot(pm.data, aes(x = correlation, group = group))
g <- g + geom_vline(xintercept = 0, color = "gray", alpha=0.5)
g <- g + geom_density(aes(color = group), size=1, trim=T)
g <- g + labs(x=sprintf("%s correlation", cormethod))
g <- g + scale_x_continuous(limit = c(-1,1))
g <- g + scale_color_brewer(name="protein-miRNA\npair", palette="Set1")
g <- g + guides(color = guide_legend(override.aes = list(fill = brewer.pal(3,"Set1")[1:2])))


plot.file <- file.path(PLOTDIR, "protein-mirna-correlation.pdf")
ggsave(plot.file, g, height=h, width=w, dpi=dpi)




## OUTPUT ####

table.out <- gm.data
