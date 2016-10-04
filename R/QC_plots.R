# Script for QC plots


## INIT ####

library(ggplot2)
library(ggbiplot)
library(reshape2)
library(dplyr)
library(cluster)
library(dendextend)
theme_set(theme_classic())

d      <- table1 # arrays/samples should be cols, probes/vars rows
groups <- table2
data.name <- param1 # e.g. "protein"
PLOTDIR   <- document.dir

id.name <- names(d)[1]
names(groups)[1] <- "Sample"
group.name <- names(groups)[2]

clust.method = "ward"

# make expression matrix
expr <- t(as.matrix(d[,-1]))
colnames(expr) <- d[,1]
rownames(expr) <- colnames(d)[-1]
z <- scale(expr)

# compute variable-wise median order of vars
med.order <- d[order(apply(expr, 2, median)),1]


## TRANSFORM DATA ####

# Make long version of data for ggplot (fuck this shit...)
dm <- melt(d, id.vars=id.name, variable.name="Sample")
# Add group var
dm <- left_join(dm, groups, by="Sample")
# Make sure variables come in alphabetical order
uvals <- sort(unique(dm[,id.name]))
#dm[,id.name] <- factor(dm[,id.name,], levels = uvals)





## BOXPLOTS ####

osize <- 0.3
h   <- 4
w   <- 10
dpi <- 800
bs  <- 12
legpos <- "none"
# if(identical(data.name, "protein")) {
#     legpos <- "top"
#     h <- h+1
# }

# Boxplot of arrays
g <- ggplot(dm, aes_string(x="Sample", y="value", fill=group.name))
g <- g + geom_boxplot(outlier.size = osize) #, coef = 100)
g <- g + scale_fill_brewer(palette="Set1")
g <- g + theme_classic(base_size = bs)
g <- g + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
g <- g + theme(legend.direction = "horizontal",
               legend.position = legpos)
g <- g + labs(y="Log2 expression", x=NULL)#,x="Tumor sample", title=sprintf("Boxplot of %s arrays",data.name))

plot.file <- file.path(PLOTDIR, sprintf("%s-sample_boxplot.png",instance.name))
ggsave(plot.file, g, height=h, width=w, dpi=dpi)

# Boxplot of variables
g <- ggplot(dm, aes_string(x=id.name, y="value"))
g <- g + geom_boxplot(outlier.size = osize) #, coef = 100)
g <- g + theme_classic(base_size = bs)
g <- g + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
g <- g + labs(y="Log2 expression", x=NULL)#, x=sprintf("%s variable",data.name))

plot.file <- file.path(PLOTDIR, sprintf("%s-variable_boxplot.png",instance.name))
ggsave(plot.file, g, height=h, width=w, dpi=dpi)

# Sorted version of variable boxplot
dm[,id.name] <- factor(dm[,id.name,], levels = med.order)
g <- ggplot(dm, aes_string(x=id.name, y="value"))
g <- g + geom_boxplot(outlier.size = osize) #, coef = 100)
g <- g + theme_classic(base_size = bs)
g <- g + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
g <- g + labs(y="Log2 expression", x=NULL)#, x=sprintf("%s variable (sorted by median expression)",data.name))

plot.file <- file.path(PLOTDIR, sprintf("%s-variable_boxplot_sorted.png",instance.name))
ggsave(plot.file, g, height=h, width=w, dpi=dpi)

# Boxplot of scaled variables
em <- melt(z, varnames=c("Sample",data.name))
g <- ggplot(em, aes_string(x=data.name, y="value"))
g <- g + geom_boxplot(outlier.size = osize) #, coef = 100)
g <- g + theme_classic(base_size = bs)
g <- g + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
g <- g + labs(y="Scaled log2 expression", x=NULL)#, x=sprintf("%s variable",data.name))

plot.file <- file.path(PLOTDIR, sprintf("%s-variable_boxplot_scaled.png",instance.name))
ggsave(plot.file, g, height=h, width=w, dpi=dpi)

# Sorted version of scaled variable boxplot
em[,data.name] <- factor(em[,data.name], levels = med.order)
g <- ggplot(em, aes_string(x=data.name, y="value"))
g <- g + geom_boxplot(outlier.size = osize) #, coef = 100)
g <- g + theme_classic(base_size = bs)
g <- g + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
g <- g + labs(y="Scaled log2 expression", x=NULL)#, x=sprintf("%s variable (sorted by unscaled median expression)",data.name))

plot.file <- file.path(PLOTDIR, sprintf("%s-variable_boxplot_scaled_sorted.png",instance.name))
ggsave(plot.file, g, height=h, width=w, dpi=dpi)






## DENSITY PLOTS ####

# Density plot of arrays
g <- ggplot(dm, aes_string(x="value", group="Sample", color="Sample")) +
    #geom_density() +
    geom_line(stat="density", alpha=0.5) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(legend.position = 'none') +
    labs(x=sprintf("Log2 %s expression",data.name))

plot.file <- file.path(PLOTDIR, sprintf("%s-sample_density.pdf",instance.name))
ggsave(plot.file, g, height=5, width=5, dpi=600)

# Density plot of variables
g <- ggplot(dm, aes_string(x="value", group=id.name, color=id.name)) +
    #geom_density() +
    geom_line(stat="density", alpha=0.5) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(legend.position = 'none') +
    labs(x=sprintf("Log2 %s expression",data.name))

plot.file <- file.path(PLOTDIR, sprintf("%s-variable_density.pdf",instance.name))
ggsave(plot.file, g, height=5, width=5, dpi=600)




## PCA PLOT ####

h <- 5
w <- 5
if(identical(data.name, "protein")) {
    legpos <- "left"
    w <- w+2
}

# copmute PCA
pca <- prcomp(expr, center=F, scale=F)
grp <- groups[match(colnames(d)[-1], groups[,1]),2]
g <- ggbiplot(pca, groups=grp, var.axes=F, scale=0)
g <- g + scale_color_brewer(name="Hospital", palette="Set1")
g <- g + theme_bw(base_size = 18)
g <- g + theme(legend.position = legpos)
g <- g + theme(axis.line.x = element_line(colour = "black"),
               axis.line.y = element_line(colour = "black"))
g <- g + labs(x="PCA1", y="PCA2")#, title=sprintf("%s",data.name))

plot.file <- file.path(PLOTDIR, sprintf("%s-pcaplot.pdf",instance.name))
ggsave(plot.file, g, height=h, width=w, dpi=300)





## CLUSTERING ####
clrs <- c(Radiumhospitalet="#E41A1C", Ullevaal="#377EB8")
dg <- as.dendrogram(agnes(expr), method=clust.method)
labels_colors(dg) <- clrs[grp][order.dendrogram(dg)]
labels(dg) <- "•"

plot.file <- file.path(PLOTDIR, sprintf("%s-clusterplot.pdf",instance.name))
pdf(file=plot.file, width=10, height=5)
plot(dg)
dev.off()





## OUTPUT ####

table.out <- dm
