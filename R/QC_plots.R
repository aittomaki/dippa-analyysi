# Script for QC plots


## INIT ####

library(ggplot2)
library(ggbiplot)
library(reshape2)
library(dplyr)
library(cluster)
library(dendextend)
theme_set(theme_classic())

df <- table1 # arrays/samples should be cols, probes/vars rows
groups <- table2
data.name <- param1 # e.g. "protein"
PLOTDIR   <- document.dir

id.name <- names(df)[1]
names(groups)[1] <- "Sample"
group.name <- names(groups)[2]

clust.method = "ward"



## TRANSFORM DATA ####

# Make long version of data for ggplot (fuck this shit...)
dm <- melt(df, id.vars=id.name, variable.name="Sample")
# Add group var
dm <- left_join(dm, groups, by="Sample")



## BOXPLOT ####

# Compute limits for y-axis
dmg <- dm %>% group_by(Sample) %>% summarise()

g <- ggplot(dm, aes(x=Sample, y=value))
g <- g + geom_boxplot(aes_string(fill=group.name), outlier.size = 0.1) #, coef = 100)
g <- g + theme_classic()
g <- g + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
g <- g + labs(y=NULL, title=sprintf("Boxplot of %s arrays",data.name))

plot.file <- file.path(PLOTDIR, sprintf("%s-boxplot.pdf",instance.name))
ggsave(plot.file, g, height=5, width=10, dpi=600)



## PCA PLOT ####

# make expression matrix
expr <- t(as.matrix(df[,-1]))
# copmute PCA
pca <- prcomp(expr, center=F, scale=F)
grp <- groups[match(colnames(df)[-1], groups[,1]),2]
g <- ggbiplot(pca, groups=grp, var.axes=F, scale=0)
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
g <- g + theme(axis.line.x = element_line(colour = "black"),
               axis.line.y = element_line(colour = "black"))
g <- g + labs(x="PCA1", y="PCA2", title=sprintf("%s",data.name))

plot.file <- file.path(PLOTDIR, sprintf("%s-pcaplot.pdf",instance.name))
ggsave(plot.file, g, height=5, width=5, dpi=600)



## CLUSTERING ####
clrs <- c(Radiumhospitalet="#F8766D",Ullevaal="#00BFC4")
dg <- as.dendrogram(agnes(expr), method=clust.method)
labels_colors(dg) <- clrs[grp][order.dendrogram(dg)]
labels(dg) <- "•"

plot.file <- file.path(PLOTDIR, sprintf("%s-clusterplot.pdf",instance.name))
pdf(file=plot.file, width=10, height=5)
plot(dg)
dev.off()


## OUTPUT ####

table.out <- dm
