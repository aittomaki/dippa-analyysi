
library(reshape2)
library(plyr)
library(ggplot2)
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))

PLOTDIR <- document.dir
ppvs <- table1
lasso <- table2

w <- 4
h <- 4
dpi <- 600


# Delta R2adj comparison
ppvs$deltaR2 <- ppvs$R2adj_full - ppvs$R2adj_gene
lasso$deltaR2 <- lasso$R2adj_lasso - lasso$R2adj_gene
ppmi <- subset(ppvs, !(is.na(n_miRNAs)) & n_miRNAs > 0)
lami <- subset(lasso, !(is.na(n_miRNAs)) & n_miRNAs > 0)
d <- rbind(data.frame(method="PPVS", deltaR2=ppmi$deltaR2, n_miRNAs=ppmi$n_miRNAs),
           data.frame(method="lasso", deltaR2=lami$deltaR2, n_miRNAs=lami$n_miRNAs))

## Combined boxplot version
# dm <- melt(d)
# dm$variable <- revalue(dm$variable, c(
#     n_miRNAs = "N~~miRNAs",
#     deltaR2 = "Delta~bar(R)^2"
# ))
# g <- ggplot(dm, aes(x=method, y=value, group=method))
# g <- g + geom_boxplot()
# g <- g + facet_grid(variable ~ ., scales="free_y", switch = "y", label=label_parsed)
# g <- g + ylab(NULL) + theme(strip.background = element_blank())
# plot.file <- file.path(PLOTDIR, "R2_comparison.pdf")
# ggsave(plot.file, g, height=h, width=w, dpi=dpi)


## Separate plots version
# R2
g <- ggplot(d, aes(x=method, y=deltaR2, group=method))
g <- g + geom_boxplot()
g <- g + ylab(expression(Delta~bar(R)^2))
plot.file <- file.path(PLOTDIR, "R2_comparison.pdf")
ggsave(plot.file, g, height=h, width=w, dpi=dpi)
# N miRNAs
g <- ggplot(d, aes(x=method, y=n_miRNAs, group=method))
g <- g + geom_boxplot()
g <- g + ylab("N miRNAs chosen")
plot.file <- file.path(PLOTDIR, "n_miRNA_comparison.pdf")
ggsave(plot.file, g, height=h, width=w, dpi=dpi)




table.out <- d
