
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
d <- rbind(data.frame(model="PPVS", deltaR2=ppvs$deltaR2), data.frame(model="lasso", deltaR2=lasso$deltaR2))
d <- cbind(face="Model fit", d)
# N miRNA comparison
ppmi <- subset(ppvs, !(is.na(n_miRNAs)) & n_miRNAs > 0)
lami <- subset(lasso, !(is.na(n_miRNAs)) & n_miRNAs > 0)
d2 <- rbind(data.frame(model="PPVS", n_miRNAs=ppvs$n_miRNAs), data.frame(model="lasso", n_miRNAs=lasso$n_miRNAs))
d2 <- cbind(face="Model size", d2)

d <- rbind(d,d2)

g <- ggplot(d, aes(x=model, y=deltaR2, group=model))
g <- g + geom_boxplot()
g <- g + ylab(expression(Delta~bar(R)^2))
plot.file <- file.path(PLOTDIR, "R2_comparison.pdf")
ggsave(plot.file, g, height=h, width=w, dpi=dpi)



g <- ggplot(d2, aes(x=model, y=n_miRNAs, group=model))
g <- g + geom_boxplot()
g <- g + ylab("N miRNAs chosen")
plot.file <- file.path(PLOTDIR, "n_miRNA_comparison.pdf")
ggsave(plot.file, g, height=h, width=w, dpi=dpi)


table.out <- d
