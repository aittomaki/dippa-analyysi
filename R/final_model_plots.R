
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
theme_set(theme_bw())

d       <- table1
d.coefs <- table2
PLOTDIR <- document.dir



## Plot chosen var num vs R2 and num signif
#dm <- subset(d, select=c(gene, n_miRNAs, n_significant_miRNAs, R2_full, R2_full_adj))
dm <- subset(d, select=c(gene, n_miRNAs, n_significant_miRNAs, R2_full))
dm$delta_R2adj <- d$R2adj_full - d$R2adj_gene
dm <- melt(dm, id.vars=c("gene","n_miRNAs"))
#dm$adjusted[dm$variable == "R2_full"] <- "no"
#dm$adjusted[dm$variable == "R2_full_adj"] <- "yes"
#dm$adjusted[dm$variable == "n_significant_miRNAs"] <- "NA"
#dm$variable[dm$variable == "R2_full_adj"] <- "R2_full"
# Pretty the variable names
dm$variable <- revalue(dm$variable, c(
    n_significant_miRNAs = "N~significant~miRNAs",
    R2_full = "R[full]^2",
    delta_R2adj = "Delta~bar(R)^2"
))
g <- ggplot(dm, aes(x = n_miRNAs, y = value))#, color=adjusted, size=adjusted))
g <- g + geom_point(alpha=0.6) + geom_line(stat="smooth", method="loess", alpha=0.3)
g <- g + geom_ribbon(stat="smooth", method="loess", alpha=0.05)##, aes(color=NULL, group=adjusted))
g <- g + facet_grid(variable ~ ., scales = "free", switch="y", labeller=label_parsed)
#g <- g + scale_size_manual(values=c(1,1,0.5), guide=F)
#g <- g + scale_color_discrete(name="R2 adjusted", breaks=c("no","yes"))
#g <- g + ggtitle(parse(text=sub("U(.*)_a(.*)", "alpha:\\2~~~gamma:\\1", thresh)))
g <- g + labs(x="N miRNAs", y=NULL)
g <- g + theme(strip.text.y = element_text(size = 12))
plot.file <- file.path(PLOTDIR, "n_miRNAs_R2s.pdf")
ggsave(plot.file, g, height=7, width=9, dpi=600)




## Plot magnitude of miRNA versus IQR and SD
# First select only miRNA coefs
d2 <- subset(d.coefs, grepl("miR", variable), select=c(gene, variable, median, IQR, sd, significant))
d2m <- melt(d2, id.vars=c("gene","variable","median","significant"), variable.name="measure", value.name="value")
g <- ggplot(d2m, aes( x = median, y = value, color = significant))
g <- g + geom_point()
g <- g + facet_grid(measure ~ ., scales = "free", switch="y")
g <- g + labs(x="median of miRNA coefficient", y=NULL)
g <- g + theme(strip.text.y = element_text(size = 12))
plot.file <- file.path(PLOTDIR, "miRNA_magnitude_IQR.pdf")
ggsave(plot.file, g, height=7, width=9, dpi=600)

# Boxplot instead of scatter
g <- ggplot(d2m, aes( x = significant, y = value, fill = significant))
g <- g + geom_boxplot()
g <- g + facet_grid(measure ~ ., scales = "free", switch="y")
g <- g + labs(x="median of miRNA coefficient", y=NULL)
g <- g + theme(strip.text.y = element_text(size = 12))
plot.file <- file.path(PLOTDIR, "miRNA_magnitude_boxplot.pdf")
ggsave(plot.file, g, height=7, width=5, dpi=600)



## Plot miRNA median vs w0 median and miRNA median vs gene median
d3 <- d.coefs %>% select(gene, variable, median, significant) %>% group_by(gene) %>% mutate(median_gene=median[variable=="gene"], median_w0=median[variable=="w0"], ratio_gene=abs(median/median_gene), ratio_w0=abs(median/median_w0)) %>% ungroup %>% filter(grepl("miR", variable))
d3m <- melt(d3, id.vars=c("gene","variable","median","significant"), variable.name="measure", value.name="value")
g1 <- ggplot(d3, aes( x = significant, y = ratio_gene, fill = significant))
g1 <- g1 + geom_boxplot() + coord_cartesian(ylim = quantile(d3$ratio_gene, c(0.1, 0.9))*1.5)
g1 <- g1 + theme(legend.position = "none") + ylab("median(miRNA) / median(mRNA)") + xlab("miRNA significant")
g2 <- ggplot(d3, aes( x = significant, y = ratio_w0, fill = significant))
g2 <- g2 + geom_boxplot() + coord_cartesian(ylim = quantile(d3$ratio_w0, c(0.1, 0.9))*1.3)
g2 <- g2 + theme(legend.position = "none")

plot.file <- file.path(PLOTDIR, "miRNA_gene_w0_magnitude.pdf")
ggsave(plot.file, arrangeGrob(g1, nrow=1), height=5, width=6, dpi=400)


table.out <- d2
