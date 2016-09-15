
## INIT ####
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
library(cowplot)
library(latex2exp)
#theme_set(theme_bw())

d       <- table1
d.coefs <- table2
PLOTDIR <- document.dir

miRs <- d.coefs %>% filter(grepl("hsa", variable))
mRNAs <- d.coefs %>% filter(variable == "gene")
mods <- d %>% filter(full_model_found == "yes" & n_miRNAs > 0)

pos_miRs <- miRs %>% group_by(gene) %>% summarise(n_pos=sum(median>0), frac_pos=sum(median>0)/n())
pos_miRs <- left_join(pos_miRs, select(d, gene, n_miRNAs, n_significant_miRNAs), by="gene")




## Plot chosen var num vs R2 and num signif ####

#dm <- subset(d, select=c(gene, n_miRNAs, n_significant_miRNAs, R2_full, R2_full_adj))
d1 <- subset(d, select=c(gene, n_miRNAs, n_significant_miRNAs, R2_full))
d1$delta_R2adj <- d$R2adj_full - d$R2adj_gene
dm <- melt(d1, id.vars=c("gene","n_miRNAs"))
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
g <- g + labs(x="N miRNAs", y=NULL) + background_grid(major = "xy", minor = "none")
g <- g + theme(strip.text.y = element_text(size = 12))
plot.file <- file.path(PLOTDIR, "n_miRNAs_R2s.pdf")
ggsave(plot.file, g, height=7, width=9, dpi=600)

g1 <- ggplot(d1, aes(x = n_miRNAs, y = n_significant_miRNAs))
g1 <- g1 + geom_point(alpha=0.6) + geom_line(stat="smooth", method="loess", alpha=0.3)
g1 <- g1 + geom_ribbon(stat="smooth", method="loess", alpha=0.05)
g2 <- ggplot(d1, aes(x = n_miRNAs, y = R2_full))
g2 <- g2 + geom_point(alpha=0.6) + geom_line(stat="smooth", method="loess", alpha=0.3)
g2 <- g2 + geom_ribbon(stat="smooth", method="loess", alpha=0.05)
g3 <- ggplot(d1, aes(x = n_miRNAs, y = delta_R2adj))
g3 <- g3 + geom_point(alpha=0.6) + geom_line(stat="smooth", method="loess", alpha=0.3)
g3 <- g3 + geom_ribbon(stat="smooth", method="loess", alpha=0.05)
plot.file <- file.path(PLOTDIR, "n_miRNAs_R2s_cow.pdf")
g123 <- plot_grid(g1, g2, g3, align = "v", ncol=1)
save_plot(plot.file, g123, ncol=1)



## Plot n chosen/signif miRNAs vs positive miRNA coefs #####
g1 <- ggplot(pos_miRs, aes(x = n_miRNAs, y = frac_pos))
g1 <- g1 + geom_hline(yintercept = 0.50, color="gray") + geom_point()
g1 <- g1 + labs(x="N miRNAs chosen", y=TeX("Fraction of w_{jk} > 0"))
plot.file <- file.path(PLOTDIR, "n_miRNAs_vs_frac_pos.pdf")
save_plot(plot.file, g1)

g1 <- ggplot(pos_miRs, aes(x = n_significant_miRNAs, y = frac_pos))
g1 <- g1 + geom_point()
g1 <- g1 + labs(x=NULL)
g2 <- ggplot(pos_miRs, aes(x = n_significant_miRNAs, y = n_pos))
g2 <- g2 + geom_point()
g12 <- plot_grid(g1, g2, align = "v", ncol = 1)



## Plot delta R2adj vs fraction significant miRNAs ####

d1$fraction_sign_miRNAs <- d1$n_significant_miRNAs / d1$n_miRNAs
g1 <- ggplot(d1, aes(x = delta_R2adj, y = fraction_sign_miRNAs))
g1 <- g1 + geom_point()
g1 <- g1 + labs(x=NULL, y="Fraction of significant miRNAs")
g2 <- ggplot(d1, aes(x = R2_full, y = fraction_sign_miRNAs))
g2 <- g2 + geom_point()
g2 <- g2 + labs(x=NULL, y=NULL)
g3 <- ggplot(d1, aes(x = delta_R2adj, y = n_significant_miRNAs))
g3 <- g3 + geom_point()
g3 <- g3 + labs(x=parse(text="Delta~bar(R)^2"), y="N significant miRNAs")
g4 <- ggplot(d1, aes(x = R2_full, y = n_significant_miRNAs))
g4 <- g4 + geom_point()
g4 <- g4 + labs(x=parse(text="R[full]^2"), y=NULL)
plot.file <- file.path(PLOTDIR, "DeltaR2adj_vs_frac_sign_miRNAs.pdf")
#ggsave(plot.file, g, height=7, width=9, dpi=600)
#ggsave(plot.file, arrangeGrob(g1, g2, g3, g4, ncol=2), height=7, width=9, dpi=400)
g1234 <- plot_grid(g1, g2, g3, g4, ncol=2, align="hv")
#save_plot(plot.file, g1234, ncol=2)
ggsave(plot.file, g1234, height=7, width=9)



## Plot magnitude of miRNA versus IQR and SD ####

# First select only miRNA coefs
d2 <- subset(d.coefs, grepl("hsa", variable), select=c(gene, variable, median, IQR, sd, significant))
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



## Plot miRNA median vs gene median ####

d3 <- d.coefs %>% select(gene, variable, median, significant) %>% group_by(gene) %>% mutate(median_gene=median[variable=="gene"], median_w0=median[variable=="w0"], ratio_gene=abs(median/median_gene), ratio_w0=abs(median/median_w0)) %>% ungroup %>% filter(grepl("hsa", variable))
d3m <- melt(d3, id.vars=c("gene","variable","median","significant"), variable.name="measure", value.name="value")
g1 <- ggplot(d3, aes( x = significant, y = ratio_gene, fill = significant))
g1 <- g1 + geom_boxplot() + coord_cartesian(ylim = quantile(d3$ratio_gene, c(0.1, 0.9))*1.5)
g1 <- g1 + theme(legend.position = "none") + ylab("median(miRNA) / median(mRNA)") + xlab("miRNA significant")
g2 <- ggplot(d3, aes( x = significant, y = ratio_w0, fill = significant))
g2 <- g2 + geom_boxplot() + coord_cartesian(ylim = quantile(d3$ratio_w0, c(0.1, 0.9))*1.3)
g2 <- g2 + theme(legend.position = "none")

# Dont use g2 since miRNA vs w0 isnt actually that interesting
plot.file <- file.path(PLOTDIR, "miRNA_gene_median.pdf")
ggsave(plot.file, arrangeGrob(g1, nrow=1), height=5, width=6, dpi=400)



## Compute summary measures ####
sry <- numeric()
sry["N model found"] <- sum(d$full_model_found == "yes")
sry["No miRNAs in found model"] <- sum(d$full_model_found == "yes" & d$n_miRNAs == 0)
sry["N model not found"] <- sum(d$full_model_found == "no")
sry["Gene model better"] <- sum(mods$R2adj_gene > mods$R2adj_full)
sry["N genes"] <- length(unique(d.coefs$gene))
sry["N miRNAs total"] <- nrow(miRs)
sry["Average N miRNAs"] <- mean(mods$n_miRNAs)
sry["Median N miRNAs"] <- median(mods$n_miRNAs)
sry["Positive miRNAs"] <- sum(miRs$median > 0)
sry["Negative miRNAs"] <- sum(miRs$median < 0)
sry["Fraction negative miRNAs"] <- sum(miRs$median < 0)/nrow(miRs)
sry["N significant miRNAs total"] <- sum(miRs$significant == "yes")
sry["Average N significant miRNAs"] <- mean(mods$n_significant_miRNAs)
sry["Median N significant miRNAs"] <- median(mods$n_significant_miRNAs)
sry["Average fraction of significant miRNAs"] <- mean(mods$n_significant_miRNAs/mods$n_miRNAs)
sry["SD of fraction of significant miRNAs"] <- sd(mods$n_significant_miRNAs/mods$n_miRNAs)
sry["Positive significant miRNAs"] <- sum(miRs$median > 0 & miRs$significant == "yes")
sry["Negative significant miRNAs"] <- sum(miRs$median < 0 & miRs$significant == "yes")
sry["Fraction negative significant miRNAs"] <- sum(miRs$median < 0 & miRs$significant == "yes")/sum(miRs$significant == "yes")
sry["N mRNA coefs"] <- nrow(mRNAs)
sry["Positive mRNA coefs"] <- sum(mRNAs$median > 0)
sry["Negative mRNA coefs"] <- sum(mRNAs$median < 0)
sry["Fraction positive mRNA coefs"] <- sry["Positive mRNA coefs"]/nrow(mRNAs)
sry["N significant mRNA coefs"] <- sum(mRNAs$significant == "yes")
sry["Positive significant mRNA coefs"] <- sum(mRNAs$median > 0 & mRNAs$significant == "yes")
sry["Negative significant mRNA coefs"] <- sum(mRNAs$median < 0 & mRNAs$significant == "yes")
sry["Fraction positive significant mRNA coefs"] <- sry["Positive significant mRNA coefs"]/sum(mRNAs$significant == "yes")

table.out <- data.frame(summary = names(sry), value = sry)
