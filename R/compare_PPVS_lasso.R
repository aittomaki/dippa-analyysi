
library(dplyr)
library(Vennerable)
library(ggplot2)
library(cowplot)

if(exists("param1")) { # Anduril
    ppvs  <- table1
    lasso <- table2
    PLOTDIR <- document.dir
} else {
    ppvs  <- read.delim("~/wrk/dippa-execute/finalModels_a0_9_g0_2/optOut2.csv")
    lasso <- read.delim("~/wrk/dippa-execute/lassoModels/optOut1")
    PLOTDIR <- "~/tmp/ppvs_vs_lasso/"
    if(!dir.exists(PLOTDIR)) dir.create(PLOTDIR)
}

# Rename badly named median column
ppvs <- dplyr::rename(ppvs, postmed=median)

# Filter for only miRNA vars
ppvs <- ppvs %>% filter(grepl("hsa", variable))
lasso <- lasso %>% filter(grepl("hsa", variable))

# Add column for combined gene-miRNA name
ppvs <- ppvs %>% mutate(pair = paste(gene,variable,sep="-"))
lasso <- lasso %>% mutate(pair = paste(gene,variable,sep="-"))

# Filter PPVS for only significant coefs
ppvs.sig <- ppvs %>% filter(significant == "yes")

# Make a Venn diagram of overlap
vd <- Venn(list(PPVS=ppvs$pair, Lasso=lasso$pair))
plot.file <- file.path(PLOTDIR, "Venn_PPVS_Lasso.pdf")
pdf(plot.file)
plot(vd)
dev.off()
vd.sig <- Venn(list(PPVS=ppvs.sig$pair, Lasso=lasso$pair))
plot.file <- file.path(PLOTDIR, "Venn_PPVS-sig_Lasso.pdf")
pdf(plot.file)
plot(vd.sig)
dev.off()


# Make joins and plot scatters on coefs, compute correlation
both <- inner_join(ppvs, lasso, by="pair")
both.sig <- inner_join(ppvs.sig, lasso, by="pair")
coefcor <- cor(both$postmed, both$coef)
coefcor.sig <- cor(both.sig$postmed, both.sig$coef)
coefspr <- cor(both$postmed, both$coef, method="spearman")
coefspr.sig <- cor(both.sig$postmed, both.sig$coef, method="spearman")

g1 <- ggplot(both, aes(x = postmed, y = coef))
g1 <- g1 + geom_abline(intercept = 0, slope = 1, color="gray")
g1 <- g1 + geom_point()
g1 <- g1 + labs(x="PPVS coef median", y="Lasso coef", title=sprintf("all, cor=%.2f",coefcor))
g2 <- ggplot(both.sig, aes(x = postmed, y = coef))
g2 <- g2 + geom_abline(intercept = 0, slope = 1, color="gray")
g2 <- g2 + geom_point()
g2 <- g2 + labs(x="PPVS coef median", y="Lasso coef", title=sprintf("signif, cor=%.2f",coefcor.sig))
g12 <- plot_grid(g1, g2, align = "h", ncol = 2)

plot.file <- file.path(PLOTDIR, "Venn_PPVS_Lasso_coef_scatter.pdf")
ggsave(plot.file, g12)

table.out <- both
