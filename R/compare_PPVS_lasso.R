
## INIT ####

library(dplyr)
library(Vennerable)
library(ggplot2)
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))

if(exists("param1")) { # Anduril
    ppvs    <- table1
    lasso   <- table2
    tarbase <- table3
    PLOTDIR <- document.dir
} else {
    ppvs  <- read.delim("~/wrk/dippa-execute/finalModels_a0_5_g0_2/optOut2")
    lasso <- read.delim("~/wrk/dippa-execute/lassoModels/optOut1")
    tarbase <- read.delim("~/wrk/dippa-execute/tarBaseRename/csv.csv")
    PLOTDIR <- "~/tmp/ppvs_vs_lasso/"
    if(!dir.exists(PLOTDIR)) dir.create(PLOTDIR)
}

# Rename badly named median column in ppvs
ppvs <- dplyr::rename(ppvs, postmed=median)
# Filter for only miRNA vars and add column for combined gene-miRNA name
ppvs <- ppvs %>% filter(grepl("hsa", variable)) %>% mutate(pair = paste(gene,variable,sep="-"))
lasso <- lasso %>% filter(grepl("hsa", variable)) %>% mutate(pair = paste(gene,variable,sep="-"))
# Filter PPVS for only significant coefs
ppvs.sig <- ppvs %>% filter(significant == "yes")

# Clean tarbase
tarbase <- tarbase %>% mutate(pair = paste(geneName, mirna, sep="-"))
# filter to positive results
tarpos <- tarbase %>% filter(positive_negative == "POSITIVE")


## VENN DIAGRAMS ####
w <- 5
h <- 5
dpi <- 600

# Make a Venn diagram of overlap
# for all PPVS & lasso
vd <- Venn(list(PPVS=ppvs$pair, Lasso=lasso$pair))
plot.file <- file.path(PLOTDIR, "Venn_PPVS_Lasso.pdf")
pdf(plot.file, width = w, height = h)
plot(vd, doWeights=F, show=list(Faces=F))
dev.off()
# for signif PPVS & lasso
vd.sig <- Venn(list(PPVS=ppvs.sig$pair, Lasso=lasso$pair))
plot.file <- file.path(PLOTDIR, "Venn_PPVS-sig_Lasso.pdf")
pdf(plot.file, width = w, height = h)
plot(vd.sig, doWeights=F, show=list(Faces=F))
dev.off()
# for PPVS & signif PPVS & lasso
vd <- Venn(list(PPVS=ppvs$pair, PPVS.sig=ppvs.sig$pair, Lasso=lasso$pair))
plot.file <- file.path(PLOTDIR, "Venn_PPVS_PPVS-sig_Lasso.pdf")
pdf(plot.file, width = w, height = h)
plot(vd, doWeights=F, show=list(Faces=F))
dev.off()



## SCATTER AND COR PLOT ####

# Join the results for matching predictions
both <- inner_join(ppvs, lasso, by="pair")
both.sig <- inner_join(ppvs.sig, lasso, by="pair")
coefcor <- cor(both$postmed, both$coef)
coefcor.sig <- cor(both.sig$postmed, both.sig$coef)
coefspr <- cor(both$postmed, both$coef, method="spearman")
coefspr.sig <- cor(both.sig$postmed, both.sig$coef, method="spearman")

# Plot scatters on coefs, compute correlation
g1 <- ggplot(both, aes(x = postmed, y = coef))
g1 <- g1 + geom_abline(intercept = 0, slope = 1, color="gray")
g1 <- g1 + geom_point(aes(shape=significant), size=2)
g1 <- g1 + scale_shape(guide="none") #scale_color_manual(values=c("black","blue"), guide="none")
g1 <- g1 + labs(x="PPVS coef median", y="Lasso coef", title=sprintf("cor=%.2f",coefcor))
g2 <- ggplot(both.sig, aes(x = postmed, y = coef))
g2 <- g2 + geom_abline(intercept = 0, slope = 1, color="gray")
g2 <- g2 + geom_point()
g2 <- g2 + labs(x="PPVS coef median", y="Lasso coef", title=sprintf("cor=%.2f",coefcor.sig))
#g12 <- plot_grid(g1, g2, align = "h", ncol = 2)

plot.file <- file.path(PLOTDIR, "PPVS_Lasso_coef_scatter.pdf")
ggsave(plot.file, g1)
plot.file <- file.path(PLOTDIR, "PPVS-sig_Lasso_coef_scatter.pdf")
ggsave(plot.file, g2)




## COMPARE TO TARBASE ####

# function for the results
comptar <- function(x, y, coef="coef") {
    # positives and negatives found
    found <- inner_join(x, y, by="pair")
    npos <- nrow(filter(found, positive_negative=="POSITIVE"))
    nneg <- nrow(filter(found, positive_negative=="NEGATIVE"))
    # sign disagreement for positives
    wneg <- sum(found[[coef]] < 0 & found$up_down == "UP" & found$positive_negative=="POSITIVE")
    wpos <- sum(found[[coef]] > 0 & found$up_down == "DOWN" & found$positive_negative=="POSITIVE")

    # positives missed
    # needs to be filtered down to mirnas in data!!!
    nmis <- nrow(anti_join(y, x, by="pair"))
    # breast found (total 7 in tarbase for my data)
    yb <- filter(y, tissue == "Breast Cancerous Tissues")
    nb <- nrow(inner_join(x, yb, by="pair"))

    c(npos, wpos, wneg, nneg, nb, nmis)
}
tarres <- data.frame(Measure = c("N found", "N wrong pos sign", "N wrong neg sign", "N negative found", "N breast found (out of 7)", "N positive missed"),
                     PPVS = comptar(ppvs, tarbase, "postmed"),
                     PPVS.sig = comptar(ppvs.sig, tarbase, "postmed"),
                     lasso = comptar(lasso, tarbase),
                     both = comptar(both, tarbase),
                     both.sig=comptar(both.sig, tarbase))


## OUTPUT ####

table.out <- tarres
