
## INIT ####

library(dplyr)
library(Vennerable)
library(ggplot2)
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))

if(exists("param1")) { # Anduril
    ppvs    <- table1
    lasso   <- table2
    tarbase    <- table3
    mirtarbase <- table4
    PLOTDIR <- document.dir
} else {
    ppvs  <- read.delim("~/wrk/dippa-execute/finalModels_a0_5_g0_2/optOut2")
    lasso <- read.delim("~/wrk/dippa-execute/lassoModels/optOut1")
    tarbase    <- read.delim("~/wrk/dippa-execute/tarBaseRename/csv.csv")
    mirtarbase <- read.delim("~/wrk/dippa-execute/mirTarBaseRename/csv.csv")
    PLOTDIR <- "~/tmp/ppvs_vs_lasso/"
    if(!dir.exists(PLOTDIR)) dir.create(PLOTDIR)
}

# Rename badly named median column in ppvs
ppvs <- dplyr::rename(ppvs, postmed=median)
# Filter for only miRNA vars and add column for combined gene-miRNA name
ppvs <- ppvs %>% filter(grepl("hsa", variable)) %>% mutate(pair = paste(gene,variable,sep="-"))
lasso <- lasso %>% filter(grepl("hsa", variable)) %>% mutate(pair = paste(gene,variable,sep="-"))
# Add a variable for size of model from which each prediction comes from
ppvs <- ppvs %>% group_by(gene) %>% mutate(model_size_PPVS=n_distinct(variable)) %>% ungroup()
lasso <- lasso %>% group_by(gene) %>% mutate(model_size_lasso=n_distinct(variable)) %>% ungroup()
# Filter PPVS for only significant coefs and lasso for same number of top coefs
ppvs.sig <- ppvs %>% filter(significant == "yes")
lasso.top <- lasso[order(abs(lasso$coef), decreasing = T), ][1:nrow(ppvs.sig),]

# Clean tarbases
tarbase <- tarbase %>% mutate(pair = paste(geneName, mirna, sep="-"))
mirtarbase <- mirtarbase %>% mutate(pair = paste(TargetGene, miRNA, sep="-"))
# Filter tarbase to positive results
tarpos <- tarbase %>% filter(positive_negative == "POSITIVE")
# Validated targets
valids <- union(tarpos$pair, mirtarbase$pair)


## VENN DIAGRAMS ####
w <- 6
h <- 6
dpi <- 1200

# Make a Venn diagram of overlap
# for all PPVS & lasso
vd <- Venn(list(PPVS=ppvs$pair, lasso=lasso$pair))
plot.file <- file.path(PLOTDIR, "Venn-PPVS-lasso.pdf")
pdf(plot.file, width = w, height = h)
plot(vd, show=list(Faces=F))
dev.off()
# for signif PPVS & lasso top
vd.sig <- Venn(list("PPVS signif."=ppvs.sig$pair, "lasso top"=lasso.top$pair))
vd.sig <- compute.Venn(vd.sig)
labs <- VennGetSetLabels(vd.sig)
labs$x <- c(labs$x[1]-0.4, labs$x[2]+0.4)
labs$y <- labs$y + 0.1
vd.sig <- VennSetSetLabels(vd.sig, labs)
plot.file <- file.path(PLOTDIR, "Venn-PPVS_sig-lasso_top.pdf")
pdf(plot.file, width = w, height = h)
plot(vd.sig, show=list(Faces=F))
dev.off()
# for PPVS & signif PPVS & lasso
vd <- Venn(list(lasso=lasso$pair, "PPVS signif."=ppvs.sig$pair, PPVS=ppvs$pair))
plot.file <- file.path(PLOTDIR, "Venn-PPVS-PPVS_sig-lasso.pdf")
#png(plot.file, width = 800, height = 800, res=150)
pdf(plot.file, width = w, height = h)
plot(vd, doWeights=F, show=list(Faces=F))
dev.off()



## SCATTER AND COR PLOT ####

# Join the results for matching predictions
both <- inner_join(ppvs, lasso, by="pair")
both$dSize <- both$model_size_PPVS - both$model_size_lasso
both$dCoef <- both$postmed - both$coef
both.sig <- inner_join(ppvs.sig, lasso.top, by="pair")
coefcor <- cor(both$postmed, both$coef)
coefcor.sig <- cor(both.sig$postmed, both.sig$coef)
coefspr <- cor(both$postmed, both$coef, method="spearman")
coefspr.sig <- cor(both.sig$postmed, both.sig$coef, method="spearman")

# Plot scatters on coefs, compute correlation
g1 <- ggplot(both, aes(x = postmed, y = coef)) +
    geom_abline(intercept = 0, slope = 1, color="gray") +
    #geom_hline(yintercept = 0, color="gray", size=0.3) +
    #geom_vline(xintercept = 0, color="gray", size=0.3) +
    geom_point(aes(shape=significant), size=2) +
    scale_shape("PPVS significant") + #scale_shape(guide="none") #scale_color_manual(values=c("black","blue"), guide="none")
    #scale_color_brewer("Diff. in model size", type="div", palette="RdBu") +
    #scale_color_gradient2() +
    #theme(panel.grid.major = element_blank()) +
    labs(x="PPVS coefficient posterior median", y="lasso coefficient", title=sprintf("cor=%.2f",coefcor))
g2 <- ggplot(both.sig, aes(x = postmed, y = coef)) +
    geom_abline(intercept = 0, slope = 1, color="gray") +
    geom_hline(yintercept = 0, color="gray") +
    geom_vline(xintercept = 0, color="gray") +
    geom_point() +
    theme(panel.grid.major = element_blank()) +
    labs(x="PPVS coefficient posterior median", y="lasso coefficient", title=sprintf("cor=%.2f",coefcor.sig))
#g12 <- plot_grid(g1, g2, align = "h", ncol = 2)
g3 <- ggplot(both, aes(x=model_size_PPVS, y=model_size_lasso)) +
    geom_abline(intercept = 0, slope = 1, color="gray") +
    geom_point(aes(shape=significant), size=2) +
    labs(x="PPVS model size", y="lasso model size")
g4 <- ggplot(both, aes(x=dSize, y=dCoef)) +
    geom_point(aes(shape=significant), size=2) +
    labs(x="model size diff", y="coef diff")


plot.file <- file.path(PLOTDIR, "Scatter-PPVS-lasso.pdf")
ggsave(plot.file, g1, height=5, width=6)
plot.file <- file.path(PLOTDIR, "Scatter-PPVS_sig-lasso_top.pdf")
ggsave(plot.file, g2)
plot.file <- file.path(PLOTDIR, "Scatter-model-sizes.pdf")
ggsave(plot.file, g3)
plot.file <- file.path(PLOTDIR, "Scatter-model-size-diff.pdf")
ggsave(plot.file, g4)


## COMPARE TO TARBASES ####
#
# # function for the results
# comptar <- function(x, y, coef="coef") {
#     # positives and negatives found
#     found <- inner_join(x, y, by="pair")
#     npos <- nrow(filter(found, positive_negative=="POSITIVE"))
#     nneg <- nrow(filter(found, positive_negative=="NEGATIVE"))
#     # sign disagreement for positives
#     wneg <- sum(found[[coef]] < 0 & found$up_down == "UP" & found$positive_negative=="POSITIVE")
#     wpos <- sum(found[[coef]] > 0 & found$up_down == "DOWN" & found$positive_negative=="POSITIVE")
#
#     # positives missed
#     # needs filtering down to mirnas in data!!!
#     nmis <- nrow(anti_join(y, x, by="pair"))
#     # breast found (total 7 in tarbase for my data)
#     yb <- filter(y, tissue == "Breast Cancerous Tissues")
#     nb <- nrow(inner_join(x, yb, by="pair"))
#
#     c(npos, wpos, wneg, nneg, nb, nmis)
# }
# tarres <- data.frame(Measure = c("N validated found", "N wrong pos sign", "N wrong neg sign", "N negative found", "N breast found (out of 7)", "N positive missed (WRONG)"),
#                      PPVS = comptar(ppvs, tarbase, "postmed"),
#                      PPVS.sig = comptar(ppvs.sig, tarbase, "postmed"),
#                      lasso = comptar(lasso, tarbase),
#                      lasso.top = comptar(lasso.top, tarbase),
#                      both = comptar(both, tarbase),
#                      both.sig = comptar(both.sig, tarbase))

compref <- function(pred, ref, pair.name="pair") {
    nfound <- sum(pred[[pair.name]] %in% ref)
    prvalid <- nfound/nrow(pred)
    c(nfound, prvalid)
}
tarres <- data.frame(Measure = c("N validated found", "Percent validated"),
                     PPVS = compref(ppvs, valids),
                     PPVS.sig = compref(ppvs.sig, valids),
                     lasso = compref(lasso, valids),
                     lasso.top = compref(lasso.top, valids))

## OUTPUT ####

table.out <- tarres
