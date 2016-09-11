# R script for computing Lasso versions of models

library(glmnet)

# Helper function for conversion from df to matrix
matrix.plz <- function(x) {
    m <- as.matrix(x[,-1])
    rownames(m) <- as.character(x[,1])
    m
}

# Input, params and output
if(exists("param1")) { # Anduril
    PLOTDIR <- document.dir
    prot    <- table1
    gene    <- table2
    mirna   <- table3
    nfolds  <- as.numeric(param1)
    lambda  <- param2
} else { # non-Anduril
    WRKDIR  <- Sys.getenv("WRKDIR")
    PLOTDIR <- file.path(WRKDIR,"dippa-execute","lassoModels")
    prot    <- read.delim(file.path(WRKDIR,"dippa-data","protein_normalized.csv"))
    gene    <- read.delim(file.path(WRKDIR,"dippa-data","gene_normalized.csv"))
    mirna   <- read.delim(file.path(WRKDIR,"dippa-data","mirna_normalized.csv"))
    if(!dir.exists(PLOTDIR)) dir.create(PLOTDIR)
    nfolds <- 10
    lambda <- "lambda.1se"
}

# Ensure matching samples and genes
samples <- as.character(prot[,1])
prot  <- matrix.plz(prot)[samples,]
genes <- colnames(prot)
gene  <- matrix.plz(gene)[samples, genes]
mirna <- matrix.plz(mirna)[samples,]

n <- length(samples)

# Result table data frame
d <- data.frame(gene=genes, R2_gene=0, R2adj_gene=0,
                R2_lasso=NA, R2adj_lasso=NA,
                gene_only_significant="", gene_lasso_included="",
                n_miRNAs=0, chosen_miRNAs="", stringsAsFactors=F)

# For each gene, build lasso model
for(i in 1:length(genes)) {
    g <- genes[i]

    # Build matrices for modeling
    y <- prot[,g]
    x <- cbind( gene[,g,drop=F], mirna )

    # Fit lasso and lm model
    cvfit <- cv.glmnet(x, y, type.measure="mse", nfolds=nfolds)

    # Get chosen vars and their coefs
    cfs <- as.matrix(coef(cvfit, s=lambda))
    chosen.vars <- rownames(cfs)[cfs[,1] != 0]
    # Drop intercept
    if("(Intercept)" %in% chosen.vars)
        chosen.vars <- chosen.vars[-1*match("(Intercept)", chosen.vars)]
    n_vars <- length(chosen.vars)
    chosen.mirnas <- grep("miR", chosen.vars)
    n_mirnas <- length(chosen.mirnas)
    gene.chosen <- ifelse(n_vars > n_mirnas, "yes", "no")

    # Get R2 and R2adj
    r2 <- cvfit$glmnet.fit$dev.ratio[cvfit[[lambda]] == cvfit$glmnet.fit$lambda]
    r2.adj <- 1 - (1-r2)*(n-1)/(n-n_vars-1)
    # Compute R2 myself and compare
    ypred <- predict(cvfit, newx=x, s=lambda)
    r2.comp <- 1 - sum((y-ypred)^2)/sum((y-mean(y))^2)
    if(r2-r2.comp > .Machine$double.neg.eps)
        warning(sprintf("GLMnet and computed R2 differ! (%f vs %f)",r2,r2.comp))


    # Convert . in miRNA names to - and *
    chosen.mirnas <- gsub("\\.", "-", sub("\\.$", "\\*", chosen.mirnas))

    # Gather results
    d[i, "chosen_mirnas"] <- paste(chosen.mirnas, collapse=",")
    d[i, "n_mirnas"] <- n_mirnas
    d[i, "R2_lasso"] <- r2
    d[i, "R2adj_lasso"] <- r2.adj
    d[i, "gene_lasso_included"] <- gene.chosen

    # Fit gene only model with simple lm
    dlm <- data.frame(P=y, G=gene[,g])
    fit <- lm(P ~ G, data=dlm)
    d[i,"R2_gene"] <- summary(fit)$r.squared
    d[i,"R2adj_gene"] <- summary(fit)$adj.r.squared
    d[i,"gene_only_significant"] <- ifelse(coef(summary(fit))["G","Pr(>|t|)"] < 0.05, "yes", "no")
}

table.out <- d
