## Getting all gene and protein IDs from Ensembl

library(biomaRt)

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

filts <- listFilters(mart, what=c("name","description","options","fullDescription","type","operation"))
atts <- listAttributes(mart)

chromosomes <- c(1:22, "X", "Y")

allgenes <- getBM(filters="chromosome_name", values=chromosomes, attributes=c("ensembl_gene_id","chromosome_name"), mart=mart)
protids <- getBM(filters="chromosome_name", values=chromosomes, attributes=c("ensembl_peptide_id","chromosome_name"), mart=mart)
protgenes <- getBM(filters="chromosome_name", values=chromosomes, attributes=c("ensembl_gene_id","ensembl_peptide_id","chromosome_name"), mart=mart)
#protgenes <- getBM(filters="with_protein_id", values=TRUE, attributes=c("ensembl_gene_id"), mart=mart)

codinggenes <- getBM(filters="ensembl_peptide_id", values=protids$ensembl_peptide_id, attributes=c("ensembl_gene_id"), mart=mart)

write.table(protids, file="~/tmp/ensembl_ids/ensembl_protein_ids.csv", row.names=F, quote=F)
write.table(codinggenes, file="~/tmp/ensembl_ids/ensembl_coding_genes.csv", row.names=F, quote=F)




#### THE FOLLOWING IS FOR GETTING INTERACTIONS FROM miRWalk FOR CODING GENES ##########

# Take a 200 gene sample from the coding genes
samplesize <- 1000
gene_sample <- data.frame(ensembl_gene_id=sample(codinggenes[,1], samplesize))
# Write sample out in 50 gene batches
for(i in 1:ceiling(samplesize/50)) {
    ind <- seq(i,samplesize,ceiling(samplesize/50))
    write.table(gene_sample[ind,], file=paste("~/tmp/ensembl_ids/gene_sample",i,".txt",sep=""), row.names=F, quote=F)
}

# Load interactions back in (downloaded manually from miRWalk website
# use dplyr for help
library(dplyr)
inter <- list()
for(i in 1:ceiling(samplesize/50))
    inter[[i]] <- read.delim(paste("~/tmp/ensembl_ids/genemirna",i,".csv",sep=""))
inter <- do.call(rbind, inter)
# Check number of genes with no interactions
n_no_inter <- samplesize - length(unique(inter$Gene))
# Compute number of interacting miRNAs per gene
mirnas_per_gene <- inter %>% group_by(Gene) %>% summarise(n_distinct(miRNA)) %>% select(2) %>% unlist(use.names=F)
# Add zeros
mirnas_per_gene <- c(mirnas_per_gene, rep(0,n_no_inter))
summary(mirnas_per_gene)
