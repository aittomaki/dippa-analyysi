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
