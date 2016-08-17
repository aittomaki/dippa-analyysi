
library(plyr)
library(reshape2)
library(xtable)
library(reporttools)

# Input
numeric_cols <- unlist(strsplit(param1, ","))
nominal_cols <- unlist(strsplit(param2, ","))
numeric_caption <- param3
nominal_caption <- param4

# Params for tables
longtbl <- FALSE
fontsize <- "small"
numeric_stats <- c("n","median","min","max","na")
cap_pos <- "top"

# Transform Grade variable
table1[,"Grade"] <- c("I","II","III","IV")[table1[,"Grade"]]
# Clean some variable values
table1$Histology <- revalue(table1$Histology, c("DCIS"="Ductal CIS"))
table1$Histology <- revalue(table1$Histology, c("PapillaryCIS"="Papillary CIS"))
table1$Multifocality <- revalue(table1$Multifocality, c("SingleTumor"="Single tumor"))
table1$ER <- revalue(table1$ER, c("neg"="negative"))
table1$ER <- revalue(table1$ER, c("pos"="positive"))
table1$PR <- revalue(table1$PR, c("neg"="negative"))
table1$PR <- revalue(table1$PR, c("pos"="positive"))
table1$HER2 <- revalue(table1$HER2, c("neg"="negative"))
table1$HER2 <- revalue(table1$HER2, c("pos"="positive"))

# Create tables
nominal_table <- tableNominal(table1[,nominal_cols], longtable=longtbl, font.size=fontsize, cap=nominal_caption, caption.placement=cap_pos)
numeric_table <- tableContinuous(table1[,numeric_cols], stats=numeric_stats, longtable=longtbl, font.size=fontsize, cap=numeric_caption, caption.placement=cap_pos)
# Add latex tables into output report
document.out <- c(nominal_table, "\n\n", numeric_table)

table.out <- table1
