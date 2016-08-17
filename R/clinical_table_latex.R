
library(plyr)
library(reshape2)
library(xtable)

numerical_cols <- unlist(strsplit(param1, ","))
nominal_cols <- unlist(strsplit(param2, ","))
numeric_cap <- param3
nominal_cap <- param4

fontsize <- "\\small"

table1 <- table1[,c(numerical_cols,nominal_cols)]

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


### Function for pretty clinical tables
prettyTable <- function(x, num_cols=NULL) {

    # Uses Plyr package
    require(plyr)
    require(reshape2)

    if(is.null(num_cols))
        num_cols <- names(x)[sapply(x, is.numeric)]

    ## Do num cols first
    d_num <- ldply(x[,num_cols,drop=F], function(x) c(median(x,na.rm=T),min(x,na.rm=T),max(x,na.rm=T),sum(is.na(x))))
    names(d_num) <- c("Variable","Median","Min","Max","Missing")
    d_num[,"Range"] <- paste(d_num[,"Min"], "--", d_num[,"Max"], sep=" ")
    d_num <- d_num[,c("Variable","Median","Range","Missing")]
    #d_num <- melt(d_num[,c("Variable","Median","Range","Missing")], id.vars="Variable", variable.name="Level")
    #d_num <- d_num[order(d_num$Variable),]
    #d_num[,1] <- ifelse(duplicated(d_num[,1]), "", d_num[,1])

    ## Other cols
    i <- !(names(x) %in% num_cols)
    # Transform to counts
    d <- ldply(x[,i,drop=F], function(x) t(rbind(names(table(x,useNA="ifany")),table(x,useNA="ifany"),paste0(round(prop.table(table(x,useNA="ifany"))*100,0),"%"))))

    # Spaces for variable names
    d[,1] <- ifelse(duplicated(d[,1]), "", d[,1])
    # missing instead of NA
    d[,2] <- as.character(d[,2])
    d[is.na(d[,2]),2] <- "missing"
    # Better column names
    names(d) <- c("Variable", "Level", "Count", "Fraction")

    # Output both tables
    list(d_num,d)
}

tables <- prettyTable(table1, num_cols=numerical_cols)

# Make a LaTeX fragment of numeric table
numtable <- xtable(tables[[1]], caption=numeric_cap, align="lllrr", digits=0)
numtable <- print(numtable, size=fontsize, caption.placement="top", print.results=F, include.rownames=F)
# Remove table ending
comb.ltable1 <- sub("   \\hline\n\\end{tabular}\n\\end{table}\n","", numtable, fixed=T)

# Make a LaTeX fragment of nominal table
nomtable <- xtable(tables[[2]], caption=nominal_cap, align="lllrr", digits=0)
nomtable <- print(nomtable, size=fontsize, caption.placement="top", print.results=F, include.rownames=F)
# Remove table beginning
comb.ltable2 <- sub("\\begin{table}[ht]\n\\centering\n\\begin{tabular}{lllrr}\n","", nomtable, fixed=T)

document.out <- c(numtable, "\n\n", nomtable)

table.out <- tables[[2]]
