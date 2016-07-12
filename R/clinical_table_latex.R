
library(plyr)
library(reshape2)
library(xtable)

numerical_cols <- unlist(strsplit(param1, ","))
nominal_cols <- unlist(strsplit(param2, ","))
caption <- param3

table1 <- table1[,c(numerical_cols,nominal_cols)]

# Transform Grade variable
table1[,"Grade"] <- c("I","II","III","IV")[table1[,"Grade"]]


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
    d[is.na(d[,2]),2] <- "NA"
    # Better column names
    names(d) <- c("Variable", "Level", "Count", "Fraction")

    # Output both tables
    list(d_num,d)
}

tables <- prettyTable(table1, num_cols=numerical_cols)

# Make a LaTeX fragment of first table piece
ltable1 <- xtable(tables[[1]], caption=caption, align="lllrr", digits=0)
ltable1 <- print(ltable1, caption.placement = "top", print.results=F, include.rownames=F)
# Remove table ending
ltable1 <- sub("   \\hline\n\\end{tabular}\n\\end{table}\n","", ltable1, fixed=T)

# Make a LaTeX fragment of second table piece
ltable2 <- xtable(tables[[2]], digits=0)
ltable2 <- print(ltable2, print.results=F, include.rownames=F)
# Remove table beginning
ltable2 <- sub("\\begin{table}[ht]\n\\centering\n\\begin{tabular}{llll}\n","", ltable2, fixed=T)

document.out <- c(ltable1, ltable2)

table.out <- tables[[2]]
