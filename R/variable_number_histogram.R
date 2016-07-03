
# Input, params and output
if(exists("param1")) { # Anduril
    varnums <- table1
    PLOTDIR <- document.dir
} else { # non-Anduril
    WRKDIR <- Sys.getenv("WRKDIR")
    DATADIR <- file.path(WRKDIR,"dippa-data")
    in.file <- file.path(DATADIR, "n_selected_variables")
    varnums <- read.delim(in.file)
    PLOTDIR <- "/home/viljami/wrk/cvresults/plots"
    if(!dir.exists(PLOTDIR)) dir.create(PLOTDIR)
}

# Required packages
library(dplyr)
library(reshape2)
library(ggplot2)
theme_set(theme_bw())

# Convert number of covariates to number of miRNA vars
varnums[,2:ncol(varnums)] <- varnums[,2:ncol(varnums)] - 1

# Make a histogram of num of vars selected for each threshold
d <- melt(varnums, id.vars="gene", value.name="N_variables", variable.name="threshold")
d$U <- sub("U(.*)_a(.*)", "\\1", d$threshold)
d$a <- sub("U(.*)_a(.*)", "\\2", d$threshold)
# Add num of NAs and zeros as a text label
numNAs <- d %>% group_by(U,a) %>% summarise(numNAs = sum(is.na(N_variables)), num0s = sum(N_variables==0, na.rm=T))
numNAs$NAlabs <- paste("#NA = ", numNAs$numNAs, sep="")
numNAs$zerolabs <- paste("#0    = ", numNAs$num0s, sep="")
numNAs <- data.frame(x=141, y=56, y0=51, numNAs)
# Construct plot
g <- ggplot(d, aes(x=N_variables))
g <- g + geom_histogram(binwidth = 5)
g <- g + facet_grid(U ~ a, labeller = label_bquote(cols=alpha: .(a), rows=U_f: .(U)))
g <- g + geom_text(aes(x, y, label=NAlabs, group=NULL), data=numNAs, size=4, color="grey40", vjust="top", hjust="left")
g <- g + geom_text(aes(x, y0, label=zerolabs, group=NULL), data=numNAs, size=4, color="grey40", vjust="top", hjust="left")
g <- g + xlab("N miRNA variables")
# Save plot
plot.file <- file.path(PLOTDIR, "ZZ_variable_number_hist.pdf")
ggsave(plot.file, g, height=7, width=9, dpi=600)

# Output (needed for Anduril)
table.out <- varnums
