

### PREPROCESS EXPR DATA FOR REGRESSION ####

# Data
d <- as.matrix(table1[,2:ncol(table1)])
rownames(d) <- table1[,1]
# Parameter
SIGMA <- as.numeric(param1)

# transpose
d2 <- t(d)
# transform data back to linear scale
d2 <- 2^d2
# scale data to [0]
mins   <- apply(d2,2,min)
ranges <- apply(d2,2,function(x)diff(range(x)))
ds <- scale(d2, center=mins, scale=ranges) + SIGMA
# back to log2 scale
table.out <- log2(ds)

# Output also mins and ranges
mins.df <- data.frame(ID=names(mins), min=mins)
write.table(mins.df, file=get.output(cf,'optOut1'), sep='\t', row.names=FALSE)
ranges.df <- data.frame(ID=names(ranges), range=ranges)
write.table(ranges.df, file=get.output(cf,'optOut2'), sep='\t', row.names=FALSE)
save(mins,ranges,table.out, file=get.output(cf,'optOut3'))
rm(optOut1, optOut2, optOut3)
