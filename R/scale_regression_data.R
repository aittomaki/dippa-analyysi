## Anduril R script for transposing and normalizing
## a data frame to [0,1] and then taking log2.
## Assumes input data is in log-scale (and converts
## and converts before normalizing). Plots the
## densities of all columns after each step.
## Also makes a density plot of z-normalized data.


# DATA
d <- as.matrix(table1[,2:ncol(table1)])
rownames(d) <- table1[,1]
# Parameter
SIGMA <- as.numeric(param1)
name <- param2
line.colors <- rainbow(ncol(d))



# TRANSPOSE (to have vars as columns instead of )
dt <- t(d)
print(paste("Variables to plot: ", as.character(ncol(dt)), sep=""))
print(paste("First var name: ", colnames(dt)[1], sep=""))

# Make density plot in log2 scale
png(file.path(document.dir,sprintf('%s-log2.png',instance.name)))
y.max <- max(apply(dt, 2, function(x) max(density(x)$y)))
x.min <- min(apply(dt, 2, function(x) min(density(x)$x)))
x.max <- max(apply(dt, 2, function(x) max(density(x)$x)))
plot(0,0, ylim=c(0,y.max), xlim=c(x.min,x.max), type='n')
for(i in 1:ncol(dt)) {
	lines(density(dt[,i]), col=line.colors[i])
}
title(paste(name, "log2", sep=" "))
dev.off()



# Z-NORMALIZE and plot
dz <- scale(dt)
# Make density plot of z-normalized
png(file.path(document.dir,sprintf('%s-znormlog2.png',instance.name)))
y.max <- max(apply(dz, 2, function(x) max(density(x)$y)))
x.min <- min(apply(dz, 2, function(x) min(density(x)$x)))
x.max <- max(apply(dz, 2, function(x) max(density(x)$x)))
plot(0,0, ylim=c(0,y.max), xlim=c(x.min,x.max), type='n')
for(i in 1:ncol(dz)) {
	lines(density(dz[,i]), col=line.colors[i])
}
title(paste(name, "z-normalized + log2", sep=" "))
dev.off()



# TRANSFORM data back to LINEAR scale
d2 <- 2^dt

# density plot for linear scale
png(file.path(document.dir,sprintf('%s-linear.png',instance.name)))
if(name == "miRNA"){
	y.max <- 6
} else {
	y.max <- max(apply(d2, 2, function(x) max(density(x)$y)))
}
x.min <- min(apply(d2, 2, function(x) min(density(x)$x)))
x.max <- max(apply(d2, 2, function(x) max(density(x)$x)))
plot(0,0, ylim=c(0,y.max), xlim=c(x.min,x.max), type='n')
for(i in 1:ncol(d2)) {
	lines(density(d2[,i]), col=line.colors[i])
}
title(paste(name, "raw scale", sep=" "))
dev.off()



# SCALE data to [0,1]
mins   <- apply(d2,2,min)
ranges <- apply(d2,2,function(x)diff(range(x)))
ds <- scale(d2, center=mins, scale=ranges) + SIGMA
# BACK TO LOG2 scale
ds <- log2(ds)

# density plot for [0,1]+log2 scale
png(file.path(document.dir,sprintf('%s-scaledlog2.png',instance.name)))
if(name == "miRNA") {
	y.max <- 6
} else {
	y.max <- max(apply(ds, 2, function(x) max(density(x)$y)))
}
x.min <- min(apply(ds, 2, function(x) min(density(x)$x)))
x.max <- max(apply(ds, 2, function(x) max(density(x)$x)))
plot(0,0, ylim=c(0,y.max), xlim=c(x.min,x.max), type='n')
for(i in 1:ncol(ds)) {
	lines(density(ds[,i]), col=line.colors[i])
}
title(paste(name, "[0,1]-scaled + log2", sep=" "))
dev.off()



# OUTPUT
table.out <- ds
document.out <- ""

# Output also mins and ranges
mins.df <- data.frame(ID=names(mins), min=mins)
write.table(mins.df, file=get.output(cf,'optOut1'), sep='\t', row.names=FALSE)
ranges.df <- data.frame(ID=names(ranges), range=ranges)
write.table(ranges.df, file=get.output(cf,'optOut2'), sep='\t', row.names=FALSE)
save(mins,ranges,table.out, file=get.output(cf,'optOut3'))
rm(optOut1, optOut2, optOut3)
