### Just testing the data, comparing published and old

DATADIR = "~/wrk/dippa-data/"
OLDDATADIR = "~/wrk/vessela_brca/data/"

# read hospital data, replace minus signs in patient IDs
hosp <- read.delim(file.path(OLDDATADIR,"OSLO2","key_array_patient.csv"))
hosp[,1] <- sub("-", ".", hosp[,1])



### PROTEIN DATA ###################

# load protein data
prot <- read.delim(file.path(DATADIR, "Oslo2-RPPA_data.csv"))
prot[1:2,1:3]
dim(prot)
# load old protein data
protold <- read.delim(file.path(OLDDATADIR, "OSLO2", "protein.txt"))
protold[1:2,1:3]
dim(protold)
# old has more proteins and samples, see which are not included in new
missingprotsamp <- setdiff(names(protold),names(prot))
missingprotsamp
missingprot.new <- setdiff(protold[,1], prot[,1])
missingprot.old <- setdiff(prot[,1], protold[,1])
missingprot.new
missingprot.old
# check proteins missing from old data
# apparently with replicated antibodies, one was chosen randomly
for(p in missingprot.old) {
    s <- substr(p, 1, 4)
    show(prot[grep(s,prot[,1]),1:5])
    show(protold[grep(s,protold[,1]),1:4])
}
# check if proteins missing from new are only phosphorylated ones
nophosph <- grep("_p", missingprot.new, invert=T, value=T)
# further remove all replicated antibodies (ending _1, _2 or _3)
grep("_[1-3]$", nophosph, invert=T, value=T)
# Bcl_xL, Caspase and PARP seem to be only ones not phosphorylated or replicated
show(prot[grep("Caspa",prot[,1]),1:4])
show(protold[grep("Caspa",protold[,1]),1:4])
show(prot[grep("Bcl",prot[,1]),1:4])
show(protold[grep("Bcl",protold[,1]),1:4])
show(prot[grep("PARP",prot[,1]),1:4])
show(protold[grep("PARP",protold[,1]),1:4])
### FOR SOME REASON THESE PROTEINS HAVE JUST BEEN DROPPED

# See which hospital missing samples were from

hosp[match(missingprotsamp, hosp[,1]),]
## THESE ARE NOT FROM AHUS, WHY DROPED???
ahus <- hosp[which(hosp$Hospital == "Ahus"),1]
any(ahus %in% names(prot))
any(ahus %in% names(protold))
### NO AHUS SAMPLES IN EITHER OLD OR NEW DATA
