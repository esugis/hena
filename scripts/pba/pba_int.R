# This script processes dataset of protein-protein inteacrtions related to brain ageing (PBA)

library(gProfileR)
## Create the folder where current results will be written
resdir <- paste("results", "pba", sep = "/")
dir.create(file.path(resdir), showWarnings  =  FALSE,  recursive  =  TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Read in the data
print("Loading the data")
pba_ppi.hs <- read.table(file = "../../data/pba/PBA_PPI_HS.txt", sep = "\t", header = T, stringsAsFactors = F)
pba_ppi.hs<-pba_ppi.hs[, c(2,3,1,4)]
print("Performing analysis")
# Convert pba_ppi.hs protein names to ENSG and bing them to the dataframe.
pba_pr <- unique(c(pba_ppi.hs[,1], pba_ppi.hs[,2]))
pba2ensg<-gconvert(pba_pr)
pba2ensg<- pba2ensg[,c(2,4)]

# Convert entrezgene IDs to ENSG
colnames(pba2ensg)<-c(".id", "Target")

# Merge for the first interactor
pba_ppi.hs.p1 = merge(pba_ppi.hs, pba2ensg, by.x = "p1", by.y = ".id", all = F)
pba_ppi.hs.p1p2 <- merge(pba_ppi.hs.p1, pba2ensg, by.x = "p2", by.y = ".id", all = F)
pba_ppi.hs.ensg <- pba_ppi.hs.p1p2[, c(5,6,3)]

# Bind additional columns
pba_ppi.hs_int <- cbind(pba_ppi.hs.ensg, interaction_type = "PPI")
pba_ppi.hs_int <- cbind(pba_ppi.hs_int, data_source = "PBA")#  evidence code for Hybrigenics experimental interactions
colnames(pba_ppi.hs_int)[c(1,2,3)] <- c("ensg1","ensg2","score")
pba_int<- pba_ppi.hs_int

# Remove duplicates
pba_int <- pba_int[!duplicated(pba_int),]

df2string<-function(df){
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
df[,3]<-as.numeric(df[,3])
return (df)}

# PBA
pba_int <- df2string(pba_int)

# Remove the duplicated undirrescted edges with the same score.
# For example ENSG1-ENSG2 0.5 and ENSG2-ENSG1 0.5
pba_int <- pba_int[!duplicated(data.frame(t(apply(pba_int[1:2], 1, sort)), pba_int$score)),]

# Save the part of the integrated dataset related to interactions in HS.
print("Writing to file")
save(pba_int, file = "pba_int.RData")
write.table(pba_int, file = "pba_int.txt", sep = "\t", quote = F, row.names = F)

setwd("../../")

