# This script processes manually curated synapse realted  dataset from IntAct

# Create the folder where current results will be written
resdir<-paste("results","intact",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Read in text file of all interactions with miscore >0.45 downloaded from IntAct and save them as RData
# Please indicate the path to the downloaded data
pring("Loading the data")
syn_intact <- read.delim("../../data/intact/synapse_intact_v_4_2_6.txt", stringsAsFactors=F, sep = "\t", header=T)

print("Performing analysis")
# Select only interactions related to human
syn_intact <- syn_intact[syn_intact$Taxid.interactor.A%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]
syn_intact <- syn_intact[syn_intact$Taxid.interactor.B%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]

# Filter out interactions with the scores <0.45
syn_intact$Confidence.value.s. <- as.numeric(sapply(strsplit(syn_intact$Confidence.value.s., split = "miscore:"),'[',2))
syn_intact <- syn_intact[syn_intact$Confidence.value.s. >= 0.45,]

# Convert interactors' names to ensg
# Get interactors' names
A <- unique(as.character(syn_intact[,1]))
B <- unique(as.character(syn_intact[,2]))

# Show the structure
A_ID <- sapply(strsplit(A,split = ":"),'[',2)
B_ID <- sapply(strsplit(B,split = ":"),'[',2)
AB <- unique(c(A_ID,B_ID))

## Convert the names to ENSG
AB2ensg<-gconvert(AB)
AB2ensg<-AB2ensg[,c(2,4)]

# Remove duplicates
AB2ensg<-AB2ensg[!duplicated(AB2ensg), ]
colnames(AB2ensg) <- c(".id", "Target")


# Non-converted IDs
AB_nc<-AB[!AB%in%AB2ensg$.id]

# Convert unconverted isoforms that have "-" in the name
AB_nc <- sapply(strsplit(AB_nc,split = "-"),'[',1)
AB_nc_original <- data.frame(cbind(AB[!AB%in%AB2ensg$.id],AB_nc))

# Comvert non-converted IDs
AB_nc2ensg<-gconvert(AB_nc)
AB_nc2ensg<-AB_nc2ensg[,c(2,4)]

# Remove duplicates
AB_nc2ensg<-AB_nc2ensg[!duplicated(AB_nc2ensg), ]
colnames(AB_nc2ensg) <- c(".id", "Target")

AB_nc2ensg_orig <- merge(AB_nc2ensg, AB_nc_original, by.x=".id", by.y="AB_nc", all=F)
AB_nc2ensg_orig <-AB_nc2ensg_orig[, 2:3]
colnames(AB_nc2ensg_orig)<-c("Target", ".id")
AB_nc2ensg_orig<-AB_nc2ensg_orig[,c(2,1)]

# Combine AB2ensg and AB_nc
AB2ensg<- rbind(AB2ensg, AB_nc2ensg_orig) 

# Interactions with uncharacterised proteins (UCPPI) part of IntAct dataset
# Place interactions with uncharacterised proteins into data structure
# Select all unconverted proteins
UP <-AB[!AB%in%AB2ensg$.id]

# Part of Synapse IntAct dataset describing PPIs with uncharacterised proteins

# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from syn_intact
syn_intact_ABVAL <- syn_intact[,c(1,2,15)]
syn_intact_ABVAL[,1] <- as.character(syn_intact_ABVAL[,1])
syn_intact_ABVAL[,2] <- as.character(syn_intact_ABVAL[,2])
syn_intact_ABVAL[,3] <- as.character(syn_intact_ABVAL[,3])
colnames(syn_intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
syn_intact_ABVAL$A <- sapply(strsplit(syn_intact_ABVAL$A, split = ":"),'[',2)
syn_intact_ABVAL$B <- sapply(strsplit(syn_intact_ABVAL$B, split = ":"),'[',2)

# Select the rows where there are unconverted proteins in "A" or "B" columns
rn_A_UP <- rownames(syn_intact_ABVAL[syn_intact_ABVAL$A%in%UP,])
rn_B_UP <- rownames(syn_intact_ABVAL[syn_intact_ABVAL$B%in%UP,])
rn_AB_UP <- unique(c(rn_A_UP,rn_B_UP))

# Select the part of the data with uncharacterised proteins
syn_intact_ucppi <- syn_intact_ABVAL[rownames(syn_intact_ABVAL)%in%rn_AB_UP,]
syn_intact_ucppi$score<-as.numeric(syn_intact_ucppi$score)

# Convert protein names(where present) to ensg
# Merge for the first interactor
syn_intact_ucppi_A <- merge(syn_intact_ucppi, AB2ensg, by.x = "A", by.y = ".id", all.x = T)

# Merge for the second interactor
syn_intact_ucppi_AB <- merge(syn_intact_ucppi_A, AB2ensg, by.x = "B", by.y = ".id", all.x = T)

# Convert to character
syn_intact_ucppi_AB$Target.x[is.na(syn_intact_ucppi_AB$Target.x)] <- as.character(syn_intact_ucppi_AB$A[is.na(syn_intact_ucppi_AB$Target.x)])
syn_intact_ucppi_AB$Target.y[is.na(syn_intact_ucppi_AB$Target.y)] <- as.character(syn_intact_ucppi_AB$B[is.na(syn_intact_ucppi_AB$Target.y)])

# Rename the columns
colnames(syn_intact_ucppi_AB) <- c("B", "A", "score", "ensg1", "ensg2")

# Select columns with converted names and score
syn_intact_int_ucppi <- syn_intact_ucppi_AB[, c(4,5,3)]

# Rename the columns
colnames(syn_intact_int_ucppi) <- c("ensg1", "ensg2", "score")
syn_intact_int_ucppi <- cbind(syn_intact_int_ucppi, interaction_type = "UCPPI")
syn_intact_int_ucppi <- cbind(syn_intact_int_ucppi, data_source = "SIA") # source name of the interacions in homo sapiens with miscore >0.45 IntAct

### Protein-protein interaction(PPI) part of IntAct dataset
# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from syn_intact

syn_intact_ABVAL <- syn_intact[, c(1,2,15)]
syn_intact_ABVAL[,1] <- as.character(syn_intact_ABVAL[,1])
syn_intact_ABVAL[,2] <- as.character(syn_intact_ABVAL[,2])
syn_intact_ABVAL[,3] <- as.character(syn_intact_ABVAL[,3])
colnames(syn_intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
syn_intact_ABVAL$A <- sapply(strsplit(syn_intact_ABVAL$A, split = ":"),'[',2)
syn_intact_ABVAL$B <- sapply(strsplit(syn_intact_ABVAL$B, split = ":"),'[',2)

# Merge for the first interactor
syn_intact_ppi_A <- merge(syn_intact_ABVAL, AB2ensg, by.x = "A", by.y = ".id", all = F)

# Merge for the second interactor
syn_intact_ppi_AB <- merge(syn_intact_ppi_A, AB2ensg, by.x = "B", by.y = ".id", all = F)

# Rename the columns
colnames(syn_intact_ppi_AB) <- c("B", "A", "score", "ensg1", "ensg2")

# Select columns 4,5,3
syn_intact_int_ppi <- syn_intact_ppi_AB[, c(4,5,3)]

# Rename the columns
colnames(syn_intact_int_ppi) <- c("ensg1", "ensg2", "score")
syn_intact_int_ppi <- cbind(syn_intact_int_ppi, interaction_type = "PPI")
syn_intact_int_ppi <- cbind(syn_intact_int_ppi, data_source = "SIA") # source name of the interacions in homo sapiens with miscore >0.45 IntAct

### Combined UCPPI and PPI
# Merge converted PPI and PCI part of IntAct
syn_intact_int <- rbind(syn_intact_int_ppi, syn_intact_int_ucppi)

# Convert factors to characters
df2string<-function(df){
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
df[,3]<-as.numeric(df[,3])
return (df)}

# Remove the duplicated undirrescted edges with the same score.
# For example ENSG1-ENSG2 0.5 and ENSG2-ENSG 0.5
# PPIs from IntAct
syn_intact_int<-df2string(syn_intact_int)
syn_intact_int <- syn_intact_int[!duplicated(syn_intact_int), ]
#Exclude interaction_type = "UCPPI" from the integrated dataset
syn_intact_int <- syn_intact_int[!syn_intact_int$interaction_type%in%"UCPPI",]

syn_intact_int <- data.frame(t(apply(syn_intact_int[,1:2], 1, sort)), syn_intact_int[3:5])
syn_intact_int<-df2string(syn_intact_int)
syn_intact_int <- syn_intact_int[!duplicated(syn_intact_int),]
colnames(syn_intact_int)[1:2]<-c("ensg1", "ensg2")

#Save the part of the integrated dataset from Synapse IntAct PPI and UCPPI for human
print("Writing the results to the file")
save(syn_intact_int, file = "syn_intact_int.RData")
write.table(syn_intact_int, file = "syn_intact_int.txt", sep = "\t", quote = F, row.names = F)

setwd("../../")
