# This script processes interactions from curated datasets related to Alzheimer's disease from IntAct with mi score >0.45

library(gProfileR)

# Create the folder where current results will be written
resdir<-paste("results","intact",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Read in text file of all interactions with miscore >0.45 downloaded from IntAct and save them as RData
# Please indicate the path to the downloaded data
print("Loading the data")
alz_intact <- read.delim("../../data/intact/alzheimers_intact_v_4_2_6.txt", stringsAsFactors=F, sep = "\t", header=T)

print("Performing analysis")
# Select only interactions related to human
alz_intact <- alz_intact[alz_intact$Taxid.interactor.A%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]
alz_intact <- alz_intact[alz_intact$Taxid.interactor.B%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]

# Filter out interactions with the scores < 0.45
alz_intact$Confidence.value.s. <- as.numeric(sapply(strsplit(alz_intact$Confidence.value.s., split = "miscore:"),'[',2))

alz_intact <- alz_intact[alz_intact$Confidence.value.s. >= 0.45,]

# Convert interactors' names to ensg
# Get interactors' names
A <- unique(as.character(alz_intact[,1]))
B <- unique(as.character(alz_intact[,2]))

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

# Part of Alzheimer's disease IntAct dataset with uncharacterised proteins

# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from alz_intact
alz_intact_ABVAL <- alz_intact[,c(1,2,15)]
alz_intact_ABVAL[,1] <- as.character(alz_intact_ABVAL[,1])
alz_intact_ABVAL[,2] <- as.character(alz_intact_ABVAL[,2])
alz_intact_ABVAL[,3] <- as.character(alz_intact_ABVAL[,3])
colnames(alz_intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
alz_intact_ABVAL$A <- sapply(strsplit(alz_intact_ABVAL$A, split = ":"),'[',2)
alz_intact_ABVAL$B <- sapply(strsplit(alz_intact_ABVAL$B, split = ":"),'[',2)

# Select the rows where there are unconverted proteins in "A" or "B" columns
rn_A_UP <- rownames(alz_intact_ABVAL[alz_intact_ABVAL$A%in%UP,])
rn_B_UP <- rownames(alz_intact_ABVAL[alz_intact_ABVAL$B%in%UP,])
rn_AB_UP <- unique(c(rn_A_UP,rn_B_UP))

# Select the part of the data with uncharacterised proteins
alz_intact_ucppi <- alz_intact_ABVAL[rownames(alz_intact_ABVAL)%in%rn_AB_UP,]
alz_intact_ucppi$score<-as.numeric(alz_intact_ucppi$score)

# Convert protein names(where present) to ensg
# Merge for the first interactor
alz_intact_ucppi_A <- merge(alz_intact_ucppi, AB2ensg, by.x = "A", by.y = ".id", all.x = T)

# Merge for the second interactor
alz_intact_ucppi_AB <- merge(alz_intact_ucppi_A, AB2ensg, by.x = "B", by.y = ".id", all.x = T)

# Convert to character
alz_intact_ucppi_AB$Target.x[is.na(alz_intact_ucppi_AB$Target.x)] <- as.character(alz_intact_ucppi_AB$A[is.na(alz_intact_ucppi_AB$Target.x)])
alz_intact_ucppi_AB$Target.y[is.na(alz_intact_ucppi_AB$Target.y)] <- as.character(alz_intact_ucppi_AB$B[is.na(alz_intact_ucppi_AB$Target.y)])

# Rename the columns
colnames(alz_intact_ucppi_AB) <- c("B", "A", "score", "ensg1", "ensg2")

# Select columns with converted names and score
alz_intact_int_ucppi <- alz_intact_ucppi_AB[, c(4,5,3)]

# Rename the columns
colnames(alz_intact_int_ucppi) <- c("ensg1", "ensg2", "score")
alz_intact_int_ucppi <- cbind(alz_intact_int_ucppi, interaction_type = "UCPPI")
alz_intact_int_ucppi <- cbind(alz_intact_int_ucppi, data_source = "ADIA") # source name of the interacions in homo sapiens with miscore >0.45 IntAct

### Protein-protein interaction(PPI) part of Alzheimer's disease IntAct dataset
# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from alz_intact
alz_intact_ABVAL <- alz_intact[, c(1,2,15)]
alz_intact_ABVAL[,1] <- as.character(alz_intact_ABVAL[,1])
alz_intact_ABVAL[,2] <- as.character(alz_intact_ABVAL[,2])
alz_intact_ABVAL[,3] <- as.character(alz_intact_ABVAL[,3])
colnames(alz_intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
alz_intact_ABVAL$A <- sapply(strsplit(alz_intact_ABVAL$A, split = ":"),'[',2)
alz_intact_ABVAL$B <- sapply(strsplit(alz_intact_ABVAL$B, split = ":"),'[',2)

# Merge for the first interactor
alz_intact_ppi_A <- merge(alz_intact_ABVAL, AB2ensg, by.x = "A", by.y = ".id", all = F)

# Merge for the second interactor
alz_intact_ppi_AB <- merge(alz_intact_ppi_A, AB2ensg, by.x = "B", by.y = ".id", all = F)

# Rename the columns
colnames(alz_intact_ppi_AB) <- c("B", "A", "score", "ensg1", "ensg2")

# Save converted data with original names and corresponding ENSG IDs
#save(alz_intact_ppi_AB, file = "alz_intact_ppi_AB.RData") #file that has fields "B","A","score","ensg1","ensg2"

# Select columns 4,5,3
alz_intact_int_ppi <- alz_intact_ppi_AB[, c(4,5,3)]

# Rename the columns
colnames(alz_intact_int_ppi) <- c("ensg1", "ensg2", "score")
alz_intact_int_ppi <- cbind(alz_intact_int_ppi, interaction_type = "PPI")
alz_intact_int_ppi <- cbind(alz_intact_int_ppi, data_source = "ADIA") # source name of the interacions in homo sapiens with miscore >0.45 IntAct

# Combined UCPPI and PPI
# Merge converted PPI and PCI part of IntAct
alz_intact_int <- rbind(alz_intact_int_ppi, alz_intact_int_ucppi)

# Convert factors to characters
df2string<-function(df){
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
df[,3]<-as.numeric(df[,3])
return (df)}

# Remove the duplicated undirrescted edges with the same score.
# For example ENSG1-ENSG2 0.5 and ENSG2-ENSG 0.5
# PPIs from IntAct
alz_intact_int<-df2string(alz_intact_int)
alz_intact_int <- alz_intact_int[!duplicated(alz_intact_int), ]

# Exclude interaction_type = "UCPPI" from the integrated dataset
alz_intact_int <- alz_intact_int[!alz_intact_int$interaction_type%in%"UCPPI",]

alz_intact_int <- data.frame(t(apply(alz_intact_int[,1:2], 1, sort)), alz_intact_int[3:5])
alz_intact_int<-df2string(alz_intact_int)
alz_intact_int <- alz_intact_int[!duplicated(alz_intact_int),]
colnames(alz_intact_int)[1:2]<-c("ensg1", "ensg2")

#Save the part of the integrated dataset from Alzheimer's disease IntAct PPI and UCPPI for human
print("Writing the results to the file")
save(alz_intact_int, file = "alz_intact_int.RData")
write.table(alz_intact_int, file = "alz_intact_int.txt", sep = "\t", quote = F, row.names = F)

setwd("../../")


