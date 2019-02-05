# This script processes inteacrtion from IntAct with mi score >0.45

# Create the folder where current results will be written
resdir<-paste("results","intact",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Read in text file of all interactions with miscore >0.45 downloaded from IntAct and save them as RData
print("Loading the data")
intact <- read.delim("../../data/intact/intact_hs_v_4_2_6.txt", stringsAsFactors = F, sep = "\t", header=T)

print("Performing analyis")
# Select only interactions related to human
intact <- intact[intact$Taxid.interactor.A%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]
intact <- intact[intact$Taxid.interactor.B%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]

# Convert interactors' names to ensg
# Get interactors' names
A <- unique(as.character(intact[,1]))
B <- unique(as.character(intact[,2]))

# Show the structure
A_ID <- sapply(strsplit(A,split = ":"),'[',2)
B_ID <- sapply(strsplit(B,split = ":"),'[',2)
AB <- unique(c(A_ID,B_ID))

# Convert the names to ENSG
AB2ensg<-gconvert(AB)
AB2ensg<-AB2ensg[,c(2,4)]

# Remove duplicates
AB2ensg<-AB2ensg[!duplicated(AB2ensg), ]
colnames(AB2ensg) <- c(".id", "Target")


# Interactions with uncharacterised proteins (UCPPI) part of IntAct dataset
# Place interactions with uncharacterised proteins into data structure
# Select all unconverted proteins
UP <-AB[!AB%in%AB2ensg$.id]

# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from intact
intact_ABVAL <- intact[,c(1,2,15)]
intact_ABVAL[,1] <- as.character(intact_ABVAL[,1])
intact_ABVAL[,2] <- as.character(intact_ABVAL[,2])
intact_ABVAL[,3] <- as.character(intact_ABVAL[,3])
colnames(intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
intact_ABVAL$A <- sapply(strsplit(intact_ABVAL$A, split = ":"),'[',2)
intact_ABVAL$B <- sapply(strsplit(intact_ABVAL$B, split = ":"),'[',2)
intact_ABVAL$score <- sapply(strsplit(intact_ABVAL$score, split = "miscore:"),'[',2)

# Select the rows where there are unconverted proteins in "A" or "B" columns
rn_A_UP <- rownames(intact_ABVAL[intact_ABVAL$A%in%UP,])
rn_B_UP <- rownames(intact_ABVAL[intact_ABVAL$B%in%UP,])
rn_AB_UP <- unique(c(rn_A_UP,rn_B_UP))

# Select the part of the data with uncharacterised proteins
intact_ucppi <- intact_ABVAL[rownames(intact_ABVAL)%in%rn_AB_UP,]
intact_ucppi$score<-as.numeric(intact_ucppi$score)

# Convert protein names(where present) to ensg
# Merge for the first interactor
intact_ucppi_A <- merge(intact_ucppi, AB2ensg, by.x = "A", by.y = ".id", all.x = T)

# Merge for the second interactor
intact_ucppi_AB <- merge(intact_ucppi_A, AB2ensg, by.x = "B", by.y = ".id", all.x = T)

# Convert to character
intact_ucppi_AB$Target.x[is.na(intact_ucppi_AB$Target.x)] <- as.character(intact_ucppi_AB$A[is.na(intact_ucppi_AB$Target.x)])
intact_ucppi_AB$Target.y[is.na(intact_ucppi_AB$Target.y)] <- as.character(intact_ucppi_AB$B[is.na(intact_ucppi_AB$Target.y)])

# Rename the columns
colnames(intact_ucppi_AB) <- c("B", "A", "score", "ensg1", "ensg2")

# Select columns with converted names and score
intact_int_ucppi <- intact_ucppi_AB[, c(4,5,3)]

# Rename the columns
colnames(intact_int_ucppi) <- c("ensg1", "ensg2", "score")
intact_int_ucppi <- cbind(intact_int_ucppi, interaction_type = "UCPPI")
intact_int_ucppi <- cbind(intact_int_ucppi, data_source = "IAH") # source name of the interacions in homo sapiens with miscore >0.45 IntAct

### Protein-protein interaction(PPI) part of IntAct dataset
# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from intact
intact_ABVAL <- intact[, c(1,2,15)]
intact_ABVAL[,1] <- as.character(intact_ABVAL[,1])
intact_ABVAL[,2] <- as.character(intact_ABVAL[,2])
intact_ABVAL[,3] <- as.character(intact_ABVAL[,3])
colnames(intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
intact_ABVAL$A <- sapply(strsplit(intact_ABVAL$A, split = ":"),'[',2)
intact_ABVAL$B <- sapply(strsplit(intact_ABVAL$B, split = ":"),'[',2)
intact_ABVAL$score <- sapply(strsplit(intact_ABVAL$score, split = "miscore:"),'[',2)

# Merge for the first interactor
intact_ppi_A <- merge(intact_ABVAL, AB2ensg, by.x = "A", by.y = ".id", all = F)

# Merge for the second interactor
intact_ppi_AB <- merge(intact_ppi_A, AB2ensg, by.x = "B", by.y = ".id", all = F)

# Rename the columns
colnames(intact_ppi_AB) <- c("B", "A", "score", "ensg1", "ensg2")

# Select columns 4,5,3
intact_int_ppi <- intact_ppi_AB[, c(4,5,3)]

# Rename the columns
colnames(intact_int_ppi) <- c("ensg1", "ensg2", "score")
intact_int_ppi <- cbind(intact_int_ppi, interaction_type = "PPI")
intact_int_ppi <- cbind(intact_int_ppi, data_source = "IAH") # source name of the interacions in homo sapiens with miscore >0.45 IntAct

### Combined UCPPI and PPI
# Merge converted PPI and PCI part of IntAct
intact_int <- rbind(intact_int_ppi, intact_int_ucppi)

# Convert factors to characters
df2string<-function(df){
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
df[,3]<-as.numeric(df[,3])
return (df)}

# Remove the duplicated undirrescted edges with the same score.
# For example ENSG1-ENSG2 0.5 and ENSG2-ENSG 0.5
# PPIs from IntAct
intact_int<-df2string(intact_int)
intact_int <- intact_int[!duplicated(intact_int), ]

# Exclude interaction_type = "UCPPI" from the integrated dataset
intact_int <- intact_int[!intact_int$interaction_type%in%"UCPPI",]

intact_int <- data.frame(t(apply(intact_int[,1:2], 1, sort)), intact_int[3:5])
intact_int<-df2string(intact_int)
intact_int <- intact_int[!duplicated(intact_int),]
colnames(intact_int)[1:2]<-c("ensg1", "ensg2")

#Save the part of the integrated dataset from IntAct PPI for human
print("Writing the results to the file")
save(intact_int, file = "intact_int.RData")
write.table(intact_int, file = "intact_int.txt", sep = "\t", quote = F, row.names = F)

setwd("../../")


