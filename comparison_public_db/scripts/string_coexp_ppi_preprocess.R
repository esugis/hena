# This script selects co-expression and physical interactions from STRING datase in Homo sapiens.

# Load the data
print("Loading STRING-DB data")
string <- read.table(file = "data_db/string_db/9606.protein.links.full.v10.5.txt", sep = " ", stringsAsFactors = F,header=T)

print("Separating co-expression and PPIs")
# Leave int1-int2-coexpression-ppi
string_coexp_ppi <- string[, c(1, 2, 8, 10)]

# Remove rows where coexpression and ppi = 0
library(dplyr)

string_coexp_ppi_filt <- string_coexp_ppi%>%
 filter(coexpression > 0 | experiments > 0)

# Remove organism prefix 9606 from protain names in col1 and col2
string_coexp_ppi_filt$protein1 <- gsub("9606.", "", string_coexp_ppi_filt$protein1)
string_coexp_ppi_filt$protein2 <- gsub("9606.", "", string_coexp_ppi_filt$protein2)

# Create a vector of protein names to convert
string_proteins <- unique(c(string_coexp_ppi_filt$protein1, string_coexp_ppi_filt$protein2))

print("Converting ENSP IDs to ENSG IDs and removing genes with pseugene biotype")
# Convert ENSP protein ids to ENSG ids

library(gProfiler)
string2ensg<-gconvert(string_proteins)
string2ensg<-string2ensg[,c(2,4)]
colnames(string2ensg) <- c(".id", "Target")

# Merge for the first interactor
string_coexp_ppi_filt_ensg1 <- merge(string_coexp_ppi_filt, string2ensg, by.x = "protein1", by.y = ".id", all = F)

string_coexp_ppi_filt_ensg12 <- merge(string_coexp_ppi_filt_ensg1, string2ensg, by.x = "protein2", by.y = ".id", all = F)

# Filter ppi>0
string_ppi_filt <- string_coexp_ppi_filt_ensg12%>%filter(experiments>0)

# Select ensg1, ensg2     
string_ppi_filt <- string_ppi_filt[, c(5,6)]
colnames(string_ppi_filt)<-c("ensg2", "ensg1")

# Remove diplicated interactions
string_ppi_filt<- string_ppi_filt[!duplicated(string_ppi_filt),]

# Save to file
print("Writing PPIs from STRING to file")
save(string_ppi_filt, file="results/string_db/string_ppi_filt.RData")
write.table(string_ppi_filt, file="results/string_db/string_ppi_filt.txt", sep = "\t", quote = F, row.names = F)

# Filter coexpression > 0 columns
string_coexp_filt <- string_coexp_ppi_filt_ensg12%>%filter(coexpression>0)                       

# Select ensg1, ensg2	  
string_coexp_filt <- string_coexp_filt[, c(5,6)]
colnames(string_ppi_filt)<-c("ensg2", "ensg1")
# Remove diplicated interactions
string_coexp_filt<- string_coexp_filt[!duplicated(string_coexp_filt),]

# Save to file
print("Writing co-expression from STRING to file")
save(string_coexp_filt, file="results/string_db/string_coexp_filt.RData")
write.table(string_coexp_filt, file="results/string_db/string_coexp_filt.txt", sep = "\t", quote = F, row.names = F)

