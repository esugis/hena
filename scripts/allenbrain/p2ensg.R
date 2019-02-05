# This script converts probeset names from allenbrain atlas whole brain microarrays to ENSG
# Convert the probe names to ensg.
# Used liraries
library(reshape2)
library(gProfileR)

# Ontologies and probe annotations are the same for 6 microarray data sets, so read files from one as below
probes <- read.csv("data/allenbrain/178236545_ds/Probes.csv")

# Keep probe_id and entrez_id  
p <- probes[,c(1,6)]
gene <- as.numeric(p$entrez_id)
gene <- unique(gene)
print ("Number of genes in the dataset")
print(length(gene))

# Convert probes to ENSG IDs
print("Converting genes to ENSG IDs")
p2ensg1 <- gconvert(gene[1:5000], numeric_ns = "ENTREZGENE_ACC")
p2ensg2 <- gconvert(gene[5001:10000], numeric_ns = "ENTREZGENE_ACC")
p2ensg3 <- gconvert(gene[10001:15000], numeric_ns = "ENTREZGENE_ACC")
p2ensg4 <- gconvert(gene[15001:20000],numeric_ns = "ENTREZGENE_ACC")
p2ensg5 <- gconvert(gene[20001:20788], numeric_ns = "ENTREZGENE_ACC")

# Combine all converted names
p2ensg <- rbind(p2ensg1,p2ensg2,p2ensg3,p2ensg4,p2ensg5)
p2ensg <- p2ensg[, c(2,4,5)]
colnames(p2ensg) <- c("entrez_id", "Target", "name")

# Remove duplicates
p2ensg <- p2ensg[!duplicated(p2ensg),]
p2ensg$entrez_id <- gsub("ENTREZGENE_ACC:", "", p2ensg$entrez_id)
p2ensg$entrez_id <- as.numeric(p2ensg$entrez_id)
p2ensg$Target <- as.character(p2ensg$Target)
p2ensg$name <- as.character(p2ensg$name)

p2ensg_orig <- merge(p, p2ensg, by.x = "entrez_id", by.y = "entrez_id", all = F)
p2ensg_orig <- p2ensg_orig[!duplicated(p2ensg_orig),]
print("Let's have a look at the converted gene names")
head(p2ensg_orig)

# Save probes and corresponding ENSG IDs to the file
print("Writing to file")
save(p2ensg_orig, file = "results/allenbrain/p2ensg_all.RData")

# Remove tmp structures
rm(p2ensg, gene, p, probes, p2ensg1, p2ensg2, p2ensg3, p2ensg4, p2ensg5)

