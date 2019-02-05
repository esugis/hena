# This script assembles the nodes attribute file from GWAS, Allen Brain atlas expression, 
# and  adds biotype of all the genes(gwas, allen brain atlas, IntAct, pba, epistasis, co-expression)

# Attach libraries
library(biomaRt)
library(tidyr)
library(dplyr)

# Create the folder where current results will be written
resdir<-paste("results","integration",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Load individual  datasets of gene attributes
print("Loading individual daasets containing node attributes")
load(file="../gwas/gwas_ensg.RData")
load(file="../ps/ps.RData")
load(file= "../allenbrain/brain_tissue_zscores_aggreg.RData")

# Make aggregated z-scored in 231 brain regions wide format
brain_z <- spread(brain_tissue_zscores_aggreg, structure_acronym,zscore)
brain_z <- brain_z[, -2]

# Load integrated interactions
print("Loading integrated dataset (list of interactions)")
load(file="integrated_int.RData")

# Add gene names and biotype to the interactors
rownames(integrated_int) <- 1:nrow(integrated_int)

# Extract all nodes
print("Extracting the list of unique nodes from the list of interactions.")
nodes<-data.frame(unique(c(integrated_int$ensg.A, integrated_int$ensg.B)))
colnames(nodes)<-"node.ensg"

############## Attributes of the nodes related to intergenic regios ###################
print("Integrating node attributes.")
# Extract nodes that correspond to Intergenic regios
nodes_igri_vect<-as.character(nodes$node.ensg[grep("\\ENSG.*-",nodes$node.ensg)] )
nodes_igri<-data.frame(nodes_igri_vect)
colnames(nodes_igri)<-"node.ensg"

# Split ensg for further work with IGRI nodes as well
nodes_igri_split <-  nodes_igri %>% separate(node.ensg, c("ENSG.1", "ENSG.2"), "-")

# Select all individual ENSG
ensg_nodes_igri <- unique(c(nodes_igri_split$ENSG.1,nodes_igri_split$ENSG.2))

# Convert to gene names
print("Fetching gene name and gene biotype for the list of nodes.")
mart <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "ensembl.org")
ensg2name_igri <- getBM(attributes=c("ensembl_gene_id","external_gene_name","gene_biotype"), filters=c("ensembl_gene_id"),values=ensg_nodes_igri, mart=mart)
colnames(ensg2name_igri) <- c("ensg", "gene_name", "biotype")

# Combined Gene Names with the original data
nodes_igri_split<- merge(nodes_igri_split, ensg2name_igri, by.x = "ENSG.1", by.y = "ensg", all.x=T)
colnames(nodes_igri_split) <- gsub("gene_name","gene_name.1", colnames(nodes_igri_split))
colnames(nodes_igri_split) <- gsub("biotype","biotype.1", colnames(nodes_igri_split))
nodes_igri_split <- merge(nodes_igri_split, ensg2name_igri, by.x = "ENSG.2", by.y = "ensg", all.x=T)
colnames(nodes_igri_split) <- gsub("\\gene_name$","gene_name.2", colnames(nodes_igri_split))
colnames(nodes_igri_split) <- gsub("biotype$","biotype.2", colnames(nodes_igri_split))

# Combine Gene names from IGRI
nodes_igri_split <- cbind(nodes_igri_split, gene_name=paste(nodes_igri_split$gene_name.1, nodes_igri_split$gene_name.2, sep="-"))
nodes_igri_split <- cbind(nodes_igri_split, biotype=paste(nodes_igri_split$biotype.1, nodes_igri_split$biotype.2, sep="-"))

# Add GWAS
nodes_igri_split <- merge(nodes_igri_split, gwas_ensg, by.x = "ENSG.1", by.y = "ensg", all.x=T)
colnames(nodes_igri_split) <- gsub("GWAS_pvalue","GWAS_pvalue.1", colnames(nodes_igri_split))
colnames(nodes_igri_split) <- gsub("SNP_id","SNP_id.1", colnames(nodes_igri_split))

nodes_igri_split <- merge(nodes_igri_split, gwas_ensg, by.x = "ENSG.2", by.y = "ensg", all.x=T)
colnames(nodes_igri_split) <- gsub("\\GWAS_pvalue$","GWAS_pvalue.2", colnames(nodes_igri_split))
colnames(nodes_igri_split) <- gsub("\\SNP_id$","SNP_id.2", colnames(nodes_igri_split))

# Create colunms corresponding to gene names
nodes_igri_split <- cbind(nodes_igri_split, GWAS_pvalue=paste(nodes_igri_split$GWAS_pvalue.1, nodes_igri_split$GWAS_pvalue.2, sep="-"))
nodes_igri_split <- cbind(nodes_igri_split, SNP_id=paste(nodes_igri_split$SNP_id.1, nodes_igri_split$SNP_id.2, sep="-"))

# Add positive selection
nodes_igri_split <- merge(nodes_igri_split, ps, by.x = "ENSG.1", by.y = "ensg", all.x=T)
colnames(nodes_igri_split) <- gsub("ps_pval","ps_pval.1", colnames(nodes_igri_split))

nodes_igri_split <- merge(nodes_igri_split, ps, by.x = "ENSG.2", by.y = "ensg", all.x=T)
colnames(nodes_igri_split) <- gsub("\\ps_pval$","ps_pval.2", colnames(nodes_igri_split))

# Create colunms corresponding to gene names
nodes_igri_split <- cbind(nodes_igri_split, ps_pvalue=paste(nodes_igri_split$ps_pval.1, nodes_igri_split$ps_pval.2, sep="-"))

# Create colunms corresponding to ensg
nodes_igri_split <- cbind(nodes_igri_split, ensg=paste(nodes_igri_split$ENSG.1, nodes_igri_split$ENSG.2, sep="-"))

# Select ensg, gene_name, biotype,GWAS_pvalue, SNP_id
nodes_igri<- nodes_igri_split[,colnames(nodes_igri_split)%in%c("ensg","gene_name","biotype", "GWAS_pvalue","SNP_id", "ps_pvalue")]

# Reorder
nodes_igri <- nodes_igri[ ,c(6,1,2,4,3,5)]
nodes_igri$biotype <- as.character(nodes_igri$biotype)

# Select genes related to "pseudo_gene" biotype
nodes_igri_pseudo <-  nodes_igri[nodes_igri$biotype%in%nodes_igri$biotype[grep("pseudo.*", nodes_igri$biotype)],]
nodes_igri_rm <- as.character(nodes_igri_pseudo$ensg)

### Remove these nodes and theirintercations from the list of nodes and the final integrated dataset
nodes_igri<- nodes_igri[!nodes_igri$ensg%in%nodes_igri_rm,]

# Bind zscores based on the expression in the brain regions
nodes_igri_z <- merge(nodes_igri,brain_z, by.x="ensg", by.y="ensg", all.x=T)
rm(nodes_igri)

############## Attributes of the normal nodes ###################

# Create data frame for all the rest nodes
nodes_norm <- data.frame(nodes[!nodes$node.ensg%in%nodes_igri_vect,])
colnames(nodes_norm)<-"node.ensg"
nodes_norm$node.ensg <- as.character(nodes_norm$node.ensg)
# Select all ensg
ensg_nodes<-as.character(nodes_norm$node.ensg)

# Convert to gene names
mart <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "ensembl.org")
ensg2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name","gene_biotype"), filters=c("ensembl_gene_id"),values=ensg_nodes, mart=mart)
colnames(ensg2name) <- c("ensg", "gene_name", "biotype")

# Combined Gene Names with the original data
nodes_norm<- merge(nodes_norm, ensg2name, by.x = "node.ensg", by.y = "ensg", all.x=T)

# Add GWAS
nodes_norm <- merge(nodes_norm, gwas_ensg, by.x = "node.ensg", by.y = "ensg", all.x=T)

# Add positive selection
nodes_norm <- merge(nodes_norm, ps, by.x = "node.ensg", by.y = "ensg", all.x=T)
colnames(nodes_norm)[6]<-c("ps_pvalue")

# Select genes related to "pseudo_gene" biotype
print("Removing nodes with pseudogene biotype.")
 nodes_norm$biotype <- as.character(nodes_norm$biotype)
nodes_norm_pseudo <-  nodes_norm[nodes_norm$biotype%in%nodes_norm$biotype[grep("pseudo.*", nodes_norm$biotype)],]
nodes_norm_rm <- as.character(nodes_norm_pseudo$node.ensg)

### Remove these nodes and theirintercations from the list of nodes and the final integrated dataset
nodes_norm <- nodes_norm[!nodes_norm$node.ensg%in%nodes_norm_rm,]
colnames(nodes_norm)[1]<-"ensg"

# Bind zscores based on the expression in the brain regions
nodes_norm_z <- merge(nodes_norm,brain_z, by.x="ensg", by.y="ensg", all.x=T)
rm(nodes_norm,brain_z)

############## Combine and save node attributes ##############
node_attributes<-rbind(nodes_norm_z, nodes_igri_z)

# Remove duplicates
node_attributes<-node_attributes[!duplicated(node_attributes), ]

# Save to file
print("Writing node attributes to file.")
write.table(node_attributes, file="node_attributes.txt", sep="\t", quote=F, row.names=F )
save(node_attributes, file="node_attributes.RData")
rm(node_attributes)

############## Remove interactions that contain genes with "pseudogene" biotype #################
# Combine nodes to be removed
print("Updating the list of interactions.")
print("Removing genes with pseudogene biotype from the integrated dataset")
nodes_rm<-unique(c(nodes_norm_rm,nodes_igri_rm))

# Remove the interactions containing the nodes to be removed
integrated_int<-integrated_int[!integrated_int$ensg.A%in%nodes_rm,]
integrated_int<-integrated_int[!integrated_int$ensg.B%in%nodes_rm,]

# Save to file
print("Writing updated integrated dataset to file.")
write.table(integrated_int, file="integrated_int.txt", sep="\t", quote=F, row.names=F )
save(integrated_int, file="integrated_int.RData") 

setwd("../../")
