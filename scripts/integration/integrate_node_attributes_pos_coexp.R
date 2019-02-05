# This script assembles the nodes attribute file from GWAS, Allen Brain atlas expression, 
# and  adds biotype of all the genes(gwas, allen brain atlas, IntAct, pba, epistasis, co-expression)
# This script create node attributes for GraphSage analysis. Only positive coexpression in brain regions is considered.

# Create the folder where current results will be written
resdir<-paste("~/absb/results","integration",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Load individual  datasets of gene attributes

load(file="~/absb/results/gwas/gwas_ensg.RData")

# Filter based on p-value < 0.05
gwas_ensg <- (gwas_ensg[gwas_ensg$GWAS_pvalue<=0.05,])

# Extract node IDs
gwas_nodes <- unique(gwas_ensg$ensg)


# Create nodes with labels alzheimer from GWAS genes and nodes from Intact Alzheimer's
## IntAct Alzheimer's related manually curated PPI dataset with all MIscores
load(file="~/absb/results/intact/alz_intact_int_no_score_filtering.RData")
alz_int_nodes <- unique(c(alz_intact_int$ensg1, alz_intact_int$ensg2))

# Combine IDs from GWAS and Intact Alzheimer's
alz_nodes <- unique(c(alz_int_nodes, gwas_nodes)) #2116
write.table(alz_nodes , file="alz_nodes.txt", sep="\t", quote=F, row.names=F)



load(file="~/absb/results/allenbrain/mtx_zscores_max_sorted.RData")
load(file="~/absb/results/ps/ps.RData")

# Load integrated interactions
load(file="~/absb/results/integration/integrated_int_pos_coexp.RData")

# Add gene names and biotype to the interactors
rownames(integrated_int) <- 1:nrow(integrated_int)

library(biomaRt)
library(tidyr)
library(dplyr)

ensg_integrated_int <- unique(c(integrated_int$ensg.A,integrated_int$ensg.B))


# Extract ensg ids

# GWAS
ensg_gwas <- gwas_nodes
length(ensg_gwas)

# Positve selection
ensg_ps <- unique(as.character(as.vector(ps$ensg)))
length(ensg_ps)

# Integrated interactions
ensg_integrated_int <- unique(c(integrated_int$ensg.A,integrated_int$ensg.B))

# Combine all ensg ids
ensg <- unique(c(ensg_integrated_int, ensg_ps, ensg_allen, ensg_gwas))
length(ensg) 

# Get gene names and biotype
# Biomart
library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="ensembl.org")
ensg2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "gene_biotype"), filters=c("ensembl_gene_id"),values=ensg, mart=mart)
colnames(ensg2name) <- c("ensg", "gene_name", "biotype")

# Save converted ensg ids and corresponding name
write.table(ensg2name, file = "ensg2name.txt", sep="\t", quote=F, row.names=F)
save(ensg2name,file = "ensg2name.RData")

# Merge individual attribute data together based on ensg
# Combine gwas and positive selection data 
node_attr <- merge(gwas_ensg, ps, by.x="ensg",by.y="ensg", all=T)

# Combine with allen brain expression
node_attr <- merge(node_attr, allen, by.x="ensg",by.y="ensg", all=T)

# Combine with gene names and biotype
node_attr <- merge(ensg2name, node_attr, by.x="ensg",by.y="ensg", all=T)

# Remove duplicates
node_attr<-node_attr[!duplicated(node_attr), ]
dim(node_attr)#  12551 4


# Merge with the rest node attributes
ensg_int <- cbind(ensg_integrated_int,ensg_integrated_int)
colnames(ensg_int) <- c("ensg_rm", "ensg")

node_attr <- merge(ensg_int, node_attr, by.x="ensg",by.y="ensg", all=T)
node_attr <- node_attr[, -2]
node_attr <- node_attr[!duplicated(node_attr), ]
dim(node_attr) #45248   420

# Add feature column alz/nonalz

nodes <- node_attr%>%
  mutate(node_type = ifelse(ensg%in%alz_nodes)==T , "nonalz", "alz"))%>%
  select(c(1,421,2:420))




# Save to file
write.table(nodes, file="node_attributes.txt", sep="\t", quote=F, row.names=F )
save(nodes, file="node_attributes.RData")








