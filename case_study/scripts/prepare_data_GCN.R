# This script combines all indiviadual interaction datasets into one and removes duplicated interaction of one type
# Only positive co-expression in brain regions is used for the analysis GraphSage
# Create the folder where current results will be written
resdir<-paste("case_study/datasets/genes_data","integration",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Co-expression from Allen Brain Atlas coexpression values are > 0.5
# IGRI are excluded
# ENSGs are kept

############### Prepare interactions ################
# Load integrated dataset
print("Loading interactions")
load(file="case_study/data/integrated_int.RData")
print("Preparing intercations")

# Exclude IGRI
int<- integrated_int[!integrated_int$interaction_type%in%"IGRI",]

# Set coexpression cutoff in ABA datasets to >0.5
int_aba_coexp<-int[int$data_source%in%c("ABA_CA1",
"ABA_CA2", "ABA_CA3", "ABA_CA4", "ABA_SptN", "ABA_subiculum","ABA_DG"),]
int_aba_coexp<-int_aba_coexp[int_aba_coexp$score>=0.5,]

# Combine filetered dataset
int<- int[!int$data_source%in%c("ABA_CA1",
"ABA_CA2", "ABA_CA3", "ABA_CA4", "ABA_SptN", "ABA_subiculum","ABA_DG"),]
int_filt <- rbind(int,int_aba_coexp)
rm(int,int_aba_coexp)

############### Prepare attributes ###############################################################
print("Loading node attributes")
load(file="case_study/data/node_attributes.RData")
ds_ensg<-unique(c(int_filt$ensg.A, int_filt$ensg.B))
node_attributes<-node_attributes[node_attributes$ensg%in%ds_ensg,]

# Remove rows with IGR
print("Removing nodes representing inter genic regions")
biotype_rm <-unique(node_attributes$biotype[grep("-",node_attributes$biotype)])
node_attributes<-node_attributes[!node_attributes$biotype%in%biotype_rm,]
ensg<-unique(node_attributes$ensg[grep("ENSG",node_attributes$ensg)])
node_attributes<-node_attributes[node_attributes$ensg%in%ensg,]
node_attributes<-node_attributes[,-c(2,3,4,6)]

# Create labels: positive class is based on GWAS and Alzheimer's ppi from IntAct
# Disease label ==1 if node has GWAS p-value
print("Adding node labels: disease, non-disease and unknown")
node_attributes$GWAS_pvalue<-as.numeric(node_attributes$GWAS_pvalue)
node_attributes[!is.na(node_attributes$GWAS_pvalue),c(2)]<-"disease"
node_attributes[is.na(node_attributes$GWAS_pvalue),c(2)]<-"unknown"
node_attributes<-node_attributes[!duplicated(node_attributes),]
colnames(node_attributes)[2]<-"label"

# Disease label if node is in the list of nodes from IntAct Alzheimer's dataset
int_alz <-int[int$data_source%in%"ADIA",]
nodes_alz<-unique(c(int_alz$ensg.A, int_alz$ensg.B))
node_attributes[node_attributes$ensg%in%nodes_alz,2]<-"disease"

# Crerate negative class based on the evolutionary studies.
gene_status<-read.csv(file="case_study/data/negative_class/gene_status.csv",sep=";",stringsAsFactors=F)
gene_nondis<- unique(gene_status[gene_status$Group%in%"END",2])

# Add non-disease class
node_attributes[node_attributes$ensg%in%gene_nondis&node_attributes$label%in%"unknown",2]<-"non-disease"

# Find genes that are more expressed in the gene regions of interest
# Import normalised expression from Allen brain datasets and apply Wilcoxon test on samples in selected regions vs all the rest.
load(file="case_study/datasets/genes_data/wlx_genes.RData")

wlx_genes<-wlx_genes[,c(2,4)]
node_attributes<-merge(node_attributes,wlx_genes, by.x="ensg", by.y="Target", all.x=T )
node_attributes<-node_attributes[, c(1,2,234,3:233)]
colnames(node_attributes)[3]<-"dis_brain_reg"

node_attributes[!is.na(node_attributes$dis_brain_reg),3]<-1
node_attributes[is.na(node_attributes$dis_brain_reg),3]<-0

################# Dataset check ##################
# Keep interaction related to the nodes from node_attributes
node_ensg<-unique(node_attributes$ensg)
int_filt_1<-int_filt[int_filt$ensg.A%in%node_ensg,]
int_filt_12<-int_filt_1[int_filt_1$ensg.B%in%node_ensg,]
ds_ensg_filt<-unique(c(int_filt_12$ensg.A, int_filt_12$ensg.B))
print("Number of nodes in the data set")
length(ds_ensg_filt)

interactions <-int_filt_12
interactions<-interactions[, colnames(interactions)%in%c("ensg.A", "ensg.B", "interaction_type")]
colnames(interactions)<-c("ensg1","ensg2","int_type")
interactions<-interactions[!duplicated(interactions),]
print("Number of interactions in the data set")
dim(interactions)[1]

################# Write to file ##################
print("Writing prepared data sets as RData and csv. Use csv files in case study.")
save(interactions, file="case_study/datasets/genes_data/interactions.RData")
write.table(interactions, file="case_study/datasets/genes_data/interactions.csv", sep=" ", quote=T, row.names = F)

#Save node attributes as csv for the case study
save(node_attributes, file="case_study/datasets/genes_data/node_attributes.RData")

# Save node attributes as csv for the case study
node_attributes$ensg<-paste('"',node_attributes$ensg,'"', sep="")
node_attributes$label<-paste('"',node_attributes$label,'"', sep="")
colnames(node_attributes)<-paste('"',colnames(node_attributes),'"',sep="")

write.table(node_attributes, file="case_study/datasets/genes_data/node_attributes.csv", sep=" ", quote=F, row.names = F)
