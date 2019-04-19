# This script creates a subset of HENA for the manipulation in Cytoscape.
# It selects positive co-expression in brain regions with r>0.7

# Create the folder where current results will be written
resdir<-paste("case_study/datasets/genes_data","integration",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Co-expression from Allen Brain Atlas coexpression values are > 0.7
# IGRI are excluded
# ENSGs are kept

############### Prepare interactions ################
# Load integrated dataset
print("Loading interactions")
#load(file="case_study/data/integrated_int.RData")
load(file="case_study/data/integrated_int.RData")
print("Preparing intercations")

# Exclude IGRI
int<- integrated_int[!integrated_int$interaction_type%in%"IGRI",]

# Set coexpression cutoff in ABA datasets to >=0.7 and <=-0.7
int_aba_coexp<-int[int$data_source%in%c("ABA_CA1",
"ABA_CA2", "ABA_CA3", "ABA_CA4", "ABA_Sptn", "ABA_subiculum","ABA_DG"),]
int_aba_coexp_pos<-int_aba_coexp[int_aba_coexp$score>=0.7,]
int_aba_coexp_neg<-int_aba_coexp[int_aba_coexp$score<=-0.7,]
int_aba_coexp <- rbind(int_aba_coexp_pos, int_aba_coexp_neg)

# Combine filetered dataset
int<- int[!int$data_source%in%c("ABA_CA1",
"ABA_CA2", "ABA_CA3", "ABA_CA4", "ABA_Sptn", "ABA_subiculum","ABA_DG"),]
int_small_net <- rbind(int,int_aba_coexp)

save(int_small_net, file="case_study/datasets/genes_data/int_small_net.RData")
# 2959910

############### Prepare attributes ###############################################################
print("Loading node attributes")
load(file="case_study/data/node_attributes.RData")
ds_ensg<-unique(c(int_small_net$ensg.A, int_small_net$ensg.B))
node_attributes_small_net<-node_attributes[node_attributes$ensg%in%ds_ensg,]

save(node_attributes_small_net, file="case_study/datasets/genes_data/node_attributes_small_net.RData")
#24825

############### Extracting 1-hop  and 2-hop neighbourhoods of Alzheimer's disease related genes in HENA"############
print("Extracting 1-hop  and 2-hop neighbourhoods of Alzheimer's disease related genes in HENA")
library(igraph)

### Make a graph from interactions dataframe
g <- graph_from_data_frame(int_small_net, directed=FALSE)

### Create a vector of node names that are associated with Alzheimer's disease
node_attributes_small_net$GWAS_pvalue<-as.numeric(node_attributes_small_net$GWAS_pvalue)
gwas_node_att <- node_attributes_small_net[!is.na(node_attributes_small_net$GWAS_pvalue),]
gwas_nodes <- unique(gwas_node_att$ensg)

# List of nodes from IntAct Alzheimer's dataset
int_alz <- int_small_net[int_small_net$data_source%in%"ADIA",]
intact_alz_nodes <-unique(c(int_alz$ensg.A, int_alz$ensg.B))

# Combine node lists
alz_nodes <- unique(c(gwas_nodes, intact_alz_nodes))

### Extract 2HOP neighbourhoods of genes associated with disease
#graph_2hop_alz_subgs <- make_ego_graph(g, 2, alz_nodes)

# Combine individual subgraphs
#g_combined_2hop <- graph_2hop_alz_subgs[[1]]
#for (i in 2:length(graph_2hop_alz_subgs)){
#    g_combined_2hop<-union(g_combined_2hop, graph_2hop_alz_subgs[[i]])
#}
#graph_2hop_alz <- g_combined_2hop

# Check the number of interactions
#gsize(graph_2hop_alz)

# Convert to dataframe from graph structure
#2hop_alz_int <- as_data_frame(graph_2hop_alz, what="edges")

# Save to file as RData and as txt
#save(2hop_alz_int, file="case_study/datasets/genes_data/2hop_alz_int.RData")
#write.table(2hop_alz_int, file="case_study/datasets/genes_data/2hop_alz_int.txt", sep="\t", quote=F, row.names = F)

# Extract and save corresponding node attributes
#2hop_nodes <-unique(c(2hop_alz_int$ensg.A, 2hop_alz_int$ensg.B))
#2hop_alz_attributes <- node_attributes_small_net[node_attributes_small_net$ensg%in%2hop_nodes,]

#save(2hop_alz_attributes, file="case_study/datasets/genes_data/2hop_alz_attributes.RData")
#write.table(2hop_alz_attributes, file="case_study/datasets/genes_data/2hop_alz_attributes.txt", sep="\t", quote=F, row.names = F)


### Extract 1HOP neighbourhoods of genes associated with disease
graph_1hop_alz_subgs <- make_ego_graph(g, 1, alz_nodes)

# Combine individual subgraphs
g_combined_1hop <- as_data_frame(graph_1hop_alz_subgs[[1]], what="edges")
g_combined_1hop<- g_combined_1hop[!duplicated(g_combined_1hop), ]
for (i in 2:length(graph_1hop_alz_subgs)){
    df_add<- as_data_frame(graph_1hop_alz_subgs[[i]] , what="edges")
    df_add <- df_add[!duplicated(df_add), ]
    g_combined_1hop<-rbind(g_combined_1hop, df_add)
    g_combined_1hop <- g_combined_1hop[!duplicated(g_combined_1hop), ]
}
colnames(g_combined_1hop)[1:2]<-c("ensg1", "ensg2")
g_combined_1hop_no_iah<-g_combined_1hop[!g_combined_1hop$data_source%in%"IAH",]

# Select interactions with max values for PPI from IntAct
g_combined_1hop_max_iah <- aggregate(score ~ ensg1 + ensg2 + interaction_type + data_source, data = g_combined_1hop[g_combined_1hop$data_source%in%"IAH",], max)

graph_1hop_alz <- rbind(g_combined_1hop_no_iah,g_combined_1hop_max_iah)

# Save to file as RData and as txt
save(alz_int_1hop, file="case_study/datasets/genes_data/alz_int_1hop.RData")
write.table(alz_int_1hop, file="case_study/datasets/genes_data/alz_int_1hop.txt", sep="\t", quote=F, row.names = F)

# Extract and save corresponding node attributes
nodes_1hop <-unique(c(alz_int_1hop$ensg.A, alz_int_1hop$ensg.B))
alz_attributes_1hop <- node_attributes_small_net[node_attributes_small_net$ensg%in%nodes_1hop,]

save(alz_attributes_1hop, file="case_study/datasets/genes_data/alz_attributes_1hop.RData")
write.table(alz_attributes_1hop, file="case_study/datasets/genes_data/alz_attributes_1hop.txt", sep="\t", quote=F, row.names = F)
###
###
as_data_frame(g, what="edges")
