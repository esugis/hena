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
load(file="results/integration/integrated_int.RData")
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


print("Extracting 1-hop neighbourhoods of Alzheimer's disease related genes in HENA")
library(igraph)

### Make a graph from interactions dataframe
g <- graph_from_data_frame(int_small_net, directed=FALSE)

#V(g)$name <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")

### TODO!!!!!Create a vector of node names that are associated with Alzheimer's disease

node_attributes_small_net$GWAS_pvalue<-as.numeric(node_attributes_small_net$GWAS_pvalue)
gwas_node_att <- node_attributes_small_net[!is.na(node_attributes_small_net$GWAS_pvalue),c(2)]
gwas_nodes <- unique(gwas_node_att$ensg)

# List of nodes from IntAct Alzheimer's dataset
int_alz <- int_small_net[int_small_net$data_source%in%"ADIA",]
intact_alz_nodes <-unique(c(int_alz$ensg.A, int_alz$ensg.B))

# Combine node lists
alz_nodes <- c(gwas_nodes, intact_alz_nodes)

### Extract 2HOP neighbourhoods of genes associated with disease
graph_2hop_alz_subgs <- make_ego_graph(g, 2, alz_nodes)

# Combine individual subgraphs
g_combined_2hop <- graph_2hop_alz_subgs[[1]]
for (i in 2:length(graph_2hop_alz_subgs)){
    g_combined_2hop<-union(g_combined_2hop, graph_2hop_alz_subgs[[i]])
}
graph_2hop_alz <- g_combined_2hop

# Check the number of interactions
gsize(graph_2hop_alz)

# Convert to dataframe from graph structure
2hop_alz_int <- as_data_frame(graph_2hop_alz, what="edges")

# Save to file as RData and as txt
save(2hop_alz_int, file="2hop_alz_int.RData")
write.table(2hop_alz_int, file="case_study/datasets/genes_data/2hop_alz_int.txt", sep="\t", quote=F, row.names = F)

### Extract 1HOP neighbourhoods of genes associated with disease
graph_1hop_alz_subgs <- make_ego_graph(g, 1, alz_nodes)

# Combine individual subgraphs
g_combined_1hop <- graph_1hop_alz_subgs[[1]]
for (i in 2:length(graph_1hop_alz_subgs)){
    g_combined_1hop<-union(g_combined_1hop, graph_1hop_alz_subgs[[i]])
}
graph_1hop_alz <- g_combined_1hop

# Check the number of interactions
gsize(graph_1hop_alz)

# Convert to dataframe from graph structure
1hop_alz_int <- as_data_frame(graph_1hop_alz, what="edges")

# Save to file as RData and as txt
save(1hop_alz_int, file="1hop_alz_int.RData")
write.table(1hop_alz_int, file="case_study/datasets/genes_data/1hop_alz_int.txt", sep="\t", quote=F, row.names = F)
