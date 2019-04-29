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

# Set coexpression cutoff in ABA datasets to >=0.8 and <=-0.8
int_aba_coexp<-int[int$data_source%in%c("ABA_CA1",
"ABA_CA2", "ABA_CA3", "ABA_CA4", "ABA_Sptn", "ABA_subiculum","ABA_DG"),]
int_aba_coexp_pos<-int_aba_coexp[int_aba_coexp$score>=0.8,]
int_aba_coexp_neg<-int_aba_coexp[int_aba_coexp$score<=-0.8,]
int_aba_coexp <- rbind(int_aba_coexp_pos, int_aba_coexp_neg)

# Combine filetered dataset
int<- int[!int$data_source%in%c("ABA_CA1",
"ABA_CA2", "ABA_CA3", "ABA_CA4", "ABA_Sptn", "ABA_subiculum","ABA_DG"),]
int_small_net <- rbind(int,int_aba_coexp)  # 541282

# Write to file
print("Writing filtered dataset to file")
save(int_small_net, file="case_study/datasets/genes_data/int_small_net.RData")
write.table(int_small_net, file="case_study/datasets/genes_data/int_small_nettxt", sep="\t", quote=F, row.names = F)

############### Prepare attributes ###############################################################
print("Loading node attributes")
load(file="case_study/data/node_attributes.RData")

print("Selecting node attributes")
ds_ensg<-unique(c(int_small_net$ensg.A, int_small_net$ensg.B))
node_attributes_small_net<-node_attributes[node_attributes$ensg%in%ds_ensg,]#24768

# Write to file
print("Writing node attributes to file")
save(node_attributes_small_net, file="case_study/datasets/genes_data/node_attributes_small_net.RData")
write.table(node_attributes_small_net, file="case_study/datasets/genes_data/node_attributes_small_net.txt", sep="\t", quote=F, row.names = F)


############### Extracting 1-hop  and 2-hop neighbourhoods of Alzheimer's disease related genes in HENA"############
#print("Extracting 1-hop  and 2-hop neighbourhoods of Alzheimer's disease related genes in HENA")
#library(igraph)

### Make a graph from interactions dataframe
#g <- graph_from_data_frame(int_small_net, directed=FALSE)

### Create a vector of node names that are associated with Alzheimer's disease
#node_attributes_small_net$GWAS_pvalue<-as.numeric(node_attributes_small_net$GWAS_pvalue)
#gwas_node_att <- node_attributes_small_net[!is.na(node_attributes_small_net$GWAS_pvalue),]
#gwas_nodes <- unique(gwas_node_att$ensg)

# List of nodes from IntAct Alzheimer's dataset
#int_alz <- int_small_net[int_small_net$data_source%in%"ADIA",]
#intact_alz_nodes <-unique(c(int_alz$ensg.A, int_alz$ensg.B))

# Combine node lists
#alz_nodes <- unique(c(gwas_nodes, intact_alz_nodes))

### Extract 1HOP neighbourhoods of genes associated with disease
#graph_1hop_alz_subgs <- make_ego_graph(g, 1, alz_nodes)

# Combine individual subgraphs
#g_combined_1hop <- as_data_frame(graph_1hop_alz_subgs[[1]], what="edges")
#g_combined_1hop <- g_combined_1hop[!duplicated(g_combined_1hop), ]
#for (i in 2:length(graph_1hop_alz_subgs)){
#    df_add<- as_data_frame(graph_1hop_alz_subgs[[i]] , what="edges")
#    df_add <- df_add[!duplicated(df_add), ]
#    g_combined_1hop<-rbind(g_combined_1hop, df_add)
#    g_combined_1hop <- g_combined_1hop[!duplicated(g_combined_1hop), ]
#}
#colnames(g_combined_1hop)[1:2]<-c("ensg1", "ensg2")
#g_combined_1hop_no_iah<-g_combined_1hop[!g_combined_1hop$data_source%in%"IAH",]

# Select interactions with max values for PPI from IntAct
#g_combined_1hop_max_iah <- aggregate(score ~ ensg1 + ensg2 + interaction_type + data_source, data = g_combined_1hop[g_combined_1hop$data_source%in%"IAH",], max)

#graph_1hop_alz <- rbind(g_combined_1hop_no_iah,g_combined_1hop_max_iah)
#alz_int_1hop<-graph_1hop_alz
# Save to file as RData and as txt
#save(alz_int_1hop, file="case_study/datasets/genes_data/alz_int_1hop.RData")
#write.table(alz_int_1hop, file="case_study/datasets/genes_data/alz_int_1hop.txt", sep="\t", quote=F, row.names = F)

# Extract and save corresponding node attributes
#nodes_1hop <-unique(c(alz_int_1hop$ensg1, alz_int_1hop$ensg2))
#alz_attributes_1hop <- node_attributes_small_net[node_attributes_small_net$ensg%in%nodes_1hop,]

#save(alz_attributes_1hop, file="case_study/datasets/genes_data/alz_attributes_1hop.RData")
#write.table(alz_attributes_1hop, file="case_study/datasets/genes_data/alz_attributes_1hop.txt", sep="\t", quote=F, row.names = F)

####### Select 100 top GWAS and IntAct Alzheimer genes as alznodes ######
### Create a vector of node names that are associated with Alzheimer's disease
#node_attributes_small_net$GWAS_pvalue<-as.numeric(node_attributes_small_net$GWAS_pvalue)
#gwas_node_att <- node_attributes_small_net[!is.na(node_attributes_small_net$GWAS_pvalue),]
#gwas_node_att_sign <- gwas_node_att[gwas_node_att$GWAS_pvalue<=0.05,1:5]
#
#gwas_node_att_sign_agg <- aggregate(GWAS_pvalue ~ ensg+gene_name, data = gwas_node_att_sign, min)
#gwas_node_att_sign_agg <- gwas_node_att_sign_agg[order(gwas_node_att_sign_agg$GWAS_pvalue),]

#gwas_nodes_sign <- unique(gwas_node_att_sign$ensg)
#alz_int_1hop[alz_int_1hop$ensg1%in%"ENSG00000136717",]

# Check rerviewers nodes
#rv_nodes<-c("ENSG00000136717","ENSG00000186567", "ENSG00000104853", "ENSG00000120885", "ENSG00000203710", "ENSG00000109920", "ENSG00000119684", "ENSG00000110079", "ENSG00000110077", "ENSG00000110077", "ENSG00000073921", "ENSG00000120899",  "ENSG00000120899", "ENSG00000073008", "ENSG00000130204", "ENSG00000100991")

# Explore HENA
#rvsubset_int_large_1 <-integrated_int[integrated_int$ensg.A%in%rv_nodes,]
#rvsubset_int_large_2 <-integrated_int[integrated_int$ensg.B%in%rv_nodes,]
#rvsubset_int_large<-rbind(rvsubset_int_large_1, rvsubset_int_large_2)
#dim(rvsubset_int_large)
#dim(rvsubset_int_large[rvsubset_int_large$interaction_type%in%"coexpression",])

# explore HENA subset with coexp > 0.7 and < -0.7
#rvsubset_small_net_1 <-int_small_net[int_small_net$ensg.A%in%rv_nodes,]
#rvsubset_small_net_2 <-int_small_net[int_small_net$ensg.B%in%rv_nodes,]
#rvsubset_small_net<-rbind(rvsubset_small_net_1, rvsubset_small_net_2)
#dim(rvsubset_small_net)
#dim(rvsubset_small_net[rvsubset_small_net$interaction_type%in%"coexpression",])

#coexp <-rvsubset_small_net[rvsubset_small_net$interaction_type%in%"coexpression",]
#dim(coexp[coexp$score>0.7,])
#summary(int_aba_coexp[int_aba_coexp$score>0,]) # 38636234
#summary(int_aba_coexp[int_aba_coexp$score<0,]) #24343985

