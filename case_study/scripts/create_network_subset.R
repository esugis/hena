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
