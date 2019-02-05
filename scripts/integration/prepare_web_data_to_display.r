# This script prepares smaller set of data to be displayed on the web and use in cytoscape
# The difference with the original dataset
# in co-expression from Allen Brain Atlas coexpression values are >|+-0.5|
# IGRI are excluded
# Create the folder where current results will be written
resdir<-paste("results","integration",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Load integrated dataset
print("Loading interactions")
load(file="results/integration/integrated_int.RData")

# Load node attributes
print("Loading node attributes")
load(file="results/integration/node_attributes.RData")

############### Prepare dataset for the visualisation and analysis in cytoscape ###################
print("Preparing dataset for the visualisation and analysis in cytoscape")
print("Preparing intercations")
# Exclude IGRI
int_cyto<- integrated_int[!integrated_int$interaction_type%in%"IGRI",]

# Set coexpression cutoff in ABA datasets to >0.5 and <-0.5
int_aba_coexp<-int_cyto[int_cyto$data_source%in%c("ABA_CA1",
"ABA_CA2", "ABA_CA3", "ABA_CA4", "ABA_SptN", "ABA_subiculum","ABA_DG"),]
int_aba_coexp_pos<-int_aba_coexp[int_aba_coexp$score>=0.5,]
int_aba_coexp_neg<-int_aba_coexp[int_aba_coexp$score<=-0.5,]
int_aba_coexp_filt<-rbind(int_aba_coexp_pos, int_aba_coexp_neg)
rm(int_aba_coexp, int_aba_coexp_pos, int_aba_coexp_neg)

# Combine filetered dataset
int_cyto<- int_cyto[!int_cyto$data_source%in%c("ABA_CA1",
"ABA_CA2", "ABA_CA3", "ABA_CA4", "ABA_SptN", "ABA_subiculum","ABA_DG"),]
int_cyto_filt <- rbind(int_cyto,int_aba_coexp_filt)
rm(int_cyto,int_aba_coexp_filt)

# Save to file
print("Saving filtered interactions to file")
write.table(int_cyto_filt, file="results/integration/integrated_int_cytoscape.txt", sep="\t", quote=F, row.names=F )
save(int_cyto_filt, file="results/integration/integrated_int_cytosccape.RData")

print("Preparing node attributes")
print("Keep ensg id, gene name, SNP id, biotype, positive selection p-value, GWAS p-value, aggregared exression in CA1, CA2, CA3, CA4, SptN, subiculum and DG brain regions related to the disease.")
# Filter attributes. Keep engs name, gene name, SNP id, positive selection p-value, GWAS p-value, aggregared exression in CA1, CA2, CA3, CA4, SptN, subiculum and DG brain regions related to the disease.

node_attributes_filt <- node_attributes[,c("ensg", "gene_name", "biotype","SNP_id","GWAS_pvalue","ps_pvalue","CA1",
"CA2", "CA3", "CA4","DG","SptN", "S" )]

node_attributes_filt$ps_pvalue<-as.numeric(node_attributes_filt$ps_pvalue)
node_attributes_filt$GWAS_pvalue<-as.numeric(node_attributes_filt$GWAS_pvalue)

print("Writing selected node attributes to file")
write.table(node_attributes_filt, file="results/integration/node_attributes_cytoscape.txt", sep="\t", quote=F, row.names=F )
save(node_attributes_filt, file="results/integration/node_attributes_cytosccape.RData")


############### Prepare dataset with attributes for displaying on the web ###################
print("Preparing dataset with attributes for displaying on the web")
# Merge interactions and node attributes based on ensg.A
int_web_A <-merge(int_cyto_filt, node_attributes_filt, by.x="ensg.A", by.y="ensg", all.x=T )
colnames(int_web_A)[6:17] <- c("gene_name.A","biotype.A","SNP_id.A", "GWAS_pvalue.A", "ps_pvalue.A","CA1.A", "CA2.A", "CA3.A", "CA4.A", "DG.A","SptN.A","S.A")

int_web_AB <-merge(int_web_A, node_attributes_filt, by.x="ensg.B", by.y="ensg", all.x=T )
colnames(int_web_AB)[18:29] <- c("gene_name.B","biotype.B","SNP_id.B", "GWAS_pvalue.B", "ps_pvalue.B","CA1.B", "CA2.B", "CA3.B", "CA4.B", "DG.B","SptN.B","S.B")

integrated_int_attr_web <- int_web_AB
integrated_int_attr_web <- integrated_int_attr_web[,c("gene_name.A","gene_name.B","ensg.A","ensg.B","score","interaction_type","data_source","biotype.A","biotype.B","SNP_id.A","SNP_id.B","GWAS_pvalue.A","GWAS_pvalue.B","ps_pvalue.A","ps_pvalue.B","CA1.A","CA1.B","CA2.A","CA2.B","CA3.A","CA3.B","CA4.A","CA4.B","DG.A","DG.B","S.A","S.B","SptN.A","SptN.B")]


# Sort by GWAS p-value in both GWAS columns
integrated_int_attr_web<-integrated_int_attr_web[order(integrated_int_attr_web$GWAS_pvalue.A,integrated_int_attr_web$GWAS_pvalue.B,integrated_int_attr_web$ps_pvalue.A,integrated_int_attr_web$ps_pvalue.B),]

# Save dataset for web app
net  <- integrated_int_attr_web
#character_vars <- lapply(net, class) == "character"
#net[, character_vars] <- lapply(net[, character_vars], as.factor)

write.table(net, file="results/integration/integrated_int_attributes_web_net_2018_11_22.txt", sep="\t", quote=F, row.names=F )
save(net, file="results/integration/integrated_int_attributes_web_net_2018_11_22.RData")




