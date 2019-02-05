# This script will prepare gene lists and find matches in tthe results provided by the classifier
# All individual gene lists are obrained from the publications listed in the articles.
# Gene names are converted to the ENSG IDs, search is performed based on the ENSG ID.

library(gProfileR)

# Load the results from the model
modres <- read.csv(file="case_study/datasets/genes_data/predictions_ranked.csv", sep=';')
modres<-modres[, -1]
# Convert ENSG ids
#moderes_ensg2g <- gconvert(unique(modres$ensg))


############ 1. Genes from MalaCards database
genecards <- read.csv(file="case_study/data/interpretation/genecards/GeneCards_gene_names.csv")
genecards_genes<-as.character(genecards$gene_name)
genecards_g2ensg<- gconvert(genecards_genes)
genecards_g2ensg <- genecards_g2ensg[,c(4,5)]
genecards_g2ensg <- genecards_g2ensg[!duplicated(genecards_g2ensg),]
colnames(genecards_g2ensg) <-c("ensg", "gene_name")
save(genecards_g2ensg, file="case_study/data/interpretation/genecards/genecards_g2ensg.RData")
#Overlap with the genes of interest
genecards_modres <- merge(modres, genecards_g2ensg, by.x = "ensg", by.y = "ensg", all=F)
genecards_modres <- genecards_modres[order(- genecards_modres$hinsage_pred),]
save(genecards_modres,file="case_study/data/interpretation/genecards_modres.RData")


############ 2. Genes from protein arrray data. Disease-Specific Autoantibody Profiles in Human Sera
# Genes extracted from the article Nagele et al. (published)
protorray_genes<-read.csv(file="case_study/data/interpretation/protoarray_alz/genes_protoarray_alz.csv")
protorray_g2ensg <- gconvert(protorray_genes$gene_name)
protorray_g2ensg <-protorray_g2ensg[, c(4,5)]
colnames(protorray_g2ensg) <-c("ensg", "gene_name")
save(protorray_g2ensg, file="case_study/data/interpretation/protoarray_alz/protorray_g2ensg.RData")
# Overlap with the genes of interest
protoarray_modres <- merge(modres, protorray_g2ensg, by.x = "ensg", by.y = "ensg", all.y = F)
save(protoarray_modres,file="case_study/data/interpretation/protoarray_modres.RData")


############ 3. Genes from GWAX  
# Genes extracted from the article Liu et al. (published)
liu<-read.csv(file="case_study/data/interpretation/liu/liu_alzgenes_sup_table4.csv")
liu_g2ensg <- gconvert(liu$gene_name)
liu_g2ensg <-liu_g2ensg[, c(4,5)]
colnames(liu_g2ensg) <-c("ensg", "gene_name")
save(liu_g2ensg, file="case_study/data/interpretation/liu/liu_g2ensg.RData")
# Overlap with the genes of interest
liu_modres <- merge(modres, liu_g2ensg, by.x = "ensg", by.y = "ensg", all.y = F)
save(liu_modres,file="case_study/data/interpretation/liu_modres.RData")


############ 4. Genes from GWAS 
# Genes extracted from the article Jansen et al.  (bioRxiv, not peer-reviewed yet)
jansen <- read.csv(file="case_study/data/interpretation/jansen/jansen_alzgenes_sup_table5.csv", sep=";")
jansen_alzgenes<-unique(as.character(jansen$nearestGene))
jansen_g2ensg<- gconvert(jansen_alzgenes)
jansen_g2ensg <- jansen_g2ensg[,c(4,5)]
jansen_g2ensg <- jansen_g2ensg[!duplicated(jansen_g2ensg),]
colnames(jansen_g2ensg) <-c("ensg", "gene_name")
save(jansen_g2ensg, file="case_study/data/interpretation/jansen/jansen_g2ensg.RData")
# Overlap with the genes of interest
jansen_modres <- merge(modres,jansen_g2ensg, by.x = "ensg", by.y = "ensg", all=F)
save(jansen_modres,file="case_study/data/interpretation/jansen_modres.RData")

############ 5. Genes from GWAS 
# Genes extracted from the article Kunkle et al. (bioRxiv, not peer-reviewed yet)
kunkle <- read.csv(file="case_study/data/interpretation/kunkle/kunkle_alzgenes.csv")
kunkle_alzgenes<-unique(as.character(kunkle$gene_name))
kunkle_g2ensg<- gconvert(kunkle_alzgenes)
kunkle_g2ensg <- kunkle_g2ensg[,c(4,5)]
kunkle_g2ensg <- kunkle_g2ensg[!duplicated(kunkle_g2ensg),]
colnames(kunkle_g2ensg) <-c("ensg", "gene_name")
save(kunkle_g2ensg, file="case_study/data/interpretation/kunkle/kunkle_g2ensg.RData")
# Overlap with the genes of interest
kunkle_modres <- merge(modres,kunkle_g2ensg, by.x = "ensg", by.y = "ensg", all=F)
save(kunkle_modres,file="case_study/data/interpretation/kunkle_modres.RData")


############ 5. Genes from GWAS 
# Genes extracted from the article Marioni et al.(published)
marioni <- read.csv(file="case_study/data/interpretation/marioni/marioni_alzgenes.csv")
marioni_alzgenes<-unique(as.character(marioni$gene_name))
marioni_g2ensg<- gconvert(marioni_alzgenes)
marioni_g2ensg <- marioni_g2ensg[,c(4,5)]
marioni_g2ensg <- marioni_g2ensg[!duplicated(marioni_g2ensg),]
colnames(marioni_g2ensg) <-c("ensg", "gene_name")
save(marioni_g2ensg, file="case_study/data/interpretation/marioni/marioni_g2ensg.RData")
# Overlap with the genes of interest
marioni_modres <- merge(modres,marioni_g2ensg, by.x = "ensg", by.y = "ensg", all=F)
save(marioni_modres,file="case_study/data/interpretation/marioni_modres.RData")


############ 6. Genes from GWAS, TWAS, splicing analysis
# Genes extracted from the article Raj et al. (published)
raj <- read.csv(file="case_study/data/interpretation/raj/raj_alzgenes_table10.csv")
raj_alzgenes<-unique(as.character(raj$gene_name))
raj_g2ensg<- gconvert(raj_alzgenes)
raj_g2ensg <- raj_g2ensg[,c(4,5)]
raj_g2ensg <- raj_g2ensg[!duplicated(raj_g2ensg),]
colnames(raj_g2ensg) <-c("ensg", "gene_name")
save(raj_g2ensg, file="case_study/data/interpretation/raj/raj_g2ensg.RData")
# Overlap with the genes of interest
raj_modres <- merge(modres,raj_g2ensg, by.x = "ensg", by.y = "ensg", all=F)
save(raj_modres,file="case_study/data/interpretation/raj_modres.RData")

# Select genes predicted by both methods that are among Alzheimer's-related independent data sets
all_results <- rbind(genecards_modres , protoarray_modres, raj_modres, kunkle_modres, jansen_modres, marioni_modres, liu_modres )
all_results<- all_results[!duplicated(all_results),]
rf_alz <- all_results[all_results$rf_pred>0.5,]
write.table(rf_alz, file="case_study/datasets/genes_data/rf_predictions_alzheimer.csv", sep=" ", quote=F, row.names = F)
hinsage_alz <- all_results[all_results$hinsage_pred>0.5,]
write.table(hinsage_alz, file="case_study/datasets/genes_data/hinsage_predictions_alzheimer.csv", sep=" ", quote=F, row.names = F)
