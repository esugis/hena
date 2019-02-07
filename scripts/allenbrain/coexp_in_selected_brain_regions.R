#######################################################################################################################################
# This script computes co-expression in CA1, CA2, CA3, CA4, DG, SptN and subiculum brain regions according to Allen Brain Atlas ontology.
# Co-expression is computed over 6 whole-brain datasets of gene expression in healthy brains.
#########################################################################################################################################
# Used libaries
library(reshape)
library(foreach)
library(plyr)
library(psych)
library(doParallel)

# Use 4 cores for parallel computations locally
registerDoParallel(cores=4)

# Use more cores on the server
#registerDoParallel(cores=30)

#############################################################################################
# This part of the script computes co-expression in CA1 brain region from Allen Brain Atlas.
# Co-expression is computed over 6 whole-brain datasets of gene expression in healthy brains.
#############################################################################################
tissue_path <- "results/allenbrain/tissues_rdata_selected/CA1/"

# Load CA1 dataset
print("Loading data")
load(file="results/allenbrain/tissues_rdata_selected/CA1/CA1_preprocessed.RData")

print("Performing co-expression analysis in CA1 region")
tissue <- "CA1"
filenametissue<- sprintf("%s.txt",tissue)
m<-t(m)

probesets <- colnames(m)
foreach(i = 1:length(probesets)) %dopar%{
    cor1gds <- c() # Corelation for one gene in one ds
    probeset <- probesets[i]
    vect <- m[,i]
    mtrx <- m[,-c(i)]
    cor1gds <- t(cor(vect, mtrx, method = "spearman", use = "pairwise.complete.obs"))
    # Compute p-values for correlation matrix
    cor_pval <- corr.p(cor1gds,length(vect),adjust="fdr",alpha=.05)$p
    # Combine probesets' names, corresponding correlations and p-values.
    cor_one_probeset <- cbind(rep(paste(probeset),length(rownames(cor1gds))),rownames(cor1gds),cor1gds,cor_pval)
    colnames(cor_one_probeset) <- c("probeset.A", "probeset.B", "rho", "pvalue")
    # Filter our rows with p-values > 0.05
    cor_one_probeset <- data.frame(cor_one_probeset, stringsAsFactors = F)
    cor_one_probeset$rho <- as.numeric(cor_one_probeset$rho)
    cor_one_probeset$pvalue <- as.numeric(cor_one_probeset$pvalue)
    cor_one_probeset <- cor_one_probeset[cor_one_probeset$pvalue<=0.05,]
    # Append write to txt file
    write.table(cor_one_probeset, file=file.path(tissue_path,filenametissue), append = T, quote=F,sep="\t", row.names=F, col.names=F)
}

# Read the computed correlations and remove duplicated pares AB BA
tmp <- read.table(file=file.path(tissue_path,filenametissue), header=F, sep = "\t")
tmp_orig<-tmp

# p-value <=0.01
tmp_01 <- tmp[tmp$V4<=0.01,]
fst_quartl_pos_01 <- unname(summary(tmp_01[tmp_01$V3>0,3])[2])
frd_quartl_neg_01 <- unname(summary(tmp_01[tmp_01$V3<0,3])[5])
tmp_01_pos <- tmp_01[tmp_01$V3>fst_quartl_pos_01,]
tmp_01_neg <- tmp_01[tmp_01$V3<frd_quartl_neg_01,]
tmp_01 <- rbind(tmp_01_pos,tmp_01_neg)
tmp<-tmp_01

  # Convert factors to characters
  df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Remove duplicated values like AB BA  
  tmp <- df2string(tmp)
  tmp<- tmp[!duplicated(tmp), ]
  tmp <- tmp[!duplicated(data.frame(t(apply(tmp[1:2], 1, sort)), tmp$V3)),]
 filtered_coexpression <- tmp
colnames(filtered_coexpression) <- c("probeset.A", "probeset.B", "cor", "cor_pval")

# Save to file
filenametissue_filt<- sprintf("%s.txt", paste(tissue,"filt", sep="_"))
#write.table(filtered_coexpression, file=file.path(tissue_path,filenametissue_filt), append = F, quote=F,sep="\t", row.names=F, col.names=F)

# Load converted probesets IDs
load(file = "results/allenbrain/p2ensg_all.RData")

# Add ensg IDs for CA1 co-expression dataset
CA1_filt <- filtered_coexpression
rm(filtered_coexpression)
CA1_filt_ensg1 <- merge(CA1_filt, p2ensg_orig, by.x = "probeset.A", by.y = "probe_id", all = F)
rm(CA1_filt)
CA1_filt_ensg12 <- merge(CA1_filt_ensg1, p2ensg_orig, by.x = "probeset.B", by.y = "probe_id", all = F)
rm(CA1_filt_ensg1)

# Select only ensg1 ensg2 score and add type of interaction and data source
colnames(CA1_filt_ensg12)[c(6,9)] <- c("ensg1", "ensg2")
CA1_filt_ensg12 <- CA1_filt_ensg12[,c("ensg1","ensg2", "cor")]
CA1 <- aggregate(cor ~ ensg1 + ensg2, data = CA1_filt_ensg12, max)
CA1 <- CA1[, c("ensg1", "ensg2", "cor")]
CA1 <- cbind(CA1, interaction_type = "coexpression")
rm(CA1_filt_ensg12 )

# Add data source name Allen Brain Atlas (ABA)
CA1 <- cbind(CA1, data_source = "ABA_CA1")

# Rename the columns
colnames(CA1) <- c("ensg1","ensg2", "score", "interaction_type","data_source")

## Remove the duplicated undirrescted pairs with the same score.
# For example, ENSG1-ENSG2 0.5 and ENSG2-ENSG1 0.5
# Convert factors to characters
df2string<-function(df){
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  df[,3]<-as.numeric(df[,3])
  return (df)}

# Co-expression in CA1
CA1 <- df2string(CA1)
CA1<- CA1[!duplicated(CA1), ]
CA1 <- CA1[!duplicated(data.frame(t(apply(CA1[1:2], 1, sort)), CA1$score)),]

# Size after self-loops and duplicated interactions were removed
CA1_coexp_int <- CA1
rm(CA1)
print("Writng results to file")
save(CA1_coexp_int, file ="results/allenbrain/CA1_coexp_int_updated_01.RData")
rm(list=ls())


#############################################################################################
# This part of the script computes co-expression in CA2 brain region from Allen Brain Atlas.
# Co-expression is computed over 6 whole-brain datasets of gene expression in healthy brains.
#############################################################################################
tissue_path <- "results/allenbrain/tissues_rdata_selected/CA2"

# Load CA2 dataset
print("Loading the data")
load(file= "results/allenbrain/tissues_rdata_selected/CA2/CA2_preprocessed.RData")

print("Performing co-expression analysis in CA2 region")
tissue <- "CA2"
filenametissue<- sprintf("%s.txt",tissue)
m<-t(m)
probesets <- colnames(m)

foreach(i = 1:length(probesets)) %dopar%{
    cor1gds <- c() # Corelation for one gene in one ds
    probeset <- probesets[i]
    vect <- m[,i]
    mtrx <- m[,-c(i)]
    cor1gds <- t(cor(vect, mtrx, method = "spearman", use = "pairwise.complete.obs"))
    # Compute p-values for correlation matrix
    cor_pval <- corr.p(cor1gds,length(vect),adjust="fdr",alpha=.05)$p
    # Combine probesets' names, corresponding correlations and p-values.
    cor_one_probeset <- cbind(rep(paste(probeset),length(rownames(cor1gds))),rownames(cor1gds),cor1gds,cor_pval)
    colnames(cor_one_probeset) <- c("probeset.A", "probeset.B", "rho", "pvalue")
    # Filter our rows with p-values > 0.05
    cor_one_probeset <- data.frame(cor_one_probeset, stringsAsFactors = F)
    cor_one_probeset$rho <- as.numeric(cor_one_probeset$rho)
    cor_one_probeset$pvalue <- as.numeric(cor_one_probeset$pvalue)
    cor_one_probeset <- cor_one_probeset[cor_one_probeset$pvalue<=0.05,]
    # Append write to txt file
    write.table(cor_one_probeset, file=file.path(tissue_path,filenametissue), append = T, quote=F,sep="\t", row.names=F, col.names=F)
}

# Read the computed correlations and remove duplicated pares AB BA
tmp <- read.table(file=file.path(tissue_path,filenametissue), header=F, sep = "\t")

# p-value <=0.01
tmp_01 <- tmp[tmp$V4<=0.01,]
fst_quartl_pos_01 <- unname(summary(tmp_01[tmp_01$V3>0,3])[2])
frd_quartl_neg_01 <- unname(summary(tmp_01[tmp_01$V3<0,3])[5])
tmp_01_pos <- tmp_01[tmp_01$V3>fst_quartl_pos_01,]
tmp_01_neg <- tmp_01[tmp_01$V3<frd_quartl_neg_01,]
tmp_01 <- rbind(tmp_01_pos,tmp_01_neg)

# Convert factors to characters
df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Remove duplicated values like AB BA
tmp <- df2string(tmp)
tmp<- tmp[!duplicated(tmp), ]
tmp <- tmp[!duplicated(data.frame(t(apply(tmp[1:2], 1, sort)), tmp$V3)),]
filtered_coexpression <- tmp
colnames(filtered_coexpression) <- c("probeset.A", "probeset.B", "cor", "cor_pval")

# Save to file
filenametissue_filt<- sprintf("%s.txt", paste(tissue,"filt", sep="_"))
#write.table(filtered_coexpression, file=file.path(tissue_path,filenametissue_filt), append = F, quote=F,sep="\t", row.names=F, col.names=F)

# Load converted probesets IDs
load(file = "results/allenbrain/p2ensg_all.RData")

# Add ensg IDs for CA2 co-expression dataset
CA2_filt <- filtered_coexpression
rm(filtered_coexpression)
CA2_filt_ensg1 <- merge(CA2_filt, p2ensg_orig, by.x = "probeset.A", by.y = "probe_id", all = F)
rm(CA2_filt)
CA2_filt_ensg12 <- merge(CA2_filt_ensg1, p2ensg_orig, by.x = "probeset.B", by.y = "probe_id", all = F)
rm(CA2_filt_ensg1)

# Select only ensg1 ensg2 score and add type of interaction and data source
colnames(CA2_filt_ensg12)[c(6,9)] <- c("ensg1", "ensg2")
CA2_filt_ensg12 <- CA2_filt_ensg12[,c("ensg1","ensg2", "cor")]
CA2 <- aggregate(cor ~ ensg1 + ensg2, data = CA2_filt_ensg12, max)
CA2 <- CA2[, c("ensg1", "ensg2", "cor")]
CA2 <- cbind(CA2, interaction_type = "coexpression")
rm(CA2_filt_ensg12)

# Add data source name Allen Brain Atlas (ABA)
CA2 <- cbind(CA2, data_source = "ABA_CA2")

# Rename the columns
colnames(CA2) <- c("ensg1","ensg2", "score", "interaction_type","data_source")

## Remove the duplicated undirrescted pairs with the same score.
# For example, ENSG1-ENSG2 0.5 and ENSG2-ENSG1 0.5
# Convert factors to characters
df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Co-expression in CA2
CA2 <- df2string(CA2)
CA2<- CA2[!duplicated(CA2), ]
CA2 <- CA2[!duplicated(data.frame(t(apply(CA2[1:2], 1, sort)), CA2$score)),]

# Self-loops and duplicated interactions were removed
CA2_coexp_int <- CA2
rm(CA2)
print("Writng results to file")
save(CA2_coexp_int, file = "results/allenbrain/CA2_coexp_int_updated_01.RData")
rm(list=ls())


#############################################################################################
# This part of the script computes co-expression in CA3 brain region from Allen Brain Atlas.
# Co-expression is computed over 6 whole-brain datasets of gene expression in healthy brains.
#############################################################################################
tissue_path <- "results/allenbrain/tissues_rdata_selected/CA3/"

# Load CA3 dataset
print("Loading the data")
load(file = "results/allenbrain/tissues_rdata_selected/CA3/CA3_preprocessed.RData")

print("Performing co-expression analysis in CA3 region")
tissue <- "CA3"
filenametissue<- sprintf("%s.txt",tissue)
m<-t(m)
probesets <- colnames(m)

foreach(i = 1:length(probesets)) %dopar%{
    cor1gds <- c() # Corelation for one gene in one ds
    probeset <- probesets[i]
    vect <- m[,i]
    mtrx <- m[,-c(i)]
    cor1gds <- t(cor(vect, mtrx, method = "spearman", use = "pairwise.complete.obs"))
    # Compute p-values for correlation matrix
    cor_pval <- corr.p(cor1gds,length(vect),adjust="fdr",alpha=.05)$p
    # Combine probesets' names, corresponding correlations and p-values.
    cor_one_probeset <- cbind(rep(paste(probeset),length(rownames(cor1gds))),rownames(cor1gds),cor1gds,cor_pval)
    colnames(cor_one_probeset) <- c("probeset.A", "probeset.B", "rho", "pvalue")
    # Filter our rows with p-values > 0.05
    cor_one_probeset <- data.frame(cor_one_probeset, stringsAsFactors = F)
    cor_one_probeset$rho <- as.numeric(cor_one_probeset$rho)
    cor_one_probeset$pvalue <- as.numeric(cor_one_probeset$pvalue)
    cor_one_probeset <- cor_one_probeset[cor_one_probeset$pvalue<=0.05,]
    # Append write to txt file
    write.table(cor_one_probeset, file=file.path(tissue_path,filenametissue), append = T, quote=F,sep="\t", row.names=F, col.names=F)
}
# Read the computed correlations and remove duplicated pares AB BA
tmp <- read.table(file=file.path(tissue_path,filenametissue), header=F, sep = "\t")
tmp_orig<-tmp

# p-value <=0.01
tmp_01 <- tmp[tmp$V4<=0.01,]
fst_quartl_pos_01 <- unname(summary(tmp_01[tmp_01$V3>0,3])[2])
frd_quartl_neg_01 <- unname(summary(tmp_01[tmp_01$V3<0,3])[5])
tmp_01_pos <- tmp_01[tmp_01$V3>fst_quartl_pos_01,]
tmp_01_neg <- tmp_01[tmp_01$V3<frd_quartl_neg_01,]
tmp_01 <- rbind(tmp_01_pos,tmp_01_neg)
tmp<-tmp_01

# Convert factors to characters
df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Remove duplicated values like AB BA
tmp <- df2string(tmp)
tmp<- tmp[!duplicated(tmp), ]
dim(tmp)
tmp <- tmp[!duplicated(data.frame(t(apply(tmp[1:2], 1, sort)), tmp$V3)),]
filtered_coexpression <- tmp
colnames(filtered_coexpression) <- c("probeset.A", "probeset.B", "cor", "cor_pval")

# Save to file
filenametissue_filt<- sprintf("%s.txt", paste(tissue,"filt", sep="_"))
#write.table(filtered_coexpression, file=file.path(tissue_path,filenametissue_filt), append = F, quote=F,sep="\t", row.names=F, col.names=F)

# Load converted probesets IDs
load(file = "results/allenbrain/p2ensg_all.RData")

# Add ensg IDs for CA1 co-expression dataset
CA3_filt <- filtered_coexpression
rm(filtered_coexpression)
CA3_filt_ensg1 <- merge(CA3_filt, p2ensg_orig, by.x = "probeset.A", by.y = "probe_id", all = F)
rm(CA3_filt)
CA3_filt_ensg12 <- merge(CA3_filt_ensg1, p2ensg_orig, by.x = "probeset.B", by.y = "probe_id", all = F)
rm(CA3_filt_ensg1)

# Select only ensg1 ensg2 score and add type of interaction and data source
colnames(CA3_filt_ensg12)[c(6,9)] <- c("ensg1", "ensg2")

CA3_filt_ensg12 <- CA3_filt_ensg12[,c("ensg1","ensg2", "cor")]
CA3 <- aggregate(cor ~ ensg1 + ensg2, data = CA3_filt_ensg12, max)
CA3 <- CA3[, c("ensg1", "ensg2", "cor")]
CA3 <- cbind(CA3, interaction_type = "coexpression")
rm(CA3_filt_ensg12 )

# Add data source name Allen Brain Atlas (ABA)
CA3 <- cbind(CA3, data_source = "ABA_CA3")

# Rename the columns
colnames(CA3) <- c("ensg1","ensg2", "score", "interaction_type","data_source")

## Remove the duplicated undirrescted pairs with the same score.
# For example, ENSG1-ENSG2 0.5 and ENSG2-ENSG1 0.5
# Convert factors to characters
df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Co-expression in CA3
CA3 <- df2string(CA3)
CA3<- CA3[!duplicated(CA3), ]
CA3 <- CA3[!duplicated(data.frame(t(apply(CA3[1:2], 1, sort)), CA3$score)),]

# Size after self-loops and duplicated interactions were removed
CA3_coexp_int <- CA3
rm(CA3)
print("Writng results to file")
save(CA3_coexp_int, file ="results/allenbrain/CA3_coexp_int_updated_01.RData")
rm(list=ls())


#############################################################################################
# This part of the script computes co-expression in CA4 brain region from Allen Brain Atlas.
# Co-expression is computed over 6 whole-brain datasets of gene expression in healthy brains.
#############################################################################################
tissue_path <- "results/allenbrain/tissues_rdata_selected/CA4"

# Load CA4 dataset
print("Loading the data")
load(file= "results/allenbrain/tissues_rdata_selected/CA4/CA4_preprocessed.RData")

print("Performing co-expression analysis in CA4 region")
tissue <- "CA4"
filenametissue<- sprintf("%s.txt",tissue)
m<-t(m)
probesets <- colnames(m)

foreach(i = 1:length(probesets)) %dopar%{
    cor1gds <- c() # Corelation for one gene in one ds
    probeset <- probesets[i]
    vect <- m[,i]
    mtrx <- m[,-c(i)]
    cor1gds <- t(cor(vect, mtrx, method = "spearman", use = "pairwise.complete.obs"))
    # Compute p-values for correlation matrix
    cor_pval <- corr.p(cor1gds,length(vect),adjust="fdr",alpha=.05)$p
    # Combine probesets' names, corresponding correlations and p-values.
    cor_one_probeset <- cbind(rep(paste(probeset),length(rownames(cor1gds))),rownames(cor1gds),cor1gds,cor_pval)
    colnames(cor_one_probeset) <- c("probeset.A", "probeset.B", "rho", "pvalue")
    # Filter our rows with p-values > 0.05
    cor_one_probeset <- data.frame(cor_one_probeset, stringsAsFactors = F)
    cor_one_probeset$rho <- as.numeric(cor_one_probeset$rho)
    cor_one_probeset$pvalue <- as.numeric(cor_one_probeset$pvalue)
    cor_one_probeset <- cor_one_probeset[cor_one_probeset$pvalue<=0.05,]
    # Append write to txt file
    write.table(cor_one_probeset, file=file.path(tissue_path,filenametissue), append = T, quote=F,sep="\t", row.names=F, col.names=F)
}

# Read the computed correlations and remove duplicated pares AB BA
tmp <- read.table(file=file.path(tissue_path,filenametissue), header=F, sep = "\t")
tmp_orig<-tmp

# p-value <=0.01
tmp_01 <- tmp[tmp$V4<=0.01,]
fst_quartl_pos_01 <- unname(summary(tmp_01[tmp_01$V3>0,3])[2])
frd_quartl_neg_01 <- unname(summary(tmp_01[tmp_01$V3<0,3])[5])
tmp_01_pos <- tmp_01[tmp_01$V3>fst_quartl_pos_01,]
tmp_01_neg <- tmp_01[tmp_01$V3<frd_quartl_neg_01,]
tmp_01 <- rbind(tmp_01_pos,tmp_01_neg)
tmp<-tmp_01

# Convert factors to characters
df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Remove duplicated values like AB BA
tmp <- df2string(tmp)
tmp<- tmp[!duplicated(tmp), ]
tmp <- tmp[!duplicated(data.frame(t(apply(tmp[1:2], 1, sort)), tmp$V3)),]
filtered_coexpression <- tmp
colnames(filtered_coexpression) <- c("probeset.A", "probeset.B", "cor", "cor_pval")

# Save to file
filenametissue_filt<- sprintf("%s.txt", paste(tissue,"filt", sep="_"))
#write.table(filtered_coexpression, file=file.path(tissue_path,filenametissue_filt), append = F, quote=F,sep="\t", row.names=F, col.names=F)

# Load converted probesets IDs
load(file = "results/allenbrain/p2ensg_all.RData")

# Add ensg IDs for CA4 co-expression dataset
CA4_filt <- filtered_coexpression
rm(filtered_coexpression)
CA4_filt_ensg1 <- merge(CA4_filt, p2ensg_orig, by.x = "probeset.A", by.y = "probe_id", all = F)
rm(CA4_filt)
CA4_filt_ensg12 <- merge(CA4_filt_ensg1, p2ensg_orig, by.x = "probeset.B", by.y = "probe_id", all = F)
rm(CA4_filt_ensg1)

# Select only ensg1 ensg2 score and add type of interaction and data source
colnames(CA4_filt_ensg12)[c(6,9)] <- c("ensg1", "ensg2")
CA4_filt_ensg12 <- CA4_filt_ensg12[,c("ensg1","ensg2", "cor")]
CA4 <- aggregate(cor ~ ensg1 + ensg2, data = CA4_filt_ensg12, max)
CA4 <- CA4[, c("ensg1", "ensg2", "cor")]
CA4 <- cbind(CA4, interaction_type = "coexpression")
rm(CA4_filt_ensg12)

# Add data source name Allen Brain Atlas (ABA)
CA4 <- cbind(CA4, data_source = "ABA_CA4")

# Rename the columns
colnames(CA4) <- c("ensg1","ensg2", "score", "interaction_type","data_source")

## Remove the duplicated undirrescted pairs with the same score.
# For example, ENSG1-ENSG2 0.5 and ENSG2-ENSG1 0.5
# Convert factors to characters
df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Co-expression in CA4
CA4 <- df2string(CA4)
CA4<- CA4[!duplicated(CA4), ]
CA4 <- CA4[!duplicated(data.frame(t(apply(CA4[1:2], 1, sort)), CA4$score)),]

# Self-loops and duplicated interactions were removed
CA4_coexp_int <- CA4
rm(CA4)

print("Writng results to file")
save(CA4_coexp_int, file =  "results/allenbrain/CA4_coexp_int_updated_01.RData")
rm(list=ls())


#############################################################################################
# This part of the script computes co-expression in SptN brain region from Allen Brain Atlas.
# Co-expression is computed over 6 whole-brain datasets of gene expression in healthy brains.
#############################################################################################
tissue_path <- "results/allenbrain/tissues_rdata_selected/SptN"

# Load Sptn dataset
print("Loading the data")
load(file="results/allenbrain/tissues_rdata_selected/SptN/SptN_preprocessed.RData")

print("Performing co-expression analysis in SptN region")
tissue <- "Sptn"
filenametissue<- sprintf("%s.txt",tissue)
m<-t(m)
probesets <- colnames(m)

foreach(i = 1:length(probesets)) %dopar%{
    cor1gds <- c() # Corelation for one gene in one ds
    probeset <- probesets[i]
    vect <- m[,i]
    mtrx <- m[,-c(i)]
    cor1gds <- t(cor(vect, mtrx, method = "spearman", use = "pairwise.complete.obs"))
    # Compute p-values for correlation matrix
    cor_pval <- corr.p(cor1gds,length(vect),adjust="fdr",alpha=.05)$p
    # Combine probesets' names, corresponding correlations and p-values.
    cor_one_probeset <- cbind(rep(paste(probeset),length(rownames(cor1gds))),rownames(cor1gds),cor1gds,cor_pval)
    colnames(cor_one_probeset) <- c("probeset.A", "probeset.B", "rho", "pvalue")
    # Filter our rows with p-values > 0.05
    cor_one_probeset <- data.frame(cor_one_probeset, stringsAsFactors = F)
    cor_one_probeset$rho <- as.numeric(cor_one_probeset$rho)
    cor_one_probeset$pvalue <- as.numeric(cor_one_probeset$pvalue)
    cor_one_probeset <- cor_one_probeset[cor_one_probeset$pvalue<=0.05,]
    # Append write to txt file
    write.table(cor_one_probeset, file=file.path(tissue_path,filenametissue), append = T, quote=F,sep="\t", row.names=F, col.names=F)
}

# Read the computed correlations and remove duplicated pares AB BA
tmp <- read.table(file=file.path(tissue_path,filenametissue), header=F, sep = "\t")
tmp_orig<-tmp

# p-value <=0.01
tmp_01 <- tmp[tmp$V4<=0.01,]
fst_quartl_pos_01 <- unname(summary(tmp_01[tmp_01$V3>0,3])[2])
frd_quartl_neg_01 <- unname(summary(tmp_01[tmp_01$V3<0,3])[5])
tmp_01_pos <- tmp_01[tmp_01$V3>fst_quartl_pos_01,]
tmp_01_neg <- tmp_01[tmp_01$V3<frd_quartl_neg_01,]
tmp_01 <- rbind(tmp_01_pos,tmp_01_neg)

tmp<-tmp_01
# Convert factors to characters
df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Remove duplicated values like AB BA
tmp <- df2string(tmp)
tmp<- tmp[!duplicated(tmp), ]
tmp <- tmp[!duplicated(data.frame(t(apply(tmp[1:2], 1, sort)), tmp$V3)),]
filtered_coexpression <- tmp
colnames(filtered_coexpression) <- c("probeset.A", "probeset.B", "cor", "cor_pval")

# Save to file
filenametissue_filt<- sprintf("%s.txt", paste(tissue,"filt", sep="_"))

# Load converted probesets IDs
load(file = "results/allenbrain/p2ensg_all.RData")

# Add ensg IDs for Sptn co-expression dataset
Sptn_filt <- filtered_coexpression
rm(filtered_coexpression)
Sptn_filt_ensg1 <- merge(Sptn_filt, p2ensg_orig, by.x = "probeset.A", by.y = "probe_id", all = F)
rm(Sptn_filt)
Sptn_filt_ensg12 <- merge(Sptn_filt_ensg1, p2ensg_orig, by.x = "probeset.B", by.y = "probe_id", all = F)
rm(Sptn_filt_ensg1)

# Select only ensg1 ensg2 score and add type of interaction and data source
colnames(Sptn_filt_ensg12)[c(6,9)] <- c("ensg1", "ensg2")
Sptn_filt_ensg12 <- Sptn_filt_ensg12[,c("ensg1","ensg2", "cor")]
Sptn <- aggregate(cor ~ ensg1 + ensg2, data = Sptn_filt_ensg12, max)
Sptn <- Sptn[, c("ensg1", "ensg2", "cor")]
Sptn <- cbind(Sptn, interaction_type = "coexpression")

# Add data source name Allen Brain Atlas (ABA)
Sptn <- cbind(Sptn, data_source = "ABA_Sptn")

# Rename the columns
colnames(Sptn) <- c("ensg1","ensg2", "score", "interaction_type","data_source")

## Remove the duplicated undirrescted pairs with the same score.
# For example, ENSG1-ENSG2 0.5 and ENSG2-ENSG1 0.5
# Convert factors to characters
df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Co-expression in Sptn
Sptn <- df2string(Sptn)
Sptn<- Sptn[!duplicated(Sptn), ]
Sptn <- Sptn[!duplicated(data.frame(t(apply(Sptn[1:2], 1, sort)), Sptn$score)),]

# Self-loops and duplicated interactions were removed
Sptn_coexp_int <- Sptn
rm(Sptn)
print("Writng results to file")
save(Sptn_coexp_int, file = "results/allenbrain/Sptn_coexp_int_updated_01.RData")
rm(list=ls())


#############################################################################################
# This part of the script computes co-expression in subiculum brain region from Allen Brain Atlas.
# Co-expression is computed over 6 whole-brain datasets of gene expression in healthy brains.
#############################################################################################
tissue_path <- "results/allenbrain/tissues_rdata_selected/S"

# Load subiculum dataset
print("Loading the data")
load(file= "results/allenbrain/tissues_rdata_selected/S/subiculum_preprocessed.RData")

print("Performning co-expression analysis in subiculum region")
tissue <- "subiculum"
filenametissue<- sprintf("%s.txt",tissue)
m<-t(m)
probesets <- colnames(m)

foreach(i = 1:length(probesets)) %dopar%{
    cor1gds <- c() # Corelation for one gene in one ds
    probeset <- probesets[i]
    vect <- m[,i]
    mtrx <- m[,-c(i)]
    cor1gds <- t(cor(vect, mtrx, method = "spearman", use = "pairwise.complete.obs"))
    # Compute p-values for correlation matrix
    cor_pval <- corr.p(cor1gds,length(vect),adjust="fdr",alpha=.05)$p
    # Combine probesets' names, corresponding correlations and p-values.
    cor_one_probeset <- cbind(rep(paste(probeset),length(rownames(cor1gds))),rownames(cor1gds),cor1gds,cor_pval)
    colnames(cor_one_probeset) <- c("probeset.A", "probeset.B", "rho", "pvalue")
    # Filter our rows with p-values > 0.05
    cor_one_probeset <- data.frame(cor_one_probeset, stringsAsFactors = F)
    cor_one_probeset$rho <- as.numeric(cor_one_probeset$rho)
    cor_one_probeset$pvalue <- as.numeric(cor_one_probeset$pvalue)
    cor_one_probeset <- cor_one_probeset[cor_one_probeset$pvalue<=0.05,]
    
    # Append write to txt file
    write.table(cor_one_probeset, file=file.path(tissue_path,filenametissue), append = T, quote=F,sep="\t", row.names=F, col.names=F)
}

tmp <- read.table(file=file.path(tissue_path,filenametissue), header=F, sep = "\t")
tmp_orig<-tmp

# p-value <=0.01
tmp_01 <- tmp[tmp$V4<=0.01,]
fst_quartl_pos_01 <- unname(summary(tmp_01[tmp_01$V3>0,3])[2])
frd_quartl_neg_01 <- unname(summary(tmp_01[tmp_01$V3<0,3])[5])
tmp_01_pos <- tmp_01[tmp_01$V3>fst_quartl_pos_01,]
tmp_01_neg <- tmp_01[tmp_01$V3<frd_quartl_neg_01,]
tmp_01 <- rbind(tmp_01_pos,tmp_01_neg)

tmp<-tmp_01
# Read the computed correlations and remove duplicated pares AB BA

# Convert factors to characters
df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Remove duplicated values like AB BA
tmp <- df2string(tmp)
tmp<- tmp[!duplicated(tmp), ]
tmp <- tmp[!duplicated(data.frame(t(apply(tmp[1:2], 1, sort)), tmp$V3)),]
filtered_coexpression <- tmp
colnames(filtered_coexpression) <- c("probeset.A", "probeset.B", "cor", "cor_pval")

# Save to file
filenametissue_filt<- sprintf("%s.txt", paste(tissue,"filt", sep="_"))

# Load converted probesets IDs
load(file = "results/allenbrain/p2ensg_all.RData")

# Add ensg IDs for subiculum co-expression dataset
subiculum_filt <- filtered_coexpression
rm(filtered_coexpression)
subiculum_filt_ensg1 <- merge(subiculum_filt, p2ensg_orig, by.x = "probeset.A", by.y = "probe_id", all = F)
rm(subiculum_filt)
subiculum_filt_ensg12 <- merge(subiculum_filt_ensg1, p2ensg_orig, by.x = "probeset.B", by.y = "probe_id", all = F)
rm(subiculum_filt_ensg1)

# Select only ensg1 ensg2 score and add type of interaction and data source
colnames(subiculum_filt_ensg12)[c(6,9)] <- c("ensg1", "ensg2")
subiculum_filt_ensg12 <- subiculum_filt_ensg12[,c("ensg1","ensg2", "cor")]
subiculum <- aggregate(cor ~ ensg1 + ensg2, data = subiculum_filt_ensg12, max)
subiculum <- subiculum[, c("ensg1", "ensg2", "cor")]
subiculum <- cbind(subiculum, interaction_type = "coexpression")
rm(subiculum_filt_ensg12)

# Add data source name Allen Brain Atlas (ABA)
subiculum <- cbind(subiculum, data_source = "ABA_subiculum")

# Rename the columns
colnames(subiculum) <- c("ensg1","ensg2", "score", "interaction_type","data_source")

## Remove the duplicated undirrescted pairs with the same score.
# For example, ENSG1-ENSG2 0.5 and ENSG2-ENSG1 0.5
# Convert factors to characters
df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Co-expression in subiculum
subiculum <- df2string(subiculum)
subiculum<- subiculum[!duplicated(subiculum), ]
subiculum <- subiculum[!duplicated(data.frame(t(apply(subiculum[1:2], 1, sort)), subiculum$score)),]

# Self-loops and duplicated interactions were removed
subiculum_coexp_int <- subiculum
rm(subiculum)
print("Writng results to file")
save(subiculum_coexp_int, file = "results/allenbrain/subiculum_coexp_int_updated_01.RData")
rm(list=ls())


#############################################################################################
# This part of the script computes co-expression in DG brain region from Allen Brain Atlas.
# Co-expression is computed over 6 whole-brain datasets of gene expression in healthy brains.
#############################################################################################
tissue_path <- "results/allenbrain/tissues_rdata_selected/DG"
# Load DG dataset
print("Loading the data")
#load(file= "results/allenbrain/tissues_rdata_selected/DG/DG_preprocessed.RData")

print("Performing co-expression analysis in DG region")
##### Try with foreach per one probeset
tissue <- "DG"
filenametissue<- sprintf("%s.txt",tissue)

# Read the computed correlations and remove duplicated pares AB BA
tmp <- read.table(file=file.path(tissue_path,filenametissue), header=F, sep = "\t")
dim(tmp)
tmp_orig<-tmp

# p-value <=0.01
tmp_01 <- tmp[tmp$V4<=0.01,]
fst_quartl_pos_01 <- unname(summary(tmp_01[tmp_01$V3>0,3])[2])
frd_quartl_neg_01 <- unname(summary(tmp_01[tmp_01$V3<0,3])[5])
if(fst_quartl_pos_01>=0.5){
    tmp_01_pos <- tmp_01[tmp_01$V3>fst_quartl_pos_01,]}else{tmp_01_pos <- tmp_01[tmp_01$V3>0,]}
if(frd_quartl_neg_01<=-0.5){
    tmp_01_neg <- tmp_01[tmp_01$V3<frd_quartl_neg_01,]}else{tmp_01_neg <- tmp_01[tmp_01$V3<0,]}

tmp_01 <- rbind(tmp_01_pos,tmp_01_neg)
dim(tmp_01)

tmp<-tmp_01

# Convert factors to characters
df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Remove duplicated values like AB BA
tmp <- df2string(tmp)
tmp<- tmp[!duplicated(tmp), ]
tmp <- tmp[!duplicated(data.frame(t(apply(tmp[1:2], 1, sort)), tmp$V3)),]
filtered_coexpression <- tmp
colnames(filtered_coexpression) <- c("probeset.A", "probeset.B", "cor", "cor_pval")

# Save to file
filenametissue_filt<- sprintf("%s.txt", paste(tissue,"filt", sep="_"))
#write.table(filtered_coexpression, file=file.path(tissue_path,filenametissue_filt), append = F, quote=F,sep="\t", row.names=F, col.names=F)

# Load converted probesets IDs
load(file = "results/allenbrain/p2ensg_all.RData")

# Add ensg IDs for DG co-expression dataset
DG_filt <- filtered_coexpression
rm(filtered_coexpression)
DG_filt_ensg1 <- merge(DG_filt, p2ensg_orig, by.x = "probeset.A", by.y = "probe_id", all = F)
rm(DG_filt)
DG_filt_ensg12 <- merge(DG_filt_ensg1, p2ensg_orig, by.x = "probeset.B", by.y = "probe_id", all = F)
rm(DG_filt_ensg1)

# Select only ensg1 ensg2 score and add type of interaction and data source
colnames(DG_filt_ensg12)[c(6,9)] <- c("ensg1", "ensg2")

DG_filt_ensg12 <- DG_filt_ensg12[,c("ensg1","ensg2", "cor")]
DG <- aggregate(cor ~ ensg1 + ensg2, data = DG_filt_ensg12, max)
DG <- DG[, c("ensg1", "ensg2", "cor")]
DG <- cbind(DG, interaction_type = "coexpression")

# Add data source name Allen Brain Atlas (ABA)
DG <- cbind(DG, data_source = "ABA_DG")
rm(DG_filt_ensg12)

# Rename the columns
colnames(DG) <- c("ensg1","ensg2", "score", "interaction_type","data_source")

## Remove the duplicated undirrescted pairs with the same score.
# For example, ENSG1-ENSG2 0.5 and ENSG2-ENSG1 0.5
# Convert factors to characters
df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Co-expression in DG
DG <- df2string(DG)
DG<- DG[!duplicated(DG), ]
DG <- DG[!duplicated(data.frame(t(apply(DG[1:2], 1, sort)), DG$score)),]

# Self-loops and duplicated interactions were removed
DG_coexp_int <- DG
rm(DG)
print("Writng results to file")
save(DG_coexp_int, file = "results/allenbrain/DG_coexp_int_updated_01.RData")
rm(list=ls())

#############################################################################################
# Additional co-expression filtering based on cutoff values
#############################################################################################

# Load individual co-expression datasets
load(file = "results/allenbrain/CA1_coexp_int_updated_01.RData")
load(file = "results/allenbrain/CA2_coexp_int_updated_01.RData")
load(file = "results/allenbrain/CA3_coexp_int_updated_01.RData")
load(file = "results/allenbrain/CA4_coexp_int_updated_01.RData")
load(file = "results/allenbrain/Sptn_coexp_int_updated_01.RData")
load(file = "results/allenbrain/subiculum_coexp_int_updated_01.RData")
load(file = "results/allenbrain/DG_coexp_int_updated_01.RData")

# Set cut-off value for positive correlation coefficient
pos<-c(CA1_coexp_int[CA1_coexp_int$score>0,3],
CA2_coexp_int[CA2_coexp_int$score>0,3],
CA3_coexp_int[CA3_coexp_int$score>0,3],
CA4_coexp_int[CA4_coexp_int$score>0,3],
Sptn_coexp_int[Sptn_coexp_int$score>0,3],
subiculum_coexp_int[subiculum_coexp_int$score>0,3],
DG_coexp_int[DG_coexp_int$score>0,3])
summary(pos)
#cutoff_pos <- 0.5
# Set cut-off based on the median of the positive correlation coeffcient
cutoff_pos<-unname(summary(pos))[3]

# Set cut-off value for negative correlation coefficient
# Set cut-off based on the median of the negative correlation coeffcient

neg<-c(CA1_coexp_int[CA1_coexp_int$score<0,3],
CA2_coexp_int[CA2_coexp_int$score<0,3],
CA3_coexp_int[CA3_coexp_int$score<0,3],
CA4_coexp_int[CA4_coexp_int$score<0,3],
Sptn_coexp_int[Sptn_coexp_int$score<0,3],
subiculum_coexp_int[subiculum_coexp_int$score<0,3],
DG_coexp_int[DG_coexp_int$score<0,3])
cutoff_neg<- unname(summary(neg))[3]

# Filter based on  positive and negative cutoff values
CA1_coexp_int_pos<-CA1_coexp_int[CA1_coexp_int$score>cutoff_pos,]
CA1_coexp_int_neg<-CA1_coexp_int[CA1_coexp_int$score<cutoff_neg,]
CA1_coexp_int<-rbind(CA1_coexp_int_pos, CA1_coexp_int_neg)

CA2_coexp_int_pos<-CA2_coexp_int[CA2_coexp_int$score>cutoff_pos,]
CA2_coexp_int_neg<-CA2_coexp_int[CA2_coexp_int$score<cutoff_neg,]
CA2_coexp_int<-rbind(CA2_coexp_int_pos, CA2_coexp_int_neg)

CA3_coexp_int_pos<-CA3_coexp_int[CA3_coexp_int$score>cutoff_pos,]
CA3_coexp_int_neg<-CA3_coexp_int[CA3_coexp_int$score<cutoff_neg,]
CA3_coexp_int<-rbind(CA3_coexp_int_pos, CA3_coexp_int_neg)

CA4_coexp_int_pos<-CA4_coexp_int[CA4_coexp_int$score>cutoff_pos,]
CA4_coexp_int_neg<-CA4_coexp_int[CA4_coexp_int$score<cutoff_neg,]
CA4_coexp_int<-rbind(CA4_coexp_int_pos, CA4_coexp_int_neg)

Sptn_coexp_int_pos<-Sptn_coexp_int[Sptn_coexp_int$score>cutoff_pos,]
SptN_coexp_int_pos<-Sptn_coexp_int_pos
Sptn_coexp_int_neg<-Sptn_coexp_int[Sptn_coexp_int$score<cutoff_neg,]
SptN_coexp_int_neg<-Sptn_coexp_int_neg
SptN_coexp_int<-rbind(SptN_coexp_int_pos, SptN_coexp_int_neg)

subiculum_coexp_int_pos<-subiculum_coexp_int[subiculum_coexp_int$score>cutoff_pos,]
subiculum_coexp_int_neg<-subiculum_coexp_int[subiculum_coexp_int$score<cutoff_neg,]
subiculum_coexp_int<-rbind(subiculum_coexp_int_pos, subiculum_coexp_int_neg)

DG_coexp_int_pos<-DG_coexp_int[DG_coexp_int$score>cutoff_pos,]
DG_coexp_int_neg<-DG_coexp_int[DG_coexp_int$score<cutoff_neg,]
DG_coexp_int<-rbind(DG_coexp_int_pos, DG_coexp_int_neg)


# Write to file as RData
save(CA1_coexp_int, file = "results/allenbrain/CA1_coexp_int_aba.RData")
save(CA2_coexp_int, file = "results/allenbrain/CA2_coexp_int_aba.RData")
save(CA3_coexp_int, file = "results/allenbrain/CA3_coexp_int_aba.RData")
save(CA4_coexp_int, file = "results/allenbrain/CA4_coexp_int_aba.RData")
save(SptN_coexp_int, file = "results/allenbrain/SptN_coexp_int_aba.RData")
save(subiculum_coexp_int, file = "results/allenbrain/subiculum_coexp_int_aba.RData")
save(DG_coexp_int, file = "results/allenbrain/DG_coexp_int_aba.RData")

dim(rbind(CA1_coexp_int,CA2_coexp_int,CA3_coexp_int,CA4_coexp_int, SptN_coexp_int,subiculum_coexp_int, DG_coexp_int))
# Write to txt file
write.table(CA1_coexp_int, file = "results/allenbrain/CA1_coexp_int_aba.txt", quote = F, sep = "\t",row.names = F)
write.table(CA2_coexp_int, file = "results/allenbrain/CA2_coexp_int_aba.txt", quote = F, sep = "\t",row.names = F)
write.table(CA3_coexp_int, file = "results/allenbrain/CA3_coexp_int_aba.txt", quote = F, sep = "\t",row.names = F)
write.table(CA4_coexp_int, file = "results/allenbrain/CA4_coexp_int_aba.txt", quote = F, sep = "\t",row.names = F)
write.table(SptN_coexp_int, file = "results/allenbrain/SptN_coexp_int_aba.txt", quote = F, sep = "\t",row.names = F)
write.table(subiculum_coexp_int, file = "results/allenbrain/subiculum_coexp_int_aba.txt", quote = F, sep = "\t",row.names = F)
write.table(DG_coexp_int, file = "results/allenbrain/DG_coexp_int_aba.txt", quote = F, sep = "\t",row.names = F)

rm(list=ls())




