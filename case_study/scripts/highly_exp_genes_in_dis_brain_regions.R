# This script is merging all 6 datasets based on the probesets names
# Then quantile normalization is applied over combined data.
# Finds genes that are expressed higher in the disease-related brain regions 
# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv
print("Reading in AllenBrain data sets. This might take some time. Be patient.")
maexp_178236545 <- read.csv(file = "data/allenbrain/178236545_ds/MicroarrayExpression.csv", header = F)

# Read ontology
onto_178236545 <- read.csv(file = "data/allenbrain/178236545_ds/MicroarrayExpression.csv", header = F)

# Read in probes names
probes_178236545 <- read.csv("data/allenbrain/178236545_ds/Probes.csv")

# Read in sample annotations
sampleann_178236545 <- read.csv("data/allenbrain/178236545_ds/SampleAnnot.csv")
rownames(maexp_178236545) <- maexp_178236545[,1]
maexp_178236545 <- maexp_178236545[,-1]
colnames(maexp_178236545) <- sampleann_178236545$structure_acronym

# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv
maexp_178238266 <- read.csv(file = "data/allenbrain/178238266_ds/MicroarrayExpression.csv", header = F)

# Read ontology
onto_178238266 <- read.csv(file = "data/allenbrain/178238266_ds/MicroarrayExpression.csv", header = F)

# Read in probes names
probes_178238266 <- read.csv("data/allenbrain/178238266_ds/Probes.csv")

# Read in sample annotations
sampleann_178238266 <- read.csv("data/allenbrain/178238266_ds/SampleAnnot.csv")
rownames(maexp_178238266) <- maexp_178238266[,1]
maexp_178238266 <- maexp_178238266[,-1]
colnames(maexp_178238266) <- sampleann_178238266$structure_acrony 

# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv
maexp_178238316 <- read.csv(file = "data/allenbrain/178238316_ds/MicroarrayExpression.csv", header = F)

# Read ontology
onto_178238316 <- read.csv(file = "data/allenbrain/178238316_ds/MicroarrayExpression.csv", header = F)

# Read in probes names
probes_178238316 <- read.csv("data/allenbrain/178238316_ds/Probes.csv")

# Read in sample annotations
sampleann_178238316 <- read.csv("data/allenbrain/178238316_ds/SampleAnnot.csv")
rownames(maexp_178238316) <- maexp_178238316[,1]
maexp_178238316 <- maexp_178238316[,-1]
colnames(maexp_178238316) <- sampleann_178238316$structure_acronym

# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv
maexp_178238359 <- read.csv(file = "data/allenbrain/178238359_ds/MicroarrayExpression.csv", header = F)

# Read ontology
onto_178238359 <- read.csv(file = "data/allenbrain/178238359_ds/MicroarrayExpression.csv", header = F)

# Read in probes names
probes_178238359 <- read.csv("data/allenbrain/178238359_ds/Probes.csv")

# Read in sample annotations
sampleann_178238359 <- read.csv("data/allenbrain/178238359_ds/SampleAnnot.csv")
rownames(maexp_178238359) <- maexp_178238359[,1]
maexp_178238359 <- maexp_178238359[,-1]
colnames(maexp_178238359) <- sampleann_178238359$structure_acronym

# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv
maexp_178238373<- read.csv(file = "data/allenbrain/178238373_ds/MicroarrayExpression.csv", header = F)

# Read ontology
onto_178238373<- read.csv(file = "data/allenbrain/178238373_ds/MicroarrayExpression.csv", header = F)

# Read in probes names
probes_178238373<- read.csv("data/allenbrain/178238373_ds/Probes.csv")

# Read in sample annotations
sampleann_178238373<- read.csv("data/allenbrain/178238373_ds/SampleAnnot.csv")
rownames(maexp_178238373) <- maexp_178238373[,1]
maexp_178238373<- maexp_178238373[,-1]
colnames(maexp_178238373) <- sampleann_178238373$structure_acronym

# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv
maexp_178238387<- read.csv(file = "data/allenbrain/178238387_ds/MicroarrayExpression.csv", header = F)

# Read ontology
onto_178238387<- read.csv(file = "data/allenbrain/178238387_ds/MicroarrayExpression.csv", header = F)

# Read in probes names
probes_178238387<- read.csv("data/allenbrain/178238387_ds/Probes.csv")

# Read in sample annotations
sampleann_178238387<- read.csv("data/allenbrain/178238387_ds/SampleAnnot.csv")
rownames(maexp_178238387) <- maexp_178238387[,1]
maexp_178238387<- maexp_178238387[,-1]
colnames(maexp_178238387) <- sampleann_178238387$structure_acronym

# Merge datasets based on the probesests
print("Merging data sets based on probesets")
# Merge maexp_178236545 and maexp_178238266
ds1_2 <- merge(maexp_178236545[,1:501], maexp_178238266[,1:407], by="row.names", all=T)
row.names(ds1_2) <- ds1_2$Row.names
ds1_2 <- ds1_2[, -1]

# Merge result with maexp_178238316
ds1_2_3 <- merge(ds1_2[,1:908], maexp_178238316[,1:529], by="row.names", all=T)
row.names(ds1_2_3) <- ds1_2_3$Row.names
ds1_2_3 <- ds1_2_3[, -1]

# Merge result with maexp_178238359 
ds1_2_3_4 <- merge(ds1_2_3[,1:1437], maexp_178238359[,1:363], by="row.names", all=T)
row.names(ds1_2_3_4) <- ds1_2_3_4$Row.names
ds1_2_3_4 <- ds1_2_3_4[, -1]

# Merge result with maexp_178238373 
ds1_2_3_4_5 <- merge(ds1_2_3_4[,1:1800], maexp_178238373[,1:893], by="row.names", all=T)
row.names(ds1_2_3_4_5) <- ds1_2_3_4_5$Row.names
ds1_2_3_4_5 <- ds1_2_3_4_5[, -1]

# Merge with maexp_178238387
ds1_2_3_4_5_6 <- merge(ds1_2_3_4_5[,1:2693], maexp_178238387[,1:946], by="row.names", all=T)
row.names(ds1_2_3_4_5_6) <- ds1_2_3_4_5_6$Row.names
ds1_2_3_4_5_6 <- ds1_2_3_4_5_6[, -1]

rm(ds1_2_3_4_5, ds1_2_3_4, ds1_2_3, ds1_2, maexp_178238387, maexp_178238373,maexp_178238359, maexp_178238316, maexp_178238266, maexp_178236545)

# Convert expression to log2
ds1_2_3_4_5_6_log <- log2(ds1_2_3_4_5_6)

# Control batch effect in the merged data
library(preprocessCore)

print("Applying normalization")
norm_brain_mtrx <- normalize.quantiles(as.matrix(ds1_2_3_4_5_6_log))
colnames(norm_brain_mtrx) <- colnames(ds1_2_3_4_5_6_log)
norm_brain_mtrx <- cbind(probe_id=row.names(ds1_2_3_4_5_6_log), norm_brain_mtrx)

# Remove .x, .y, etc from the column names.
colnames(norm_brain_mtrx) <- gsub("\\..*", "", colnames(norm_brain_mtrx))
rownames(norm_brain_mtrx)<-norm_brain_mtrx[,1]
norm_brain_mtrx <- apply(norm_brain_mtrx,2, as.numeric)
rownames(norm_brain_mtrx)<-row.names(ds1_2_3_4_5_6_log)

########### Wilcoxon test #########################################################

p.val=c()
selcols <- c("CA1", "CA2", "CA3", "CA4", "SptN", "DG", "S")

data<-norm_brain_mtrx
for(i in 1: nrow(data)){
    u <- wilcox.test(unlist(data[i,c(colnames(data)%in%selcols)]),unlist(data[i,-c(colnames(data)%in%selcols)]), alternative = "greater")$p.value
    p.val <- c(p.val, u)
}
p.val.adj <- p.adjust(p.val, method = "fdr")
genes <- cbind(rownames(data), p.val.adj)
genes <- as.data.frame(genes)
genes[,1] <- as.character(genes[,1])
genes[,2] <- as.numeric(as.character(genes[,2]))
genes <- genes[genes$p.val.adj <= 0.05,]
wlx_genes_allsam <- genes
dim(wlx_genes_allsam) 


#################
#Load converted probe IDs to ENSG
load(file = "case_study/data/p2ensg_all.RData")

p2ensg_orig$entrez_id <- as.character(p2ensg_orig$entrez_id)
wlx_genes_allsam2engs <- merge(wlx_genes_allsam,p2ensg_orig, by.x="V1", by.y="probe_id", all.x=T)
dim(wlx_genes_allsam2engs)

# Aggregate p-values related to the same gene with multiple porbesets
#Select probeset with the most conservative p-value
library(dplyr)

wlx_genes_allsam2engs_maxp <- wlx_genes_allsam2engs %>%
group_by(Target,name) %>%
slice(which.max(p.val.adj))

wlx_genes_allsam2engs_maxp <- as.data.frame(wlx_genes_allsam2engs_maxp)

dim(wlx_genes_allsam2engs_maxp)
wlx_genes<-wlx_genes_allsam2engs_maxp
rm(wlx_genes_allsam2engs_maxp)
save(wlx_genes, file="case_study/datasets/genes_data/wlx_genes.RData")
