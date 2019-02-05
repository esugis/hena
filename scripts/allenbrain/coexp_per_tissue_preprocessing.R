# This script
# - combines 6 individual datasets from Allen Brain Atlas.
# - normalizes the data.
# - outputs subsets of preprocessed data for brain regions related to the Alzheimer's disease.

# Used libaries
library(reshape)
library(foreach)


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
#norm_brain_mtrx <- cbind(probe_id=row.names(ds1_2_3_4_5_6_log), norm_brain_mtrx)
rownames(norm_brain_mtrx) <- row.names(ds1_2_3_4_5_6_log)

# Remove temporary data
#rm(ds1_2_3_4_5_6_log, ds1_2_3_4_5_6)

# Remove .x, .y, etc from the column names.
colnames(norm_brain_mtrx) <- gsub("\\..*", "", colnames(norm_brain_mtrx))
norm_brain_mtrx <- apply(norm_brain_mtrx,2, as.numeric)
rownames(norm_brain_mtrx) <- row.names(ds1_2_3_4_5_6_log)
save(norm_brain_mtrx, file="results/allenbrain/norm_brain_mtrx_selected_regions.RData")

# Filter out probes based on SD
SD <- apply(norm_brain_mtrx, 1, sd, na.rm = T)
norm_brain_mtrx_filt <- norm_brain_mtrx[SD >= 0.29, ]


# Path to the results  
pathRdata <- "results/allenbrain/tissues_rdata_selected"
dir.create(file.path(pathRdata), showWarnings = FALSE, recursive = TRUE)

# A number of tissues relevant to the disease were selected based on the domain expert input
# Allen Brain atlas brain tissue ontology was used for defning the regions and their IDs
# DG, CA1, CA2, CA3, CA4, S (subiculum), SptN (septal nuclei)
# Tissue IDs for left and right part are separated.
# They will be grouped together in individual correlation analysis.

#selected_tissues<-c(4258, 4267,4254, 4263, 4255, 4264, 4256, 4265, 4257, 4266, 4251, 4260, 4301, 4304)

# Subset one tissue
# DG
tissue <- "DG"
print("current tissue")
print(tissue)

# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)
m<- norm_brain_mtrx_filt[, colnames(norm_brain_mtrx_filt)%in%tissue]

# Write data
print("Writing preprocessed DG data file")
file_name <- "DG_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))

# CA1
tissue <- "CA1"
print("current tissue")
print(tissue)

# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)
m<- norm_brain_mtrx_filt[, colnames(norm_brain_mtrx_filt)%in%tissue]

# Write example data
print("Writing preprocessed CA1 data file")
file_name <- "CA1_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))


# CA2
tissue <- "CA2"

# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)
m<- norm_brain_mtrx_filt[, colnames(norm_brain_mtrx_filt)%in%tissue]

# Write example data
print("Writing preprocessed CA2 data file")
file_name <- "CA2_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))

# CA3
tissue <- "CA3"

# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)
m<- norm_brain_mtrx_filt[, colnames(norm_brain_mtrx_filt)%in%tissue]

# Write example data
print("Writing preprocessed CA3 data file")
file_name <- "CA3_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))

# CA4
tissue <- "CA4"

# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)
m<- norm_brain_mtrx_filt[, colnames(norm_brain_mtrx_filt)%in%tissue]

# Write example data
print("Writing preprocessed CA4 data file")
file_name <- "CA4_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))


# subiculum
tissue <- "S"

# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)
m<- norm_brain_mtrx_filt[, colnames(norm_brain_mtrx_filt)%in%tissue]

# Write example data
print("Writing preprocessed subiculum data file")
file_name <- "subiculum_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))

#SptN
tissue <- "SptN"

# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)
m<- norm_brain_mtrx_filt[, colnames(norm_brain_mtrx_filt)%in%tissue]

# Write example data
print("Writing preprocessed SptN data file")
file_name <- "SptN_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))



  



 
