# Converts the affys to ensg. 
# Not all were convertable. See below additional convertion and fitlering out affys that are not mapped to ensg.
# Set working director

# Create the folder where current results will be written
resdir<-"results/adn/all_probes/"
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Used libarries
library(foreach); library(doMC); cores=10 ; registerDoMC(cores);
library(stringr);
library("R.utils")
library(ncdf4);
library(gProfileR)


# Load the probeset names in each of the datasets after("-" and "_" in the original names were substituted with "_")
print("Loading the data")
load(file = "results/adn/all_probes/E_GEOD_18309_all_probes.RData")
load(file = "results/adn/all_probes/E_GEOD_28146_all_probes.RData")
load(file = "results/adn/all_probes/E_GEOD_29652_all_probes.RData")
load(file = "results/adn/all_probes/E_GEOD_4757_all_probes.RData")
load(file = "results/adn/all_probes/E_GEOD_5281_all_probes.RData")
load(file = "results/adn/all_probes/E_MEXP_2280_all_probes.RData")

# Path to the folder where affy2ensg results will be written
path_affy_ensg <- "results/adn/all_probes/"

# Select unique names in all the lists
gene_list <- unique(c(E_GEOD_18309_pr, E_GEOD_28146_pr, E_GEOD_29652_pr, E_GEOD_4757_pr, E_GEOD_5281_pr, E_MEXP_2280_pr))

# Convert gene_list to ENSG
affy2ensg <- gconvert(gene_list)
affy2ensg <- affy2ensg[,c(2,4)]
colnames(affy2ensg)<-c(".id", "Target")

# Remove duplicates
affy2ensg <- affy2ensg[!duplicated(affy2ensg), ]

# Save the result
save(affy2ensg, file = "results/adn/all_probes/affy2ensg.RData")

# Check if the id was included after substituting probe name containig "/" "-" to "_".

# Extract original probesets' names from E-GEOD-18309
E_GEOD_18309<- nc_open("data/adn/E-GEOD-18309.nc");
ncvar_get(E_GEOD_18309, "MetadataOrder")

# Select only Alzheimer related samples based on the additional info in the ArrayExpress
data_E_GEOD_18309 <- ncvar_get(E_GEOD_18309, "data")
rownames(data_E_GEOD_18309) <- ncvar_get(E_GEOD_18309, "MetadataOrder")
colnames(data_E_GEOD_18309) <- ncvar_get(E_GEOD_18309,"gene")
data_E_GEOD_18309 <- t(data_E_GEOD_18309)
nc_close(E_GEOD_18309)
colnames(data_E_GEOD_18309)[1:3] <- "alz"
colnames(data_E_GEOD_18309)[7:9] <- "norm"
data_E_GEOD_18309 <- data_E_GEOD_18309[,c(1:3,7:9)]

# Filter  out rows with SD values less then 0.29
SD <- apply(data_E_GEOD_18309, 1, sd, na.rm = T)
data_E_GEOD_18309_filt <- data_E_GEOD_18309[SD >= 0.29, ]
m <- t(data_E_GEOD_18309_filt)
ds_genes <- colnames(m)
E_GEOD_18309_pr_init <- ds_genes

# Load gsub corrected variant where "-" and "\" are changed to "_" to use probeset name as afile name. E_GEOD_18309_pr
load(file = "results/adn/all_probes/E_GEOD_18309_all_probes.RData")

# Identify how namy probes have names containing "-" and "/"
E_GEOD_18309_no <- E_GEOD_18309_pr_init[!E_GEOD_18309_pr_init%in%E_GEOD_18309_pr]
length(E_GEOD_18309_no)

# Extract original probesets' names from E_GEOD_28146
E_GEOD_28146 <- nc_open("data/adn/E-GEOD-28146.nc");
ncvar_get(E_GEOD_28146, "disease_status")
data_E_GEOD_28146 <- ncvar_get(E_GEOD_28146, "data")
rownames(data_E_GEOD_28146) <- ncvar_get(E_GEOD_28146, "MetadataOrder")
colnames(data_E_GEOD_28146) <- ncvar_get(E_GEOD_28146,"gene")
data_E_GEOD_28146 <- t(data_E_GEOD_28146)
nc_close(E_GEOD_28146)
colnames(data_E_GEOD_28146)[1:22] <- "alz"
colnames(data_E_GEOD_28146)[23:30] <- "norm"

# Filter  out rows with SD values less then 0.29
SD <- apply(data_E_GEOD_28146, 1, sd, na.rm = T)
data_E_GEOD_28146_filt <- data_E_GEOD_28146[SD >= 0.29, ]
m <- t(data_E_GEOD_28146_filt)
ds_genes <- colnames(m)
E_GEOD_28146_pr_init <- ds_genes

# Load gsub corrected variant where "-" and "\" are changed to "_" to use probeset name as afile name. 
load(file = "results/adn/all_probes/E_GEOD_28146_all_probes.RData")

# Identify how namy probes have names containing "-" and "/"
E_GEOD_28146_no <- E_GEOD_28146_pr_init[!E_GEOD_28146_pr_init%in%E_GEOD_28146_pr]
length(E_GEOD_28146_no)#

# Extract original probesets' names from E-GEOD-29652
E_GEOD_29652<- nc_open("data/adn/E-GEOD-29652.nc");
data_E_GEOD_29652 <- ncvar_get(E_GEOD_29652, "data")
rownames(data_E_GEOD_29652) <- ncvar_get(E_GEOD_29652, "MetadataOrder")
colnames(data_E_GEOD_29652) <- ncvar_get(E_GEOD_29652,"gene")
data_E_GEOD_29652 <- t(data_E_GEOD_29652)
nc_close(E_GEOD_29652)
colnames(data_E_GEOD_29652)[1:18] <- "alz"

# Filter  out rows with SD values less then 0.29
SD <- apply(data_E_GEOD_29652, 1, sd, na.rm = T)
data_E_GEOD_29652_filt <- data_E_GEOD_29652[SD >= 0.29, ]
m <- t(data_E_GEOD_29652_filt)
ds_genes <- colnames(m)
E_GEOD_29652_pr_init <- ds_genes

# Load gsub corrected variant where "-" and "\" are changed to "_" to use probeset name as afile name.
load(file = "results/adn/all_probes/E_GEOD_29652_all_probes.RData")

# Identify how namy probes have names containing "-" and "/"
E_GEOD_29652_no <- E_GEOD_29652_pr_init[!E_GEOD_29652_pr_init%in%E_GEOD_29652_pr]
length(E_GEOD_29652_no)

# Extract original probesets' names from E_GEOD_4757
E_GEOD_4757<- nc_open("data/adn/E-GEOD-4757.nc");
data_E_GEOD_4757 <- ncvar_get(E_GEOD_4757, "data")
rownames(data_E_GEOD_4757) <- ncvar_get(E_GEOD_4757, "Phenotype")
colnames(data_E_GEOD_4757) <- ncvar_get(E_GEOD_4757,"gene")
data_E_GEOD_4757 <- t(data_E_GEOD_4757)
nc_close(E_GEOD_4757)
colnames(data_E_GEOD_4757) <- gsub("normal","norm", colnames(data_E_GEOD_4757))
colnames(data_E_GEOD_4757) <- gsub("neurofibriallary tangle","alz", colnames(data_E_GEOD_4757))
SD <- apply(data_E_GEOD_4757, 1, sd, na.rm = T)
data_E_GEOD_4757_filt <- data_E_GEOD_4757[SD >= 0.29, ]
m <- t(data_E_GEOD_4757_filt)
ds_genes <- colnames(m)
E_GEOD_4757_pr_init <- ds_genes

# Load gsub corrected variant where "-" and "\" are changed to "_" to use probeset name as afile name.
load(file = "results/adn/all_probes/E_GEOD_4757_all_probes.RData")

# Identify how namy probes have names containing "-" and "/"
E_GEOD_4757_no <- E_GEOD_4757_pr_init[!E_GEOD_4757_pr_init%in%E_GEOD_4757_pr]
length(E_GEOD_4757_no)

# Extract original probesets' names from E_GEOD_5281
E_GEOD_5281<- nc_open("data/adn/E-GEOD-5281.nc");
data_E_GEOD_5281 <- ncvar_get(E_GEOD_5281 , "data")
rownames(data_E_GEOD_5281 ) <- ncvar_get(E_GEOD_5281 , "DiseaseState")
colnames(data_E_GEOD_5281 ) <- ncvar_get(E_GEOD_5281 ,"gene")
data_E_GEOD_5281 <- t(data_E_GEOD_5281)
nc_close(E_GEOD_5281)
colnames(data_E_GEOD_5281) <- gsub("normal","norm", colnames(data_E_GEOD_5281))
colnames(data_E_GEOD_5281) <- gsub("Alzheimer's Disease","alz", colnames(data_E_GEOD_5281))

# Filter  out rows with SD values less then 0.29
SD <- apply(data_E_GEOD_5281, 1, sd, na.rm = T)
data_E_GEOD_5281_filt <- data_E_GEOD_5281[SD >= 0.29, ]
m <- t(data_E_GEOD_5281_filt)
ds_genes <- colnames(m)
E_GEOD_5281_pr_init <- ds_genes

# Load gsub corrected variant where "-" and "\" are changed to "_" to use probeset name as afile name.
load(file = "adn/all_probes/E_GEOD_5281_all_probes.RData")

# Identify how namy probes have names containing "-" and "/"
E_GEOD_5281_no <- E_GEOD_5281_pr_init[!E_GEOD_5281_pr_init%in%E_GEOD_5281_pr]
length(E_GEOD_5281_no)

# Extract original probesets' names from E_MEXP_2280
E_MEXP_2280<- nc_open("data/adn/E-MEXP-2280.nc");
data_E_MEXP_2280 <- ncvar_get(E_MEXP_2280 , "data")
rownames(data_E_MEXP_2280) <- ncvar_get(E_MEXP_2280 , "DiseaseState")
colnames(data_E_MEXP_2280) <- ncvar_get(E_MEXP_2280 ,"gene")
data_E_MEXP_2280 <- t(data_E_MEXP_2280)
nc_close(E_MEXP_2280)
data_E_MEXP_2280 <- data_E_MEXP_2280[,colnames(data_E_MEXP_2280)%in%c("Alzheimer's disease", "normal")]
colnames(data_E_MEXP_2280) <- gsub("normal","norm", colnames(data_E_MEXP_2280))
colnames(data_E_MEXP_2280) <- gsub("Alzheimer's Disease","alz", colnames(data_E_MEXP_2280))

# Filter  out rows with SD values less then 0.29
SD <- apply(data_E_MEXP_2280, 1, sd, na.rm = T)
data_E_MEXP_2280_filt <- data_E_MEXP_2280[SD >= 0.29, ]
m <- t(data_E_MEXP_2280_filt)
ds_genes <- colnames(m)
E_MEXP_2280_pr_init <- ds_genes

# Load gsub corrected variant where "-" and "\" are changed to "_" to use probeset name as afile name.
load (file = "results/adn/all_probes/E_MEXP_2280_all_probes.RData")
E_MEXP_2280_no <- E_MEXP_2280_pr_init[!E_MEXP_2280_pr_init%in%E_MEXP_2280_pr]
length(E_MEXP_2280_no)

# Find all affys that differ from original names
initial_probes <- unique(c(E_MEXP_2280_pr_init,E_GEOD_5281_pr_init,E_GEOD_4757_pr_init,E_GEOD_18309_pr_init,E_GEOD_28146_pr_init,E_GEOD_29652_pr_init))
#save(initial_probes, file = "results/adn/all_probes/initial_probes_6ds.RData")

# Create a dataframe with original names and the names with substitured "-" and "/" with "_"
gsub_probes <- gsub("-", "_",initial_probes)
gsub_probes <- gsub("-", "_",gsub_probes)

# Bind original name sand converted names
test <- cbind(initial_probes, gsub_probes)
test <- as.data.frame(test)

# Convert to character
test$initial_probes <- as.character(test$initial_probes)
test$gsub_probes <- as.character(test$gsub_probes)

# Find nonmatching rows
test_non_eq_rows <- test[!test$initial_probes==test$gsub_probes, ]

# Convert the missing IDs to ensg.
# Load the converted
load(file = "results/adn/all_probes/affy2ensg.RData")

# Identify the missing ones
# Convert affy ids to lower case
affy2ensg[,1] <- tolower(affy2ensg[,1])
head(affy2ensg)

missing_ids <- test$initial_probes[!test$initial_probes%in%affy2ensg$.id]

# Convert missing_ids  to ENSG
affy2ensg_missing <- gconvert(missing_ids)
affy2ensg_missing <- affy2ensg_missing[, c(2,4)]
colnames(affy2ensg_missing)<-c(".id", "Target")
affy2ensg_missing <- affy2ensg_missing[!duplicated(affy2ensg_missing), ]

# Dimentions
dim(affy2ensg_missing)
affy2ensg_missing[,1] <- tolower(affy2ensg_missing[,1])

# Bind with the previously converted data
affy2ensg_all <- rbind(affy2ensg, affy2ensg_missing)
affy2ensg_all[,2] <- as.character(affy2ensg_all[,2])

save(affy2ensg_all, file = "results/adn/all_probes/affy2ensg_all.RData")
write.table(affy2ensg_all, file = "results/adn/all_probes/affy2ensg_all.txt", sep = "\t", quote = F, row.names = F)

# Select the part of the original dataset, that has been converted to ensgs
test$initial_probes <- tolower(test$initial_probes)
test$gsub_probes <- tolower(test$gsub_probes) 
selected_affys_ensg <- merge(test,affy2ensg_all, by.x = "initial_probes", by.y = ".id", all = F)
dim(selected_affys_ensg) 
selected_affys_ensg <- selected_affys_ensg[!duplicated(selected_affys_ensg), ]
dim(selected_affys_ensg)

print("Writing results.")
# Save to file
save(selected_affys_ensg, file = "results/adn/all_probes/selected_affys_ensg.RData")
write.table(selected_affys_ensg, file = "results/adn/all_probes/selected_affys_ensg.txt", sep = "\t", quote = F, row.names = F)

# Unique selected affy probes for the analysis
un_sel_gsub_pr <- unique(selected_affys_ensg$gsub_probes)
length(un_sel_gsub_pr)
save(un_sel_gsub_pr, file = "results/adn/all_probes/un_sel_gsub_pr.RData")
