# This script:
# Extracts Alzheimer’s related and healthy samples from the data sets.
# Filters out probesets with SD < 0.29.
# Calculating the co-expression between  all probesets in each of the data sets using Spearman correlation coefficient.
# For each individual probeset in each of the datasets script creates 2 separate files in .txt and .RData formats.
# Files are named after the probeset.
# Created files contain the names of the correlated probesets and the corresponding Spearman coefficient.


# Used libarries
library(foreach); library(doMC); cores=10 ; registerDoMC(cores);
library(stringr);
library("R.utils");
library(ncdf4);

# Open dataset in NetCDF format
# Please indicate the path to the saved .nc file, e.g. as shown below
print("Loading data.")
E_MEXP_2280 <- nc_open("data/adn/E-MEXP-2280.nc");

print("Performing analysis.")
# Extract only Alzheimer and healthy saples from E_GEOD_5281
# List of variable in E_GEOD_5281 
#names(E_MEXP_2280$var)
#ncvar_get(E_MEXP_2280 , "DiseaseState")

# Select only Alzheimer's disease and healthy  samples
data_E_MEXP_2280 <- ncvar_get(E_MEXP_2280 , "data")
rownames(data_E_MEXP_2280) <- ncvar_get(E_MEXP_2280 , "DiseaseState")
colnames(data_E_MEXP_2280) <- ncvar_get(E_MEXP_2280 ,"gene")
data_E_MEXP_2280 <- t(data_E_MEXP_2280)

# Close NetCDF file
nc_close(E_MEXP_2280)

# Rename the samples
data_E_MEXP_2280 <- data_E_MEXP_2280[,colnames(data_E_MEXP_2280)%in%c("Alzheimer's disease", "normal")]
colnames(data_E_MEXP_2280) <- gsub("normal","norm", colnames(data_E_MEXP_2280))
colnames(data_E_MEXP_2280) <- gsub("Alzheimer's Disease","alz", colnames(data_E_MEXP_2280))

# Filter  out rows with SD values less then 0.29
SD <- apply(data_E_MEXP_2280, 1, sd, na.rm = T)
data_E_MEXP_2280_filt <- data_E_MEXP_2280[SD >= 0.29, ]

m<-t(data_E_MEXP_2280_filt)
length(ds_genes <- colnames(m))

# Path to the foldet where results will be saved
pathRdata <- "results/adn/all_probes/rdata/E_MEXP_2280/"
pathtxt <- "results/adn/all_probes/txt/E_MEXP_2280/"

# Create directories
dir.create(file.path(pathRdata),showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(pathtxt),showWarnings = FALSE, recursive = TRUE)

# Substitute "-" and "/" in probesets' names with "_"
colnames(m) <- gsub("-", "_",colnames(m))
colnames(m) <- gsub("/", "_",colnames(m))

# Probesets  
ds_genes <- colnames(m)

# Compute Sperman correlation between expression profiles of probesets "all against all".

foreach(i = 1:length(ds_genes)) %dopar%{
  cor1gds <- c() # Corelation for one gene in one ds
   gene<-ds_genes[i];
   vect <- m[,i];
   cor1gds <- t(cor(vect, m, method = "spearman"))
   # Save it as RData
   filename <- sprintf("%s.txt", gene);
   filedata <- sprintf("%s.RData", gene);
   pathname <- file.path(pathtxt, filename);
   pathdata <- file.path(pathRdata, filedata);
   write.table(cor1gds, file = pathname, sep = "\t", quote=F, row.names = F);
   save(cor1gds,file = pathdata)
 }

# Exract the list of probesets
E_MEXP_2280_pr <- colnames(m)
print("Writing to file.")
save(E_MEXP_2280_pr, file = "results/adn/all_probes/E_MEXP_2280_all_probes.RData")

