# This script:
# Extracts Alzheimer’s related and healthy samples from the data sets.
# Filters out probesets with SD < 0.29.
# Calculating the co-expression between all probesets in each of the data sets using Spearman correlation coefficient.
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
E_GEOD_4757 <- nc_open("data/adn/E-GEOD-4757.nc");

print("Performing analysis.")
# Extract only Alzheimer and healthy saples from E_GEOD_4757
# List of variable in E_GEOD_4757
#names(E_GEOD_4757$var)
#ncvar_get(E_GEOD_4757, "Phenotype")

# Select only Alzheimer's disease and healthy samples
data_E_GEOD_4757 <- ncvar_get(E_GEOD_4757, "data")
rownames(data_E_GEOD_4757) <- ncvar_get(E_GEOD_4757, "Phenotype")
colnames(data_E_GEOD_4757) <- ncvar_get(E_GEOD_4757,"gene")
data_E_GEOD_4757 <- t(data_E_GEOD_4757)

# Close NetCDF file
nc_close(E_GEOD_4757)

# Rename the samples
colnames(data_E_GEOD_4757) <- gsub("normal","norm", colnames(data_E_GEOD_4757))
colnames(data_E_GEOD_4757) <- gsub("neurofibriallary tangle","alz", colnames(data_E_GEOD_4757))

# Filter out rows with SD values less then 0.29
SD <- apply(data_E_GEOD_4757, 1, sd, na.rm = T)
data_E_GEOD_4757_filt <- data_E_GEOD_4757[SD >= 0.29, ]

m <- t(data_E_GEOD_4757_filt)

# Path to the foldet where results will be saved
pathRdata <- "results/adn/all_probes/rdata/E_GEOD_4757/"
pathtxt <- "results/adn/all_probes/txt/E_GEOD_4757/"

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
  cor1gds <- c()#corelation for one gene in one ds
  gene <- ds_genes[i];
  vect <- m[,i];
   cor1gds <- t(cor(vect,m,method="spearman"))
   # Save the result as txt and RData
   filename <- sprintf("%s.txt",gene);
   filedata <- sprintf("%s.RData",gene);
   pathname <- file.path(pathtxt, filename);
   pathdata <- file.path(pathRdata, filedata);
   write.table(cor1gds, file = pathname, sep = "\t", quote = F, row.names = F);
   save(cor1gds,file = pathdata)
 }

# Exract the list of probesets
E_GEOD_4757_pr <- colnames(m)
print("Writing results to file")
save(E_GEOD_4757_pr, file = "results/adn/all_probes/E_GEOD_4757_all_probes.RData")


