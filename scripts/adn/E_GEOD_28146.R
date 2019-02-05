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
E_GEOD_28146 <-nc_open("data/adn/E-GEOD-28146.nc");

print("Performing analysis.")
# Extract only Alzheimer and healthy saples from E_GEOD_28146
# List of variable in E_GEOD_28146
#names(E_GEOD_28146$var)
# List metadata related to disease status
#ncvar_get(E_GEOD_28146, "disease_status")

# Select only Alzheimer's disease and healthy samples
data_E_GEOD_28146 <- ncvar_get(E_GEOD_28146, "data")
rownames(data_E_GEOD_28146) <- ncvar_get(E_GEOD_28146, "MetadataOrder")
colnames(data_E_GEOD_28146) <- ncvar_get(E_GEOD_28146, "gene")
data_E_GEOD_28146 <- t(data_E_GEOD_28146)

# Close NetCDF file
nc_close(E_GEOD_28146)

# Rename the samples
colnames(data_E_GEOD_28146)[1:22] <- "alz"
colnames(data_E_GEOD_28146)[23:30] <- "norm"

# Filter out rows with SD values less then 0.29
SD <- apply(data_E_GEOD_28146, 1, sd, na.rm = T)
data_E_GEOD_28146_filt <- data_E_GEOD_28146[SD >= 0.29, ]

m <- t(data_E_GEOD_28146_filt)
length(ds_genes <- colnames(m))

# Path to the foldet where results will be saved
pathRdata <- "results/adn/all_probes/rdata/E_GEOD_28146/"
pathtxt <- "results/adn/all_probes/txt/E_GEOD_28146/"

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
   gene=ds_genes[i];
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
E_GEOD_28146_pr <- colnames(m)
print("Writing to file.")
save(E_GEOD_28146_pr, file = "results/adn/all_probes/E_GEOD_28146_all_probes.RData")


