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
E_GEOD_18309<- nc_open("data/adn/E-GEOD-18309.nc");

print("Performing analysis.")
# List of variables in the dataset
#names(E_GEOD_18309$var)
#ncvar_get(E_GEOD_18309, "MetadataOrder")

# Take only Alzheimer related samples based on the additional info in the ArrayExpress
data_E_GEOD_18309 <- ncvar_get(E_GEOD_18309, "data")
rownames(data_E_GEOD_18309) <- ncvar_get(E_GEOD_18309, "MetadataOrder")
colnames(data_E_GEOD_18309) <- ncvar_get(E_GEOD_18309,"gene")
data_E_GEOD_18309 <- t(data_E_GEOD_18309)

# Close NetCDF file
nc_close(E_GEOD_18309)

# Extract and rename only Alzheimer and healthy saples   
colnames(data_E_GEOD_18309)[1:3] <- "alz"
colnames(data_E_GEOD_18309)[7:9] <- "norm"
data_E_GEOD_18309 <- data_E_GEOD_18309[,c(1:3,7:9)]

# Filter out rows with SD values less then 0.29
SD <- apply(data_E_GEOD_18309, 1, sd, na.rm = T)
data_E_GEOD_18309_filt <- data_E_GEOD_18309[SD >= 0.29, ]

# Gene names
m <- t(data_E_GEOD_18309_filt)

# Path to the foldet where results will be saved
pathRdata <- "results/adn/all_probes/rdata/E_GEOD_18309/"
pathtxt <- "results/adn/all_probes/txt/E_GEOD_18309/"

# Create directories
dir.create(file.path(pathRdata),showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(pathtxt),showWarnings = FALSE, recursive = TRUE)

# Substitute "-" and "/" in the probesets names with "_"
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
E_GEOD_18309_pr <- colnames(m)
print("Write to file.")
save(E_GEOD_18309_pr, file = "results/adn/all_probes/E_GEOD_18309_all_probes.RData")


