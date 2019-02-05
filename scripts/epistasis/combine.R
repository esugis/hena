# Script combines all epistatic data into one ds.

# Create the folder where current results will be written
resdir <- "results/epistasis/"
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

print("Loading individual results")
# Load ADNI VER data
load(file = "epi_adni_ver_int.RData")

# Load TGEN data
load(file = "epi_tgen_int.RData")

# Load HBTRC data
load(file = "epi_hbtrc_int.RData")

# Load ADNI cognitive data
load(file = "epi_adni_cog_int.RData")

# Combine individual epistatic data together
print("Combining the results")
epistasis_all_int <- rbind(epi_adni_ver_int, epi_tgen_int, epi_hbtrc_int, epi_adni_cog_int)
dim(epistasis_all_int)

print("Writing to file")
save(epistasis_all_int, file = "epistasis_all_int.RData")
write.table(epistasis_all_int, file = "epistasis_all_int.txt", sep = "\t", quote = F, row.names = F)
dim(epistasis_all_int)

setwd("../../")
