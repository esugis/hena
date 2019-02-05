print("ANALYSIS OF GENE CO-EXPRESSION IN ALZHEIMER'S DISEASE AND NORMAL SAMPLES")
print("BE AWARE THAT THE COMPUTATION OF PAIRWISE CO-EXPRESSION IS TIME AND MEMORY CONSUMING!")
print("DUE TO THE LARGE DATASET SIZE RUN THIS SCRIPT ON THE SERVER.")

print("Data preprocessing")
print("Preprocessing dataset 1 out of 6. Computing pairwise co-expression.")
source("scripts/adn/E_GEOD_4757.R")
# Clear workspace
rm(list = ls())

print("Preprocessing dataset 2 out of 6. Computing pairwise co-expression.")
source("scripts/adn/E_GEOD_5281.R")
# Clear workspace
rm(list = ls())

print("Preprocessing dataset 3 out of 6. Computing pairwise co-expression.")
source("scripts/adn/E_GEOD_18309.R")
# Clear workspace
rm(list = ls())

print("Preprocessing dataset 4 out of 6. Computing pairwise co-expression.")
source("scripts/adn/E_GEOD_28146.R")
# Clear workspace
rm(list = ls())

print("Preprocessing dataset 5 out of 6. Computing pairwise co-expression.")
source("scripts/adn/E_GEOD_29652.R")
# Clear workspace
rm(list = ls())

print("Preprocessing dataset 6 out of 6.Computing pairwise co-expression.")
source("scripts/adn/E_MEXP_2280.R")
# Clear workspace
rm(list = ls())

print("Ranking the co-expressed values for each probeset in each of the datasets. Aggregating the ranks in all the datasets.")
source("scripts/adn/RRA_probesets.R")
# Clear workspace
rm(list = ls())

print("Converting probe IDs to ENSG IDs")
source("scripts/adn/affy2ensg.R")
# Clear workspace
rm(list = ls())

print("Assembling together calculated RRA scores for all the probes.")
source("scripts/adn/coexp_int.R")
# Clear workspace
rm(list = ls())

print("Removing “self loops”(co-expression of gene with itself) from the co-expression dataset.")
source("scripts/adn/coexp2undirrected_selfloops_rm.R")
# Clear workspace
rm(list = ls())

print("Almost there. Final tuning. Adding columns interaction_type, data_source.")
source("scripts/adn/final_adn.R")
# Clear workspace
rm(list = ls())
