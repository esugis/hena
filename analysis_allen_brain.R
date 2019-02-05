print("ANALYSIS OF GENE EXPRESSION IN WHOLE-BRAIN DATASETS FROM ALLEN BRAIN ATLAS")
print("BE AWARE THAT THESE COMPUTATIONS ARE MEMORY AN TIME CONSUMING. DUE TO THE LARGE DATASET SIZE AND LARGE NUMBER OF COMPUTATIONS PLEASE RUN THIS ANALYSIS ON THE SERVER.")
print("Converting probe ids to ENSG IDs")
source("scripts/allenbrain/p2ensg.R")

print("Computing z-scores based on gene expression in brain regions over all six brains")
source("scripts/allenbrain/expression_all_regions_zscores.R")

# Clear workspace
rm(list = ls())

print("Computing co-expression in Alzheimer's disease related brain regions")
print("Preprocessing")
source("scripts/allenbrain/coexp_per_tissue_preprocessing.R")

# Clear workspace
rm(list = ls())

print("Computing co-expression in Alzheimer's disease related brain regions.")
print("Be aware that it takes lots of computational time!")
source("scripts/allenbrain/coexp_in_selected_brain_regions.R")

# Clear workspace
rm(list = ls())
