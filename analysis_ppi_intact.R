print("ANALYSIS OF PROTEIN-PROTEIN INTERACTIONS FROM INTACT DATASETS")

print("Analysing dataset 1 out of 3. Analysis of PPI interactions in human.")
source("scripts/intact/intact.R")

print("Analysing dataset 2 out of 3. Analysis of expert curated PPIs related to Alzheimer's disease deposited in IntAct database.")
source("scripts/intact/alz_intact.R")

print("Analysing dataset 3 out of 3. Analysis of automatically curated synaptic PPIs deposited in IntAct.")
source("scripts/intact/synapse_intact.R")
