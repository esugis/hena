print("INTEGRATION OF THE INDIVIDUAL DATASETS OF BIOLOGICAL INTERACTIONS AND NODE ATTRIBUTES.")

print("Integrating the resulted individual datasets of biological interactions, i.e. co-expression interactions, epistatic interactions and protein-protein interactions.")
source("scripts/integration/integrate_interactions.R")

#Clean workspace
rm(list = ls())
print("Integrating biological attributes from the corresponding results of the analysis, i.e GWAS, positve selection, aggregated gene expression in brain regions, gene biotype and gene name. ")
source("scripts/integration/integrate_node_attributes.R")

print("Plotting dataset statisctics")
source("scripts/integration/plot_stats.R")

print("Prepare dataset for cytoscape and demonstration on the web)
source("scripts/integration/prepare_web_data_to_display.r")
