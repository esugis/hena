print("ANALYSIS of EPISTASIS")

print("Analysing dataset 1 out of 4")
source("scripts/epistasis/epi_hbtrc.R")

print("Analysing dataset 2 out of 4")
source("scripts/epistasis/epi_adni.R")

print("Analyssing dataset 3 out of 4")
source("scripts/epistasis/epi_tgen.R")

print ("Analysing dataset 4 out of 4")
source("scripts/epistasis/epi_adni_cog.R")

print ("Combining analysed data into one dataset")
source("scripts/epistasis/combine.R")

