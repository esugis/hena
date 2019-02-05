# This script combines all indiviadual interaction datasets into one and removes duplicated interaction of one type

# Create the folder where current results will be written
resdir <- paste("results","integration",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Load individual interaction datasets
print("Loading individual datasets of biological interactions.")
# PPIs associated with brain ageing (PBA)
load(file="results/pba/pba_int.RData")

# Combined epistatic interactions ADNI ventrical volume, ADNI cognitive traits, TGEN, HBTRC
load(file="results/epistasis/epistasis_all_int.RData")

# Intact hyman PPIs with MIscore >=0.45
load(file="results/intact/intact_int.RData")

## IntAct Alzheimer's related manually curated PPI dataset with MIscore >=0.45
load(file="results/intact/alz_intact_int.RData")

# IntAct Synapse related automatically curated PPI dataset with MIscore >=0.45
load(file="results/intact/syn_intact_int.RData")

# Co-expression dataset with removed self-loops RRA score <= 0.00001 
load(file="results/adn/integration/adn_coexp_int.RData")

# Co-expression in DG
load(file="results/allenbrain/DG_coexp_int_aba.RData")

# Co-expression	in CA1
load(file="results/allenbrain/CA1_coexp_int_aba.RData")

# Co-expression in CA2
load(file="results/allenbrain/CA2_coexp_int_aba.RData")

# Co-expression in CA3
load(file="results/allenbrain/CA3_coexp_int_aba.RData")

# Co-expression in CA4
load(file="results/allenbrain/CA4_coexp_int_aba.RData")

# Co-expression in subiculum (S)
load(file="results/allenbrain/subiculum_coexp_int_aba.RData")

# Co-expression in SptN
load(file="results/allenbrain/SptN_coexp_int_aba.RData")

# Create one DF from separate datasets
print("Combining individual datasets.")
integrated_int <- rbind(pba_int, epistasis_all_int, intact_int, alz_intact_int, syn_intact_int, adn_coexp_int, DG_coexp_int, CA1_coexp_int, CA2_coexp_int, CA3_coexp_int, CA4_coexp_int, subiculum_coexp_int, SptN_coexp_int)

rm(pba_int, epistasis_all_int, intact_int, alz_intact_int, syn_intact_int, adn_coexp_int, DG_coexp_int, CA1_coexp_int, CA2_coexp_int, CA3_coexp_int, CA4_coexp_int, subiculum_coexp_int, SptN_coexp_int)

integrated_int <- integrated_int[!duplicated(integrated_int),]
colnames(integrated_int)[1:2]<-c("ensg.A","ensg.B")

# Save sorted integrated dataset to the files.
print("Writing integrated dataset to file.")
save(integrated_int, file = "results/integration/integrated_int.RData")
write.table(integrated_int, file = "results/integration/integrated_int.txt", sep="\t", quote=F, row.names=F)
