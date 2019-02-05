# Tris script extracts separate interac

library(dplyr)
# Load the data
print("Loading data")
load("../results/integration/integrated_int.RData")
intg_int_all<-integrated_int

# Extract positive coexpression
print("Extracting positive co-expression")
integrated_int_coexp_all <- intg_int_all[intg_int_all$interaction_type%in%"coexpression",]
integrated_int_coexp_pos <- integrated_int_coexp_all[integrated_int_coexp_all$score>0,c(1,2)]
resdir <- "results/integrated_int"
dir.create(file.path(resdir), showWarnings = FALSE, recursive = TRUE)

write.table(integrated_int_coexp_pos, file= "results/integrated_int/integrated_int_coexp_pos.txt", quote=F, row.names=F, col.names=F, sep="\t")
#rm(integrated_int_coexp)
#rm(integrated_int_coexp_all)

print("Extracting co-expression")
integrated_int_coexp <- intg_int_all[intg_int_all$interaction_type%in%"coexpression",c(1,2)]

print("Saving co-expression interactions as a separate file")
resdir <- "results/integrated_int"
dir.create(file.path(resdir), showWarnings = FALSE, recursive = TRUE)

write.table(integrated_int_coexp, file= "results/integrated_int/integrated_int_coexp.txt", quote=F, row.names=F, col.names=F, sep="\t")
rm(integrated_int_coexp)


# Extract PPIs
print("Extracting PPIs")
integrated_int_ppi <- intg_int_all[intg_int_all$interaction_type%in%"PPI",c(1,2)]

print("Saving ppis as a separate file")
write.table(integrated_int_ppi, file= "results/integrated_int/integrated_int_ppi.txt", quote=F, row.names=F, col.names=F, sep="\t")
rm(integrated_int_ppi)

# Extract epistasis
print("Extracting epistatic interactions, excluding intergenic region interactions")
integrated_int_epi <- intg_int_all[intg_int_all$interaction_type%in%"epistasis",c(1,2)]

print("Saving epistatic interactions as a separate file")
write.table(integrated_int_epi, file= "results/integrated_int/integrated_int_epi.txt", quote=F, row.names=F, col.names=F, sep="\t")
rm(integrated_int_epi)

