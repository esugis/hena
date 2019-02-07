# Script preprocesses file containing epistatic interactions in cognitive traits from ADNI cohort.
library(gProfileR)

# Read the file related to ADNI cohort cognitive traits
print("Loading the data")
epi <- read.table(file = "data/epistasis/ADNI_CT_epistasis.txt", header = T, stringsAsFactors = FALSE)

print("Performing analysis")
# Add ADNI general cohort code
epi$source <- paste("ADNI", epi$source, sep = "_")

# Number of cognitive traits
length(unique(epi$source))

# Create the folder where current results will be written
resdir <- "results/epistasis/"
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# See the structure of the data
str(epi)

# Select all the interactions that has intergenic regions (IGR)
epi_a <- epi[grep("-ENSG", epi$bin_1), ]
epi_b <- epi[grep("-ENSG", epi$bin_2), ]

# Number of IGR genes among interactors A
dim(epi_a)
dim(epi_b)

# Select the rows containing IGR
rows_cut_a <- row.names(epi_a)
rows_cut_b <- row.names(epi_b)
rows_cut <- unique(c(rows_cut_a, rows_cut_b))
length(rows_cut) 

# Select only bin_1, bin_2, source, bonf_pvalue from rows with IGR 
epi_adni_igri_tmp <- epi[rownames(epi)%in%rows_cut, c(1, 2, 9, 10)]

# Bind columns with interaction type, direction, evidence code.
# Interaction type named intergenic region interaction (IGRI)
epi_adni_igri <- cbind(epi_adni_igri_tmp[,c(1,2,4)], interaction_type = "IGRI")

# Data source
epi_adni_igri <- cbind(epi_adni_igri, data_source = epi_adni_igri_tmp$source) 

# df that has intercators as intergenic regions in one of the columns and correcponding scores
colnames(epi_adni_igri) <- c("ensg1","ensg2","score","interaction_type","data_source")

# Save with the indentifier related to cognitive traits
epi_adni_cog_igri <- epi_adni_igri

dim(epi_adni_cog_igri)

# Combine with the main data frame 

# Select only part related to  bin_1, bin_2, source, bonf_pvalue from the main df
epi_adni_cog_tmp <- epi[!row.names(epi)%in%rows_cut, c(1, 2, 9, 10)]
dim(epi_adni_cog_tmp) 

# Bind columns with interaction type, direction, evidence code.
# Interaction type intergenic region interaction
epi_adni_cog <- cbind(epi_adni_cog_tmp[,c(1,2,4)], interaction_type = "epistasis")
# Data source
epi_adni_cog <- cbind(epi_adni_cog, data_source = epi_adni_cog_tmp$source)
# Combine
colnames(epi_adni_cog) <- c("ensg1","ensg2","score","interaction_type","data_source")

# Convert gene ids and ensg id to tha latest Ensembl version
# Extract current ensg ids 
curr_id <- unique(c(epi_adni_cog$ensg1,epi_adni_cog$ensg2 ))
length(curr_id)

# Convert ids using biomart
library(biomaRt)
mart.pr <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "ensembl.org")
curr_id2ensg <- getBM(attributes = c("ensembl_gene_id"), filters=c("ensembl_gene_id"), values = curr_id, mart = mart.pr)

curr_id2ensg <-as.character(as.vector(curr_id2ensg$ensembl_gene_id))

epi_adni_cog_2ensg <- epi_adni_cog[epi_adni_cog$ensg1%in%curr_id2ensg,]
epi_adni_cog_2ensg <- epi_adni_cog_2ensg[epi_adni_cog_2ensg$ensg2%in%curr_id2ensg,]

# Size of the dataset with old ENSG IDs(ver 74) converted to the latest ENSG version
dim(epi_adni_cog_2ensg)
# Size of the dataset with old ENSG IDs(ver74)
dim(epi_adni_cog)

# Find differences between ENSG IDs in the latest Ensembl version and Ensembl ver 74
# All ENSG IDs in the latest version
ensg_adni_cog_ensg<-unique(c(as.character(epi_adni_cog_2ensg$ensg1),as.character(epi_adni_cog_2ensg$ensg2)))
length(ensg_adni_cog_ensg)

# All ENSG IDs in ver 74
ensg_adni_cog_ver74<-unique(c(as.character(epi_adni_cog$ensg1),as.character(epi_adni_cog$ensg2)))
length(ensg_adni_cog_ver74)

# ENSG IDs present in ver 74 and not present in the latest version
ensg_lv_vs_74_adni_cog<-ensg_adni_cog_ver74[!ensg_adni_cog_ver74%in%ensg_adni_cog_ensg]
length(ensg_lv_vs_74_adni_cog)

# Find differences in IGRI
epi_adni_cog_igri_ensg<-unique(c(as.character(epi_adni_cog_igri$ensg1),as.character(epi_adni_cog_igri$ensg2)))
length(epi_adni_cog_igri_ensg)#
library(stringr)
curr_igri_ensg<-unique(unlist(str_split(epi_adni_cog_igri_ensg,"-")))
# Number of the current version ENSG ids
length(curr_igri_ensg)

#Convert gene ids and ensg id to the latest Ensembl using biomart
mart.pr <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "www.ensembl.org")
epi_adni_cog_igri_ensg2ensg <- getBM(attributes = c("ensembl_gene_id"), filters=c("ensembl_gene_id"), values = curr_igri_ensg, mart = mart.pr)
epi_adni_cog_igri_ensg_lv <-as.character(as.vector(epi_adni_cog_igri_ensg2ensg$ensembl_gene_id))

# Old ensg IDs
epi_adni_cog_igri_ensg_ver74 <- curr_igri_ensg

# Find the differebce between old and new
ensg_lv_vs_74_adni_cog_igri<-epi_adni_cog_igri_ensg_ver74[!epi_adni_cog_igri_ensg_ver74%in%epi_adni_cog_igri_ensg_lv]
length(ensg_lv_vs_74_adni_cog_igri)

# Find overlap between missing IDs in IGRI and new gene IDs
ensg_lv_vs_74_adni_cog_all<-unique(c(ensg_lv_vs_74_adni_cog, ensg_lv_vs_74_adni_cog_igri))
length(ensg_lv_vs_74_adni_cog_all)

# Exclude the rows with old ids that are not mapped to the latest Ensembl version 
excluderow<-c()
for(i in 1:length(rownames(epi_adni_cog_igri))){
 pair1 <- unique(unlist(str_split(as.character(epi_adni_cog_igri[i,1]),"-")))
 p1<-c(pair1%in%ensg_lv_vs_74_adni_cog_all)
 c1<-length(p1[p1==T])
 pair2 <- unique(unlist(str_split(as.character(epi_adni_cog_igri[i,2]),"-")))
 p2<-c(pair2%in%ensg_lv_vs_74_adni_cog_all)
 c2<-length(p2[p2==T])
 if(c1>0){
  excluderow<-c(excluderow,i)
  }else if(c2>0){
        excluderow<-c(excluderow,i)
        }
}

epi_adni_cog_igri_lv <- epi_adni_cog_igri[!rownames(epi_adni_cog_igri)%in%excluderow,]

dim(epi_adni_cog_igri_lv)

# Combine dataframes of IGRI and the rest genes interactions. Write final ds to the file
epi_adni_cog_lv <- epi_adni_cog_2ensg
epi_adni_cog_int <- rbind(epi_adni_cog_lv, epi_adni_cog_igri_lv)
dim(epi_adni_cog_int)


# Remove the duplicated undirrescted edges with the same score.
# For example ENSG1-ENSG2 0.5 and ENSG2-ENSG 0.5
# Convert factors to characters
df2string<-function(df){
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
df[,3]<-as.numeric(df[,3])
return (df)}

# Epistatic interactions
epi_adni_cog_int <- df2string(epi_adni_cog_int)
str(epi_adni_cog_int)
dim(epi_adni_cog_int)

# Remove the duplicated undirrescted edges with the same score.
# For example ENSG1-ENSG2 0.5 and ENSG2-ENSG 0.5
epi_adni_cog_int <- epi_adni_cog_int[!duplicated(data.frame(t(apply(epi_adni_cog_int[1:2], 1, sort)), epi_adni_cog_int[,c(3,5)])),]
# New size
dim(epi_adni_cog_int)

# Clean dataset
rownames(epi_adni_cog_int) <- 1:nrow(epi_adni_cog_int)
rm_ensg <- c(epi_adni_cog_int$ensg1,epi_adni_cog_int$ensg2)
epi_adni_cog_int <- epi_adni_cog_int[!epi_adni_cog_int$ensg1%in%rm_ensg[grep("^-.*",rm_ensg)],]
epi_adni_cog_int <- epi_adni_cog_int[!epi_adni_cog_int$ensg2%in%rm_ensg[grep("^-.*",rm_ensg)],]

epi_adni_cog_int <- epi_adni_cog_int[!epi_adni_cog_int$ensg1%in%rm_ensg[grep("\\ENSG.*-$",rm_ensg)],]
epi_adni_cog_int <- epi_adni_cog_int[!epi_adni_cog_int$ensg2%in%rm_ensg[grep("\\ENSG.*-$",rm_ensg)],]

# Save the part of the integrated dataset related to cognitive traits studies in ADNI cohort
print("Writing results to the file")
save(epi_adni_cog_int, file = "epi_adni_cog_int.RData")
write.table(epi_adni_cog_int,file = "epi_adni_cog_int.txt",sep = "\t", quote = F, row.names = F)

setwd("../../")



