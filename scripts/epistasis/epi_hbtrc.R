## Script preprocesses file containing epistatic interactions in HBTRC cohort.
library(gProfileR)

# Read the file related to HTBRC cohort
print("Loading the data ")
epi <- read.table(file="data/epistasis/HBTRC_epistasis.tsv", header=T)

# Create the folder where current results will be written
resdir <- "results/epistasis/"
dir.create(file.path(resdir),showWarnings  =  FALSE, recursive  =  TRUE)

# Set created directory as working dirrectory
setwd(resdir)

print("Performning analysis")
epi_a<- epi[grep("LRG_", epi$ENSG_A), ]
epi_b<- epi[grep("LRG_", epi$ENSG_B), ]
rows_cut=row.names(epi_b)

# Convert LRG to ENSG
B=as.character(epi_b$ENSG_B)
LRG2ensg_hbtrc = gconvert(B)
LRG2ensg_hbtrc <- LRG2ensg_hbtrc[,c(2,4)]
colnames(LRG2ensg_hbtrc) <-c(".id", "Target")

# Remove duplicates
LRG2ensg_hbtrc=LRG2ensg_hbtrc[!duplicated(LRG2ensg_hbtrc), ]

# Megre
epi_hbtrc_lrg=merge(epi_b,LRG2ensg_hbtrc, by.x="ENSG_B", by.y=".id", all=F)

# Select only ensg1 ensg2 and score
epi_hbtrc_lrg_fin=epi_hbtrc_lrg[,c(2,8,3)]

# Data head
epi_hbtrc_lrg_fin 

# Bind columns with interaction_type, data_source
epi_hbtrc_lrg_fin=cbind(epi_hbtrc_lrg_fin, interaction_type="epistasis")
epi_hbtrc_lrg_fin=cbind(epi_hbtrc_lrg_fin, data_source="HBTRC")

# Rename the coluns
colnames(epi_hbtrc_lrg_fin)=c("ensg1","ensg2","score","interaction_type","data_source")

# Combine with the main data frame 
epi_cut=epi[!row.names(epi)%in%rows_cut,]
epi_cut_fin=epi_cut[,c(1,2,3)]
epi_cut_fin=cbind(epi_cut_fin,interaction_type="epistasis")
epi_cut_fin=cbind(epi_cut_fin, data_source="HBTRC")
colnames(epi_cut_fin)=c("ensg1","ensg2","score","interaction_type","data_source")
epi_hbtrc=epi_cut_fin

# Convert gene ids and ensg id to tha latest Ensembl version
epi_hbtrc_ensg12ensg <- gconvert(epi_hbtrc$ensg1)
epi_hbtrc_ensg12ensg <- epi_hbtrc_ensg12ensg[,c(2,4)]
colnames(epi_hbtrc_ensg12ensg)<-c(".id","Target")

# Remove duplicates
epi_hbtrc_ensg12ensg <- epi_hbtrc_ensg12ensg[!duplicated(epi_hbtrc_ensg12ensg), ]

# Convert ids
epi_hbtrc_ensg22ensg <- gconvert(epi_hbtrc$ensg2)
epi_hbtrc_ensg22ensg <- epi_hbtrc_ensg22ensg[,c(2,4)]
colnames(epi_hbtrc_ensg22ensg) <- c(".id", "Target")

# Remove duplicates
epi_hbtrc_ensg22ensg <- epi_hbtrc_ensg22ensg[!duplicated(epi_hbtrc_ensg22ensg), ]

# Merge by ensg1
epi_hbtrc_2ensg <- merge(epi_hbtrc,epi_hbtrc_ensg12ensg, by.x="ensg1", by.y=".id", all=F)

# Merge by ensg2
epi_hbtrc_2ensg <- merge(epi_hbtrc_2ensg,epi_hbtrc_ensg22ensg, by.x="ensg2", by.y=".id", all=F)

# Find differences between ENSG IDs in Ensembl ver 90 and Ensembl ver 74
# All ENSG IDs in ver 90
ensg_hbtrc_ver90 <- unique(c(as.character(epi_hbtrc_2ensg$ensg1),as.character(epi_hbtrc_2ensg$ensg2)))


# All ENSG IDs in ver 74
ensg_hbtrc_ver74 <- unique(c(as.character(epi_hbtrc$ensg1),as.character(epi_hbtrc$ensg2)))

# ENSG IDs present in ver 74 and not present in ver 90
ensg_90_vs_74_hbtrc <- ensg_hbtrc_ver74[!ensg_hbtrc_ver74%in%ensg_hbtrc_ver90]

epi_hbtrc <- epi_hbtrc_2ensg[,c(6,7,3,4,5)]
colnames(epi_hbtrc) <- c("ensg1","ensg2","score","interaction_type","data_source")
epi_hbtrc <- epi_hbtrc[!duplicated(epi_hbtrc), ]

# Combine dataframes of LRG genes interactions and the rest
epi_hbtrc_int <- rbind(epi_hbtrc,epi_hbtrc_lrg_fin)

# Remove the duplicated undirrescted edges with the same score.
# For example ENSG1-ENSG2 0.5 and ENSG2-ENSG1 0.5

# Convert factors to characters
df2string<-function(df){
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
df[,3]<-as.numeric(df[,3])
return (df)}

epi_hbtrc_int <- df2string(epi_hbtrc_int)
epi_hbtrc_int <- epi_hbtrc_int[!duplicated(data.frame(t(apply(epi_hbtrc_int[1:2], 1, sort)), epi_hbtrc_int[,c(3,5)])),]

# Save the part of the integrated dataset related to hbtrc cohort
print("Writing to file")
save(epi_hbtrc_int, file="epi_hbtrc_int.RData")
write.table(epi_hbtrc_int,file="epi_hbtrc_int.txt",sep="\t", quote=F, row.names=F)

setwd("../../")

