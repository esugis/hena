# This script puts together calculated coexpr scores for all affy probes.

# Create the folder where current results will be written
resdir <- "results/adn/integration/"
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Used libarries
library(foreach); library(doMC); cores=10 ; registerDoMC(cores);
library(stringr);
library("R.utils")
library(gProfileR)

# Set the path to the calculated RRA scores
pathRdata <- "results/adn/all_probes/scores/rdata"

# Load unique selected affy probes
print("Loading unique selected affy probes.")
load(file = "results/adn/all_probes/un_sel_gsub_pr.RData")
load(file = "results/adn/all_probes/selected_affys_ensg.RData")

print("Performing analysis")
# Remove unwanted symbols
un_sel_gsub_pr_fnames <- gsub("/","_",un_sel_gsub_pr)
un_sel_gsub_pr_fnames <- gsub("(affx*)(\\w*)(\\_at)", "\\U\\1\\U\\2\\L\\3", un_sel_gsub_pr_fnames, perl = TRUE)

# Interactions between affys, converted to ensg and stored in one df
# Length of the selected unique probesets
ll <- length(un_sel_gsub_pr)
for(i in 1:ll){

# Get the name
gene<-un_sel_gsub_pr[i]
genefname <- un_sel_gsub_pr_fnames[i]
print(gene)

# Set the path to the file
filedata <- sprintf("%s.RData",genefname);
pathdata <- file.path(pathRdata, filedata);

# Load the file
load(file = pathdata)
colnames(ar_gene)[1] <-"affy2"

# Select the rows that are in the un_sel_gsub_p
ar_gene <- ar_gene[ar_gene$affy2%in%un_sel_gsub_pr,]
ar_gene <- cbind(affy1 = sprintf("%s",gene),ar_gene)

# Convert to lower case
ar_gene$affy2 <- tolower(ar_gene$affy2)
ar_gene$affy1 <- tolower(ar_gene$affy1)

# Merge with corresponding ensg ids
ar_gene_merges <- merge(ar_gene,selected_affys_ensg, by.x = "affy1", by.y = "gsub_probes",all = F)
ar_gene_merged <- merge(ar_gene_merges,selected_affys_ensg, by.x = "affy2", by.y = "gsub_probes",all = F)
ar_gene_merged <- ar_gene_merged[,c(2,1,6,8,3,4)]

# Cut off p-value > 0.00001 instead of 0.05
ar_gene_merged <- ar_gene_merged[ar_gene_merged$adj.pval<=0.00001,]

coexp <- ar_gene_merged[,3:6]
colnames(coexp)[1:2] <- c("ensg1","ensg2")
coexp <- coexp[!duplicated(coexp), ]
coexp$gene <- paste(coexp$ensg1, coexp$ensg2, sep = "_")
coexp_num <- coexp[,c(4,5)]

if(nrow(coexp_num)>0){

# Aggregate the scores in case few probesets match to the same gene
coexp_num <- aggregate(. ~ gene, coexp_num, max)
coexp_num$ensg1 <- unlist(lapply(strsplit(as.character(coexp_num$gene), "_"), "[", 1))
coexp_num$ensg2 <- unlist(lapply(strsplit(as.character(coexp_num$gene), "_"), "[", 2))
coexp_num <- coexp_num[,c(3,4,2)]
coexp_num <- cbind(coexp_num, interaction_type = "coexpression")

# Alzheimer's disease and normal samples
coexp_num <- cbind(coexp_num, data_source = "ADN")
colnames(coexp_num) <- c("ensg1","ensg2", "score", "interaction_type","data_source")
coexp_full_ann <- cbind(coexp_num,datasets = "E_GEOD_18309,E_GEOD_28146,E_GEOD_29652,E_GEOD_4757,E_GEOD_5281,E_MEXP_2280")
coexp_short <- coexp_num

# Add colnames
colnames(coexp_short) <- c("ensg1","ensg2","score", "interaction_type","data_source")
coexp_short$interaction_type <- as.character(as.vector(coexp_short$interaction_type))
coexp_short$data_source <- as.character(as.vector(coexp_short$data_source))

# Coexp table for the integration
write.table(coexp_short,file = "results/adn/integration/coexp_signif_int.txt",append = TRUE,sep = "\t", quote = F, row.names = F,col.names = F)
}}

# Check
alzcoexp_int <- read.table(file = "results/adn/integration/coexp_signif_int.txt",sep = "\t",header = T)

# Save in RData format
print("Writing to file.")
save(alzcoexp_int, file = "results/adn/integration/alzcoexp_int.RData")

