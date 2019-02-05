# This script:
# Reads in all the coexpression matrixes calculated for each dataset
# Gets the list of genes from each of the filtered datasets
# For each gene loads correlation vector from each datasets
# Creates the ranks in each vector.
# Aggregates the ranks in vectors

# Path to the folder where selults will be stored
pathRdata <- "results/adn/all_probes/scores/rdata/"
dir.create(file.path(pathRdata),showWarnings = FALSE, recursive = TRUE)

# Load the probesets names
print("Loading data.")
load(file <- "results/adn/all_probes/E_GEOD_18309_all_probes.RData")
load(file <- "results/adn/all_probes/E_GEOD_28146_all_probes.RData")
load(file <- "results/adn/all_probes/E_GEOD_29652_all_probes.RData")
load(file <- "results/adn/all_probes/E_GEOD_4757_all_probes.RData")
load(file <- "results/adn/all_probes/E_GEOD_5281_all_probes.RData")
load(file <- "results/adn/all_probes/E_MEXP_2280_all_probes.RData")

# Select unique probesets
genes <- unique(c(E_GEOD_18309_pr, E_GEOD_28146_pr, E_GEOD_29652_pr, E_GEOD_4757_pr, E_GEOD_5281_pr, E_MEXP_2280_pr))
length(genes)

print("Performing analysis and writing results.")
# Paths to the coexpression results in each dataset
pathRdata_E_GEOD_18309 <- "results/adn/all_probes/rdata/E_GEOD_18309/"
pathRdata_E_GEOD_28146 <- "results/adn/all_probes/rdata/E_GEOD_28146/"
pathRdata_E_GEOD_29652 <- "results/adn/all_probes/rdata/E_GEOD_29652/"
pathRdata_E_GEOD_4757 <- "results/adn/all_probes/rdata/E_GEOD_4757/"
pathRdata_E_GEOD_5281 <- "results/adn/all_probes/rdata/E_GEOD_5281/"
pathRdata_E_MEXP_2280 <- "results/adn/all_probes/rdata/E_MEXP_2280/"

library(foreach); library(doMC); cores=10 ; registerDoMC(cores);
foreach(i = 1:length(genes)) %dopar%{
gene <- genes[i]
	filedata <- sprintf("%s.RData",gene);
        pathdata_E_GEOD_18309 <- file.path(pathRdata_E_GEOD_18309, filedata);
	pathdata_E_GEOD_28146 <- file.path(pathRdata_E_GEOD_28146, filedata);
	pathdata_E_GEOD_29652 <- file.path(pathRdata_E_GEOD_29652, filedata);
	pathdata_E_GEOD_4757 <- file.path(pathRdata_E_GEOD_4757, filedata);
	pathdata_E_GEOD_5281 <- file.path(pathRdata_E_GEOD_5281, filedata);
	pathdata_E_MEXP_2280 <- file.path(pathRdata_E_MEXP_2280, filedata);


 	if ((gene%in%E_GEOD_18309_pr)==T){
        load(file=pathdata_E_GEOD_18309)
	cor1gds_E_GEOD_18309 <- cbind(E_GEOD_18309_pr,cor1gds)
	list1 <- cor1gds_E_GEOD_18309[order(-as.numeric(cor1gds_E_GEOD_18309[,2])),1] 
	 } else{list1=0} 

	if((gene%in%E_GEOD_28146_pr)==T) {
	load(file=pathdata_E_GEOD_28146)
        cor1gds_E_GEOD_28146 <- cbind(E_GEOD_28146_pr,cor1gds)
        list2 <- cor1gds_E_GEOD_28146[order(-as.numeric(cor1gds_E_GEOD_28146[,2])),1]
	}else{list2=0}
 
        if((gene%in%E_GEOD_29652_pr)==T) {
        load(file=pathdata_E_GEOD_29652)
        cor1gds_E_GEOD_29652 <- cbind(E_GEOD_29652_pr,cor1gds)
        list3 <- cor1gds_E_GEOD_29652[order(-as.numeric(cor1gds_E_GEOD_29652[,2])),1]
        }else{list3=0}

	if((gene%in%E_GEOD_4757_pr)==T) {
        load(file=pathdata_E_GEOD_4757)
        cor1gds_E_GEOD_4757 <- cbind(E_GEOD_4757_pr,cor1gds)
        list4 <- cor1gds_E_GEOD_4757[order(-as.numeric(cor1gds_E_GEOD_4757[,2])),1]
       	}else{list4=0}

	if((gene%in%E_GEOD_5281_pr)==T) {
        load(file=pathdata_E_GEOD_5281)
        cor1gds_E_GEOD_5281 <- cbind(E_GEOD_5281_pr,cor1gds)
        list5 <- cor1gds_E_GEOD_5281[order(-as.numeric(cor1gds_E_GEOD_5281[,2])),1]
        }else {list5=0}
	
        if((gene%in%E_MEXP_2280_pr)==T) {
        load(file=pathdata_E_MEXP_2280)
        cor1gds_E_MEXP_2280 <- cbind(E_MEXP_2280_pr,cor1gds)
        list6 <- cor1gds_E_MEXP_2280[order(-as.numeric(cor1gds_E_MEXP_2280[,2])),1]
        }else {list6=0}

# Create the list of ranked probestes lists
glist <- list(list1, list2, list3, list4, list5, list6)
glist <- glist[!glist%in%0]

# Aggregate the ranks
library(RobustRankAggreg)
r <- rankMatrix(glist)
ar_gene <- aggregateRanks(rmat = r, method = "RRA")

# Correction for multiple testing
ar_gene$adj.pval <- p.adjust(ar_gene$Score, method = "fdr")
head(ar_gene)
filename <- sprintf("%s.txt",gene);
filedata <- sprintf("%s.RData",gene);

# Path name to save the results
pathdata_scores <- "results/adn/all_probes/scores/rdata/"
pathdata <- file.path(pathdata_scores, filedata);
save(ar_gene,file=pathdata)

} 




