library(data.table)
library(tidyverse)

#---------------------------------#
# load data
#_________________________________#
# data
path_to_data = "../datasets/genes_data/"
dt_nodes <- fread(paste(path_to_data, "node_attributes.csv", sep=''), header=T)

#---------------------------------#
# load hops data
#_________________________________#
flist <- list.files(path = path_to_data, pattern="hops_")
hops_list <- list()
for(file in flist){
  dt <- fread(paste(path_to_data, file, sep=''), header=F)
  name_comp <- str_split(file, pattern='_')[[1]]
  depth <- str_split(name_comp[3], pattern='.csv')[[1]][1]
  edge_type = name_comp[2]
  hop_count <-  paste(edge_type, depth, 'hop_count', sep='_')
  hop_total <-  paste(edge_type, depth, 'hop_total', sep='_')
  data_name <- paste(edge_type, "hop", depth, sep="_")
  colnames(dt) <- c('ensg', hop_count, hop_total)
  dt <- mutate(dt, ratio = get(hop_count)/get(hop_total))
  colnames(dt)[4] <- paste(edge_type, depth, "hop_ratio", sep='_')
  assign(data_name, dt)
  hops_list[[data_name]] <- dt
}
dt_nodes_hops <- hops_list %>%
  Reduce(function(dtf1, dtf2) full_join(dtf1,dtf2,by="ensg"), .)

dt_nodes_hops <- left_join(dt_nodes, dt_nodes_hops, by=c('ensg'))
dt_nodes_hops[is.na(dt_nodes_hops)] <- 0

dt_nodes_hops$label <- ifelse(dt_nodes_hops$label=='disease', 1, 
                              ifelse(dt_nodes_hops$label=='non-disease', 0, 2))

write.table(dt_nodes_hops[,c(1,2,235:252)], paste(path_to_data, "nodes_hops.csv", sep=''), 
            col.names = T, row.names = F, sep=',')


