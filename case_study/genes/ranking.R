
library(data.table)
fr_rank <- fread(file = "../datasets/genes_data/rf_predictions_graph_set.csv")
hinsage_rank <- fread(file = "../datasets/genes_data/hinsage_graph_features.csv")
hinsage_rank <- hinsage_rank[, -c(1,2)]

library(dplyr)
hinsage_rank<- arrange(hinsage_rank, desc(pred))
colnames(hinsage_rank)[1] <- 'hinsage_pred' 
colnames(hinsage_rank)[2] <- 'ensg'
unknown_cases <- left_join(fr_rank, hinsage_rank, by='ensg')

unknown_cases <- arrange(unknown_cases, desc(hinsage_pred))

unknown_cases <- unknown_cases %>% ungroup() %>%
  mutate(hinsage_rank=rank(desc(hinsage_pred))) %>% 
  mutate(rf_rank=rank(desc(rf_pred)))

write.csv(unknown_cases, file="../datasets/genes_data/predictions_ranked.csv")
