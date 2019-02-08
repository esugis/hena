library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(igraph)
library(h2o)
library(randomForest)
library(ggplot2)
library(ModelMetrics)
path = "../datasets/genes_data/"

#---------------------------------#
# Functions
#---------------------------------#
# autoencoder function
autoencoder_fn <- function(data){
  dt_pos <- filter(data, label=='disease')
  dt_unknown <- filter(data, label == 'unknown')
  dt_neg <- filter(data, label == 'non-disease')
  
  h2o.init(nthreads = -1)
  dt_ho_p <- as.h2o(dt_pos)
  dt_ho_n <- as.h2o(dt_neg)
  dt_ho_u <- as.h2o(dt_unknown)
  
  features <- setdiff(colnames(dt_ho_n), c('ensg', 'label'))
  
  splits_neg <- h2o.splitFrame(dt_ho_n, 
                               ratio = c(0.5), 
                               seed = 5622)
  
  train_neg  <- splits_neg[[1]]
  test_neg <- splits_neg[[2]]
  
  model_neg <- h2o.deeplearning(x = features,
                                training_frame = train_neg,
                                model_id = "model_neg",
                                autoencoder = TRUE,
                                reproducible = TRUE, #slow - turn off for real problems
                                ignore_const_cols = FALSE,
                                seed = 562,
                                hidden = c(100, 10, 100), 
                                epochs = 800,
                                activation = "Tanh") 
  print(paste("train mse on negative cases is ", h2o.mse(model_neg), sep=''))
  anomaly_neg <- h2o.anomaly(model_neg, test_neg) %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    mutate(ensg = as.vector(test_neg[, 'ensg']), 
           node_type = as.vector(test_neg[, 'label']))
  
  anomaly_uknown <- h2o.anomaly(model_neg, dt_ho_u) %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    mutate(ensg = as.vector(dt_ho_u[, 'ensg']), 
           node_type = as.vector(dt_ho_u[, 'label']))
  
  anomaly_pos <- h2o.anomaly(model_neg, dt_ho_p) %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    mutate(ensg = as.vector(dt_ho_p[, 'ensg']), 
           node_type = as.vector(dt_ho_p[, 'label']))
  
  anomaly <- bind_rows(anomaly_neg, anomaly_uknown, anomaly_pos)
  return(anomaly)
}
# random forest helper function
random_forest_helper <- function(data){
  set.seed(534976)
  known_cases <- filter(data, label != 'unknown')
  known_cases$label <- ifelse(known_cases$label=='disease', 1, 0)
  unknown_cases <- filter(data, label == 'unknown')
  idx_train <- sample(nrow(known_cases), 0.5*nrow(known_cases))
  train <- known_cases[idx_train,]
  test <- known_cases[-idx_train,]
  mrf <- randomForest(data=train[,-1], as.factor(label)~., importance=T, do.trace=T, na.action=na.omit)
  test$rf_pred <- predict(mrf, newdata=test[,-1], type='prob')[,2]
  test$rf_pred_b <- ifelse(test$rf_pred > 0.5, 1, 0)
  unknown_cases$rf_pred <- predict(mrf, newdata=unknown_cases[,-1], type='prob')[,2]
  data.name <- deparse(substitute(data))
  
  cm <- table(real=test$label, predicted=test$rf_pred_b)
  auc <- auc(test$label, test$rf_pred)
  prec <- ppv(test$label, test$rf_pred, cutoff=0.5)
  rec <- recall(test$label, test$rf_pred, cutoff=0.5)
  f1_score <- f1Score(test$label, test$rf_pred, cutoff=0.5)
  
  print(paste("Results for", data.name, sep=' '))
  print("Confusion matrix on a test set:")
  print(cm)
  print(paste("AUC:", auc, sep=' '))
  print(paste("precision:",prec , sep=' '))
  print(paste("recall:", rec, sep=' '))
  print(paste("f1:", f1_score, sep=' '))
  return(list(test, unknown_cases, c(data.name, cm, auc, prec, rec, f1_score)))
}
label_class_conversion <- function(data){
  data$class <- ifelse(data$label=='disease', 1, ifelse(data$label=='non-disease', 0, 2))
  data$label <- NULL
  return(data)
}
#---------------------------------#
# load data
#_________________________________#
dt_interactions <- fread(paste(path, "interactions.csv", sep=''), header=T)
dt_nodes <- fread(paste(path, "node_attributes.csv", sep=''), header=T)
nodes_hops <- fread(paste(path, "nodes_hops.csv", sep=''), header=T)
nodes_hops[is.na(nodes_hops)] <- 0
#---------------------------------#
# embeddings to features
#_________________________________#
flist <- list.files(path = path, pattern="graphsage_embs")
embs_list <- list()
for(file in flist){
  dt <- fread(paste(path, file, sep=''), header=T)
  dt <- dt[,-1]
  name_comp <- str_split(file, pattern='_')[[1]]
  edge_type <- name_comp[1]
  colnames(dt) <-  c(paste(edge_type, colnames(dt)[-c(257)], sep="_"), colnames(dt)[c(257)]) 
  data_name <- paste(edge_type, "embs", sep="_")
  assign(data_name, dt)
  embs_list[[data_name]] <- dt
}
dt_embs <- embs_list %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="ensg"), .)
dt_embs[is.na(dt_embs)] <- 0

#---------------------------------#
# dataset with all fetutes
#_________________________________#
nodes_hops$label <- NULL
dt_all <- left_join(dt_nodes, nodes_hops, by= "ensg") %>%
  left_join(dt_embs, by="ensg") 

# reducing graph and nodes to those that have non-empty biological features
# reducing graph and nodes to those that have non-empty biological features
#dt_reduced <- dt_all
#dt_reduced$is_empty <- rowSums(dt_reduced[,c(3:234)])
#dt_reduced <- filter(dt_reduced, is_empty!=0)
#dt_reduced$is_empty <- NULL
dt_reduced <- dt_all
#---------------------------------#
# 3 sets of features
#_________________________________#
# rf requires column names to be of a proper format
colnames(dt_reduced)[4:11] <- paste("brain",colnames(dt_reduced)[4:11], sep='_')
colnames(dt_reduced) <- str_replace_all(colnames(dt_reduced),'- | ','_')
colnames(dt_reduced) <- str_replace_all(colnames(dt_reduced),'-','_')

dt_biological <- dt_reduced[,c(1:234)] 
dt_graph <- dt_reduced %>% select(ensg, label, contains("epi"), contains('ppi'), contains('coexp'))
# all - dt_reduced 

write.table(label_class_conversion(dt_reduced), paste(path, "all_features.csv", sep=''), sep=',', col.names=T, row.names = F)
write.table(label_class_conversion(dt_biological), paste(path, "biological_features.csv", sep=''), sep=',', col.names=T, row.names = F)
write.table(label_class_conversion(dt_graph), paste(path, "graph_features.csv", sep=''), sep=',', col.names=T, row.names = F)

#---------------------------------#
# autoencoder for 3 feature sets
#_________________________________#
anomaly_biological <- autoencoder_fn(data=dt_biological) %>% mutate(feature_set="biological")
anomaly_graph <- autoencoder_fn(data=dt_graph) %>% mutate(feature_set="graph-related")
anomaly_all <- autoencoder_fn(data=dt_reduced) %>% mutate(feature_set="all")

anomaly <- bind_rows(anomaly_biological, anomaly_graph, anomaly_all)

mean_mse <- anomaly %>%
  group_by(feature_set, node_type) %>%
  summarise(mean = mean(Reconstruction.MSE), 
            std=sd(Reconstruction.MSE), 
            count=n()
            )
print("Autoencoder mse for all feature sets:")
print(mean_mse)
p_autoencoder <- ggplot(anomaly, aes(x=Reconstruction.MSE, fill=node_type)) + 
  geom_density(position="dodge", alpha=0.3) +
  theme_bw(base_size=18) + facet_wrap(~feature_set, scales='free')
p_autoencoder

png(paste(path, "reconstruction_mse.png", sep=''), width = 500, height=700)
print(p_autoencoder)
dev.off()

#---------------------------------#
# random forest
#_________________________________#
rf_biological <- random_forest_helper(data=dt_biological)
rf_graph <- random_forest_helper(data=dt_graph)
rf_all <- random_forest_helper(data=dt_reduced)

# save the unknown prediction for graph features
predictions_rf <- rf_graph[[2]] %>%
  select(ensg, label, rf_pred)

write.table(predictions_rf, paste(path, "rf_predictions_graph_set", sep=''), 
              sep=',', col.names = T, row.names = F)
