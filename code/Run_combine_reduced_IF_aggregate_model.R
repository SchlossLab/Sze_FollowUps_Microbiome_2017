###Aggregate all the data from Reduced IF model
### Get the best and worst models and correct mtry for new model
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "ggplot2", "reshape2", "gridExtra", "scales", 
           "wesanderson", "caret", "pROC"))


# Load needed data
load("exploratory/IF_reduced_RF_model_Imp_OTU.RData")

# Set up relevent environment variables
imp_vars_list <- list()
run_info_list <- list()
run_predictions <- list()
best_tune <- list()
probs_predictions <- list()
n <- 100

best_model_data <- as.data.frame(matrix(nrow = 100, ncol = 6))


for(i in 1:n){
  
  if(i == 1){
    write.csv(eighty_twenty_splits, "results/tables/reduced_IF_test_data_splits.csv", 
              row.names = F)
    write.csv(test_data_imps, "results/tables/reduced_IF_test_tune_data.csv", 
              row.names = F)
  }
  
  probs_predictions[[paste("run_", i, sep = "")]] <- 
    predict(test_tune_list[[paste("data_split", i, sep = "")]], 
            test_test_data, type = 'prob')
  
  best_tune[paste("run_", i, sep = "")] <- test_tune_list[[paste(
    "data_split", i, sep = "")]]$bestTune
  
  run_info_list[[paste("run_", i, sep = "")]] <- 
    test_tune_list[[paste("data_split", i, sep = "")]]$results 
  
  imp_vars_list[[paste("run_", i, sep = "")]] <- 
    varImp(test_tune_list[[paste("data_split", i, sep = "")]], 
           scale = FALSE)$importance %>% 
    mutate(Variable = rownames(.)) %>% arrange(desc(Overall))
  
  
  run_predictions[[paste("run_", i, sep = "")]] <- test_predictions[[paste(
    "data_split", i, sep = "")]]
  
  best_model_data[i, ] <- filter(run_info_list[[i]], 
                                 mtry == best_tune[[i]]) %>% select(-mtry)
  colnames(best_model_data) <- colnames(select(run_info_list[[i]], -mtry))
  rownames(best_model_data)[i] <- paste("run_", i, sep = "")
  
}

write.csv(
  mutate(best_model_data, run = rownames(best_model_data), 
         best_mtry = t(as.data.frame.list(best_tune))), 
  "results/tables/reduced_IF_ROC_model_summary.csv", row.names = F)


# Get Ranges of 100 10-fold 20 times CV data (worse, best)
best_run <- as.numeric(strsplit((
  mutate(best_model_data, run = rownames(best_model_data)) %>% 
    filter(ROC == max(best_model_data$ROC)) %>% 
    select(run))[1,], "_")[[1]][2])

worse_run <- as.numeric(strsplit((mutate(best_model_data, 
                                         run = rownames(best_model_data)) %>% 
                                    filter(ROC == min(best_model_data$ROC)) %>% 
                                    select(run))[1,], "_")[[1]][2])

best_split <- eighty_twenty_splits[, best_run]
worse_split <- eighty_twenty_splits[, worse_run]

roc_data_list <- list(
  best_roc = roc(
    test_data[-best_split, ]$lesion ~ 
      probs_predictions[[best_run]][, "Yes"]), 
  worse_roc = roc(test_data[-worse_split, ]$lesion ~ 
                    probs_predictions[[worse_run]][, "Yes"]))

# Get sensitivity and specificity for test data (best and worse)
test_roc_data <- cbind(
  sensitivities = c(roc_data_list[["best_roc"]]$sensitivities, 
                    roc_data_list[["worse_roc"]]$sensitivities), 
  specificities = c(roc_data_list[["best_roc"]]$specificities, 
                    roc_data_list[["worse_roc"]]$specificities), 
  run = c(
    rep("best_roc", length(roc_data_list[["best_roc"]]$sensitivities)), 
    rep("worse_roc", length(roc_data_list[["worse_roc"]]$sensitivities))))

write.csv(test_roc_data, 
          "results/tables/reduced_IF_test_data_roc.csv", row.names = F)






