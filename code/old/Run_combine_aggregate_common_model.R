### Reduced model finalized information
### Combine different runs and get best model
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "reshape2", "gridExtra", "scales", "caret", "pROC"))


# Set up relevent environment variables
run_info_list <- list()
run_predictions <- list()
best_tune <- list()
probs_predictions <- list()
n <- 100

best_model_data <- as.data.frame(matrix(nrow = 100, ncol = 6))

for(i in 1:n){
  
  load(paste("exploratory/Commonfeatures_RF_model_", i, ".RData", sep=""))
  
  if(i == 1){
    write.csv(eighty_twenty_splits, "results/tables/common_test_data_splits.csv", 
              row.names = F)
    write.csv(test_data_imps, "results/tables/common_test_tune_data.csv", 
              row.names = F)
  }
  
  probs_predictions[[paste("run_", i, sep = "")]] <- 
    predict(test_tune_list[[paste("data_split", i, sep = "")]], 
            test_test_data, type = 'prob')
  
  rm(list = setdiff(ls(), 
                    c("test_tune_list", "test_predictions", "best_tune", 
                      "best_model_data", "imp_vars_list", "run_info_list", 
                      "run_predictions", "n", "i", "probs_predictions", "all_runs_list")))
  
  best_tune[paste("run_", i, sep = "")] <- test_tune_list[[paste(
    "data_split", i, sep = "")]]$bestTune
  
  run_info_list[[paste("run_", i, sep = "")]] <- 
    test_tune_list[[paste("data_split", i, sep = "")]]$results 
  
  run_predictions[[paste("run_", i, sep = "")]] <- test_predictions[[paste(
    "data_split", i, sep = "")]]
  
  best_model_data[i, ] <- filter(run_info_list[[i]], 
                                 mtry == best_tune[[i]]) %>% select(-mtry)
  colnames(best_model_data) <- colnames(select(run_info_list[[i]], -mtry))
  rownames(best_model_data)[i] <- paste("run_", i, sep = "")
  
  rm(test_predictions, test_tune_list)
  
}

# Write out ROC summary table

write.csv(
  mutate(best_model_data, run = rownames(best_model_data), 
         best_mtry = t(as.data.frame.list(best_tune))), 
  "results/tables/common_ROC_model_summary.csv", row.names = F)


# Get Ranges of 100 10-fold 20 times CV data (worse, best)
actual_data <- read.csv("results/tables/common_test_tune_data.csv", header = T)

data_splits <- read.csv("results/tables/common_test_data_splits.csv", 
                        header = T, stringsAsFactors = F)

best_run <- as.numeric(strsplit((
  mutate(best_model_data, run = rownames(best_model_data)) %>% 
    filter(ROC == max(best_model_data$ROC)) %>% 
    select(run))[1,], "_")[[1]][2])

worse_run <- as.numeric(strsplit((mutate(best_model_data, 
                                         run = rownames(best_model_data)) %>% 
                                    filter(ROC == min(best_model_data$ROC)) %>% 
                                    select(run))[1,], "_")[[1]][2])

best_split <- data_splits[, best_run]
worse_split <- data_splits[, worse_run]

roc_data_list <- list(
  best_roc = roc(
    actual_data[-best_split, ]$lesion ~ 
      probs_predictions[[best_run]][, "Yes"]), 
  worse_roc = roc(actual_data[-worse_split, ]$lesion ~ 
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
          "results/tables/common_test_data_roc.csv", row.names = F)


# Create AUC data table

auc_data_table <- as.data.frame(matrix(
  nrow = 2, ncol = 4, dimnames = list(
    nrows = c("best", "worse"), ncols = c("AUC", "ROC_cv", "Sens_cv", "Spec_cv"))))

auc_data_table[, "AUC"] <- c( 
  roc_data_list[["best_roc"]]$auc, 
  roc_data_list[["worse_roc"]]$auc)

auc_data_table[, "ROC_cv"] <- c( 
  best_model_data[paste("run_", best_run, sep = ""), "ROC"], 
  best_model_data[paste("run_", worse_run, sep = ""), "ROC"])

auc_data_table[, "Sens_cv"] <- c( 
  best_model_data[paste("run_", best_run, sep = ""), "Sens"], 
  best_model_data[paste("run_", worse_run, sep = ""), "Sens"])

auc_data_table[, "Spec_cv"] <- c( 
  best_model_data[paste("run_", best_run, sep = ""), "Spec"], 
  best_model_data[paste("run_", worse_run, sep = ""), "Spec"])


write.csv(auc_data_table, "results/tables/common_auc_summary.csv")


