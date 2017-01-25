### Combine models from the Training runs
### Specifically identify mean and SD of model 
### and pull out most important OTUs used within it
## Marc Sze


### Note:
## Random Forest: from the R package: “For each tree, the prediction accuracy 
## on the out-of-bag portion of the data is recorded. Then the same is done after permuting 
## each predictor variable. The difference between the two accuracies are then averaged over all trees, 
## and normalized by the standard error. For regression, the MSE is computed on the out-of-bag data for each 
## tree, and then the same computed after permuting a variable. The differences are averaged and normalized 
## by the standard error. If the standard error is equal to 0 for a variable, the division is not done.”

## Cutoff of 0.5 (default was used for this) for RF model

###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "reshape2", "scales", "caret", "pROC"))

# Set up relevent environment variables
imp_vars_list <- list()
run_info_list <- list()
run_predictions <- list()
best_tune <- list()
probs_predictions <- list()
n <- 100

best_model_data <- as.data.frame(matrix(nrow = 100, ncol = 6))

for(i in 1:n){
  
  load(paste("exploratory/RF_model_", i, ".RData", sep=""))
  
  if(i == 1){
    write.csv(eighty_twenty_splits, "results/tables/test_data_splits.csv", 
      row.names = F)
    write.csv(test_data, "results/tables/test_tune_data.csv", 
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
  
  rm(test_predictions, test_tune_list)
  
}

# Write out ROC summary table

write.csv(
  mutate(best_model_data, run = rownames(best_model_data), 
    best_mtry = t(as.data.frame.list(best_tune))), 
  "results/tables/ROC_model_summary.csv", row.names = F)

# Calculate number of times an OTU is in the top 10% of overall importance

OTU_appearance_table <- as.data.frame(data_frame(
  Variable = imp_vars_list[["run_1"]]$Variable) %>% 
mutate(total_appearance = 0))

rownames(OTU_appearance_table) <- OTU_appearance_table$Variable

for(j in 1:length(imp_vars_list)){
  
  tempVars <- imp_vars_list[[j]][c(1:round(length(
    rownames(imp_vars_list[[j]]))*0.10)), ][, "Variable"]
  
  for(i in tempVars){
    
    OTU_appearance_table[i, "total_appearance"] <- 
    OTU_appearance_table[i, "total_appearance"] + 1
  }
}

OTU_appearance_table <- arrange(OTU_appearance_table, 
  desc(total_appearance))

# Keep Those over 50% of the total 100 runs of 80/20 splits
OTU_appearance_table <- filter(OTU_appearance_table, total_appearance > 50)

# Write out the important variables to a table
write.csv(OTU_appearance_table, 
  "results/tables/rf_wCV_imp_vars_summary.csv", row.names = F)

# Collect the mean and SD for the MDA of the most important variables

top_vars_MDA <- lapply(imp_vars_list, function(x) 
  x[order(x[, "Variable"]), ] %>% filter(Variable %in% OTU_appearance_table$Variable))

top_vars_MDA_by_run <- as.data.frame(matrix(nrow = length(OTU_appearance_table$Variable), 
                                            ncol = length(imp_vars_list), 
                                            dimnames = list(
                                              nrow = top_vars_MDA[["run_1"]]$Variable[
                                                order(top_vars_MDA[["run_1"]]$Overall, decreasing = T)], 
                                              ncol = paste("run_", seq(1:100), sep = ""))))

for(i in 1:length(top_vars_MDA_by_run)){
  
  tempData <- top_vars_MDA[[i]]
  rownames(tempData) <- tempData$Variable
  top_vars_MDA_by_run[, i] <- tempData[rownames(top_vars_MDA_by_run), "Overall"]
  rm(tempData)
}

# "1" pulls the value of mean or sd from the data frame
MDA_vars_summary <- cbind(
  mean_MDA = t(summarise_each(as.data.frame(t(top_vars_MDA_by_run)), funs(mean)))[, 1], 
  sd_MDA = t(summarise_each(as.data.frame(t(top_vars_MDA_by_run)), funs(sd)))[, 1], 
  variable = rownames(top_vars_MDA_by_run))

write.csv(MDA_vars_summary[order(MDA_vars_summary[, "mean_MDA"], decreasing = TRUE), ], 
          "results/tables/lesion_model_top_vars_MDA_Summary.csv", row.names = F)

lesion_model_top_vars_MDA_full_data <- 
  mutate(top_vars_MDA_by_run, variables = rownames(top_vars_MDA_by_run)) %>% 
  melt(id = c("variables"))

write.csv(lesion_model_top_vars_MDA_full_data, 
          "results/tables/lesion_model_top_vars_MDA_full_data.csv", row.names = F)
  
# Pull out middle(ish) model from runs and use that in the prediction of lesion in 

middle_run <- as.numeric(
  strsplit((best_model_data[order(desc(best_model_data$ROC)), ] %>% 
    mutate(run = rownames(.)) %>% 
    slice(length(rownames(best_model_data))/2) %>% 
    select(run))[1,], "_")[[1]][2])

# Get Ranges of 100 10-fold 20 times CV data (worse, middle, best)
actual_data <- read.csv("results/tables/full_test_data.csv", header = T, row.names = 1)

data_splits <- read.csv("results/tables/test_data_splits.csv", 
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
middle_split <- data_splits[, middle_run]

roc_data_list <- list(
  best_roc = roc(
    actual_data[-best_split, ]$lesion ~ 
    probs_predictions[[best_run]][, "Yes"]), 
  middle_roc = roc(actual_data[-middle_split, ]$lesion ~ 
    probs_predictions[[middle_run]][, "Yes"]), 
  worse_roc = roc(actual_data[-worse_split, ]$lesion ~ 
    probs_predictions[[worse_run]][, "Yes"]))

# Build data table for figure 3
test_roc_data <- cbind(
  sensitivities = c(roc_data_list[["best_roc"]]$sensitivities, 
    roc_data_list[["middle_roc"]]$sensitivities, 
    roc_data_list[["worse_roc"]]$sensitivities), 
  specificities = c(roc_data_list[["best_roc"]]$specificities, 
    roc_data_list[["middle_roc"]]$specificities, 
    roc_data_list[["worse_roc"]]$specificities), 
  run = c(rep("best_roc", 
    length(roc_data_list[["best_roc"]]$sensitivities)), 
  rep("middle_roc", length(roc_data_list[["middle_roc"]]$sensitivities)), 
  rep("worse_roc", length(roc_data_list[["worse_roc"]]$sensitivities))))

write.csv(test_roc_data, 
  "results/tables/test_data_roc.csv", row.names = F)

# Create AUC data table for figure 3

auc_data_table <- as.data.frame(matrix(
  nrow = 3, ncol = 4, dimnames = list(
    nrows = c("best", "middle", "worse"), ncols = c("AUC", "ROC_cv", "Sens_cv", "Spec_cv"))))

auc_data_table[, "AUC"] <- c( 
  roc_data_list[["best_roc"]]$auc, 
  roc_data_list[["middle_roc"]]$auc, 
  roc_data_list[["worse_roc"]]$auc)

auc_data_table[, "ROC_cv"] <- c( 
  best_model_data[paste("run_", best_run, sep = ""), "ROC"], 
  best_model_data[paste("run_", middle_run, sep = ""), "ROC"], 
  best_model_data[paste("run_", worse_run, sep = ""), "ROC"])

auc_data_table[, "Sens_cv"] <- c( 
  best_model_data[paste("run_", best_run, sep = ""), "Sens"], 
  best_model_data[paste("run_", middle_run, sep = ""), "Sens"], 
  best_model_data[paste("run_", worse_run, sep = ""), "Sens"])

auc_data_table[, "Spec_cv"] <- c( 
  best_model_data[paste("run_", best_run, sep = ""), "Spec"], 
  best_model_data[paste("run_", middle_run, sep = ""), "Spec"], 
  best_model_data[paste("run_", worse_run, sep = ""), "Spec"])


write.csv(auc_data_table, "results/tables/auc_summary.csv")


# Keep everything but roc_data_list in memory
rm(list = setdiff(ls(), "roc_data_list"))
save.image("exploratory/rocs.RData")












