###Aggregate all the data from this IF model
### Get the best and worst models and correct mtry for new model
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "ggplot2", "reshape2", "gridExtra", "scales", 
           "wesanderson", "caret", "pROC"))


# Load needed data
load("exploratory/RF_model_Imp_OTU.RData")

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
    write.csv(eighty_twenty_splits, "data/process/tables/IF_test_data_splits.csv", 
              row.names = F)
    write.csv(test_data, "data/process/tables/IF_test_tune_data.csv", 
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
  "data/process/tables/IF_ROC_model_summary.csv", row.names = F)

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
          "data/process/tables/IF_test_data_roc.csv", row.names = F)

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
          "data/process/tables/IF_rf_wCV_imp_vars_summary.csv", row.names = F)

# Collect the mean and SD for the MDA of the most important variables

top_vars_MDA <- lapply(imp_vars_list, function(x) 
  x[order(x[, "Variable"]), ] %>% filter(Variable %in% OTU_appearance_table$Variable))

top_vars_MDA_by_run <- as.data.frame(matrix(nrow = length(OTU_appearance_table$Variable), 
                                            ncol = length(imp_vars_list), 
                                            dimnames = list(
                                              nrow = top_vars_MDA[["run_1"]]$Variable, 
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
          "data/process/tables/IF_model_top_vars_MDA_Summary.csv", row.names = F)


IF_model_top_vars_MDA_full_data <- 
  mutate(top_vars_MDA_by_run, variables = rownames(top_vars_MDA_by_run)) %>% 
  melt(id = c("variables"))

write.csv(IF_model_top_vars_MDA_full_data, 
          "data/process/tables/IF_model_top_vars_MDA_full_data.csv", row.names = F)







