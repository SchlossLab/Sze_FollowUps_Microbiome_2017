### Create IF model using the important variables only
### Does this model perform similarily to the full model
## Marc Sze


#Load needed libraries
source('code/functions.R')

loadLibs(c("randomForest", "dplyr", "scales", "caret"))

#Load needed data
load("exploratory/srn_randomization_treatment_model.RData")

rm(good_metaf, good_metaf_select, shared, test_data, 
   fitControl, i, nzv)


#################################################################################
#                                                                               #
#                                                                               #
#               Model Training and Parameter Tuning                             #
#                                                                               #
#################################################################################



#Create Overall specifications for model tuning
fitControl <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  ## repeated twenty times
  p = 0.8, 
  classProbs = TRUE, 
  savePredictions = TRUE, 
  summaryFunction = twoClassSummary)

set.seed(3457)

#Set up lists to store the data
runs_list <- list()
avg_auc <- c()
sd_auc <- c()
min_auc <- c()
max_auc <- c()


for(i in 1:100){
  
  # Get the MDA counts from most to least
  test_imp <- varImp(test_tune_list[[i]], scale = F)$importance %>% 
    mutate(otu = rownames(.)) %>% 
    arrange(desc(Overall))
  
  #Create data table with only reduced features (impvars only)
  vars_to_keep <- slice(test_imp, 1:10)$otu
  test_data_imps <- select(stored_data[[i]], lesion, one_of(vars_to_keep))
  
  train_test_data <- test_data_imps
  stored_aucs <- c()
  
  for(j in 1:100){
    
    #Train the model
    runs_list[[paste("data_split", j, sep = "")]] <- 
      train(lesion ~ ., data = train_test_data, 
            method = "rf", 
            ntree = 100, 
            trControl = fitControl, 
            metric = "ROC", 
            na.action = na.omit, 
            verbose = FALSE)
    
    stored_aucs <- c(stored_aucs, ifelse(runs_list$results$ROC < 0.5, 
                                         invisible(1-runs_list$results$ROC), 
                                         invisible(runs_list$results$ROC)))
    
  }
  
  sd_auc <- c(sd_auc, sd(stored_aucs))
  avg_auc <- c(avg_auc, mean(stored_aucs, na.rm = TRUE))
  min_auc <- c(min_auc, min(stored_aucs))
  max_auc <- c(max_auc, max(stored_aucs))
  
  
}

final_results <- cbind(runs = c(1:100), avg_auc = avg_auc, min_auc = min_auc, max_auc = max_auc)


# Save image with data and relevant parameters
write.csv(final_results, "data/process/tables/srn_rand_treat_reduced_summary.csv", row.names = F)



