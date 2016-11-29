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
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", "wesanderson", "caret"))


# Set up relevent environment variables
imp_vars_list <- list()
run_info_list <- list()
run_predictions <- list()
best_tune <- list()
probs_predictions <- list()
model_cutoff <- list()
n <- 100

best_model_data <- as.data.frame(matrix(nrow = 100, ncol = 6))

for(i in 1:n){
  
  load(paste("exploratory/RF_model_", i, ".RData", sep=""))
  
  if(i == 1){
    write.csv(eighty_twenty_splits, "results/tables/test_data_splits.csv", row.names = F)
  }
  
  probs_predictions[[paste("run_", i, sep = "")]] <- 
    predict(test_tune_list[[paste("data_split", i, sep = "")]], test_test_data, type = 'prob')
  
  rm(list = setdiff(ls(), c("test_tune_list", "test_predictions", "best_tune", "best_model_data", 
                            "imp_vars_list", "run_info_list", "run_predictions", "n", "i", 
                            "probs_predictions")))
  
  
  best_tune[paste("run_", i, sep = "")] <- test_tune_list[[paste("data_split", i, sep = "")]]$bestTune
  
  run_info_list[[paste("run_", i, sep = "")]] <- 
    test_tune_list[[paste("data_split", i, sep = "")]]$results 
  
  imp_vars_list[[paste("run_", i, sep = "")]] <- 
    varImp(test_tune_list[[paste("data_split", i, sep = "")]], scale = FALSE)$importance %>% 
    mutate(Variable = rownames(.)) %>% arrange(desc(Overall))
  
  
  run_predictions[[paste("run_", i, sep = "")]] <- test_predictions[[paste("data_split", i, sep = "")]]
  
  best_model_data[i, ] <- filter(run_info_list[[i]], mtry == best_tune[[i]]) %>% select(-mtry)
  colnames(best_model_data) <- colnames(select(run_info_list[[i]], -mtry))
  rownames(best_model_data)[i] <- paste("run_", i, sep = "")
  
  rm(test_predictions, test_tune_list)
  
}

# Write out ROC summary table

write.csv(
  mutate(best_model_data, run = rownames(best_model_data), best_mtry = t(as.data.frame.list(best_tune))), 
  "results/tables/ROC_model_summary.csv", row.names = F)

# Calculate number of times an OTU is in the top 10% of overall importance

OTU_appearance_table <- as.data.frame(data_frame(Variable = imp_vars_list[["run_1"]]$Variable) %>% 
                                        mutate(total_appearance = 0))
rownames(OTU_appearance_table) <- OTU_appearance_table$Variable

for(j in 1:length(imp_vars_list)){
  
  tempVars <- imp_vars_list[[j]][c(1:round(length(rownames(imp_vars_list[[j]]))*0.10)), ][, "Variable"]
  
  for(i in tempVars){
    
    OTU_appearance_table[i, "total_appearance"] <- OTU_appearance_table[i, "total_appearance"] + 1
  }
}

OTU_appearance_table <- arrange(OTU_appearance_table, desc(total_appearance))

# Keep Those over 50% of the total 100 runs of 80/20 splits
OTU_appearance_table <- filter(OTU_appearance_table, total_appearance > 50)

# Write out the important variables to a table
write.csv(OTU_appearance_table, "results/tables/rf_wCV_imp_vars_summary.csv", row.names = F)

# Pull out best model and middle(ish) model from runs and use that in the prediction of lesion in 

best_run <- as.numeric(strsplit((mutate(best_model_data, run = rownames(best_model_data)) %>% 
  filter(ROC == max(best_model_data$ROC)) %>% select(run))[1,], "_")[[1]][2])

middle_run <- as.numeric(strsplit((best_model_data[order(desc(best_model_data$ROC)), ] %>% 
                                     mutate(run = rownames(.)) %>% 
                                     slice(length(rownames(best_model_data))/2) %>% 
                                              select(run))[1,], "_")[[1]][2])

#load in best model
load(paste("exploratory/RF_model_", middle_run, ".RData", sep=""))

#load in metadata to get IDs to select shared file by
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)

#create vector in correct order

samples <- c(good_metaf$initial, good_metaf$followUp)

#read in shared file and keep only the samples that are needed
shared <- read.delim('data/process/final.0.03.subsample.shared', header=T, sep='\t') %>% select(-label, -numOtus)
rownames(shared) <- shared$Group
shared <- shared[as.character(samples), ]

# one sample has no follow up fit so need to remove that
samplesToRemove <- filter(good_metaf, is.na(fit_followUp)) %>% select(initial, followUp)

# Update good_metaf
updated_metaf <- filter(good_metaf, initial != samplesToRemove[, "initial"])


# Remove the sample
shared <- filter(shared, Group != samplesToRemove[, "initial"], Group != samplesToRemove[, "followUp"])

# Keep only OTUs in test data
shared <- select(shared, Group, one_of(colnames(train_test_data)))


# Keep only needed things
rm(list = setdiff(ls(), c("test_tune_list", "test_predictions", "best_tune", "best_model_data", 
                          "imp_vars_list", "run_info_list", "run_predictions", "n", "updated_metaf",  
                          "probs_predictions", "OTU_appearance_table", "middle_run", "shared")))

# Load in test data set
follow_up_data <- read.csv("results/tables/follow_up_prediction_table.csv")

test3 <- cbind(select(follow_up_data, -contains("Otu0")), select(shared, -Group))


# Make predictions on samples with follow up
initial_predictions <- predict(test_tune_list[[paste("data_split", middle_run, sep = "")]], 
                               newdata = test3[1:(length(rownames(follow_up_data))/2), ], type='prob')

followup_predictions <- predict(test_tune_list[[paste("data_split", middle_run, sep = "")]], newdata = 
                                  test3[((length(rownames(follow_up_data))/2)+1):
                                          length(rownames(follow_up_data)), ], type='prob')

# Create data table needed for figure 4

probability_data_table <- cbind(No = c(initial_predictions[, "No"], followup_predictions[, "No"]), 
                                Yes = c(initial_predictions[, "Yes"], followup_predictions[, "Yes"]), 
                                sampleType = c(rep("initial", length(rownames(initial_predictions))), 
                                               rep("followup", length(rownames(followup_predictions)))), 
                                disease_free = rep(updated_metaf$Disease_Free, 2), 
                                diagnosis = rep(updated_metaf$Diagnosis, 2), 
                                Dx_Bin = rep(updated_metaf$Dx_Bin, 2), 
                                chemo = rep(updated_metaf$chemo_received, 2), 
                                rads = rep(updated_metaf$radiation_received, 2), 
                                EDRN = rep(updated_metaf$EDRN, 2))

write.csv(probability_data_table, "results/tables/follow_up_probability_summary.csv", row.names = F)

## Modify Run_Figure4 to update with the new information




### Complete Comparison R file














