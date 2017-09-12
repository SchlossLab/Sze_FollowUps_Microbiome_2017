### Reduced model finalized model generation
### get the best mtry and build final model with all data on reduced data set
## Marc Sze


#Load needed libraries
source('code/functions.R')

loadLibs(c("dplyr", "caret","scales", "wesanderson", "randomForest", "pROC"))

# Read in necessary data frames
test_data <- read.csv("data/process/tables/crc_full_test_data.csv", header = TRUE, row.names = 1)
split_data_results <- read.csv("data/process/tables/crc_ROC_model_summary.csv", header = TRUE, stringsAsFactors = F)
test_data_roc <- read.csv("data/process/tables/crc_test_data_roc.csv", header = TRUE, stringsAsFactors = F)

# Get best mtry to use
mtry_table <- table(split_data_results$mtry)
maximized_mtry <- unique(filter(split_data_results, test_auc == max(test_auc))[, "mtry"])

# Create Random Forest model

set.seed(3457)
full_model <- randomForest(lesion ~ ., data = test_data, mtry = maximized_mtry, ntree = 2000)

full_model_probs_predictions <- full_model$votes

full_model_roc <- roc(test_data$lesion ~ full_model_probs_predictions[, "Yes"])

# Add to needed data tables and write them back out
auc_data_table <- rbind(auc_data_table, full = c(full_model_roc$auc, NA, NA, NA))
write.csv(auc_data_table, "data/process/tables/crc_reduced_auc_summary.csv")
write.csv(cbind(lesion = as.character(test_data$lesion), full_model_probs_predictions), 
          "data/process/tables/crc_reduced_lesion_test_data_probs_summary.csv", row.names = F)

test_data_roc <- rbind(test_data_roc, 
                       cbind(sensitivities = full_model_roc$sensitivities, 
                             specificities = full_model_roc$specificities, 
                             run = rep("full_roc", length(full_model_roc$sensitivities))))

write.csv(test_data_roc, 
          "data/process/tables/crc_reduced_test_data_roc.csv", row.names = F)


#load in metadata to get IDs to select shared file by
good_metaf <- read.csv("data/process/mod_metadata/good_metaf_final.csv", 
                       stringsAsFactors = F, header = T)

#create vector in correct order

samples <- c(good_metaf$initial, good_metaf$followUp)

#read in shared file and keep only the samples that are needed
shared <- read.delim('data/process/final.0.03.subsample.shared', 
                     header=T, sep='\t') %>% select(-label, -numOtus)
rownames(shared) <- shared$Group
shared <- shared[as.character(samples), ]

# Keep only OTUs in test data
OTUs_to_keep <- colnames(select(test_data, -lesion))
shared <- select(shared, Group, one_of(OTUs_to_keep))

# Create test data set for follow up samples
good_metaf$followup_call[good_metaf$Disease_Free == "n"] <- "Yes"
good_metaf$followup_call[good_metaf$Disease_Free != "n"] <- "No"
test_follow_up_data <- cbind(lesion = c(rep("Yes", length(good_metaf$lesion)), good_metaf$followup_call), 
                             dx = good_metaf$dx, 
                             sampleType = c(rep("initial", length(good_metaf$lesion)), 
                                            rep("followup", length(good_metaf$lesion))), 
                             select(shared, -Group))


# Make predictions on samples with follow up
initial_predictions <- predict(full_model, 
                               newdata = (filter(test_follow_up_data, 
                                                 sampleType == "initial" & dx == "cancer") %>% 
                                            select(-dx, -sampleType)), type='prob')

followup_predictions <- predict(full_model, 
                                newdata = (filter(test_follow_up_data, 
                                                 sampleType == "followup" & dx == "cancer") %>% 
                                            select(-dx, -sampleType)), type='prob')


# Create data table needed for figure 4

probability_data_table <- cbind(
  No = c(initial_predictions[, "No"], followup_predictions[, "No"]), 
  Yes = c(initial_predictions[, "Yes"], followup_predictions[, "Yes"]), 
  sampleType = c(rep("initial", length(rownames(initial_predictions))), 
                 rep("followup", length(rownames(followup_predictions)))), 
  disease_free = rep(filter(good_metaf, dx == "cancer")[, "Disease_Free"], 2), 
  diagnosis = rep(filter(good_metaf, dx == "cancer")[, "Diagnosis"], 2), 
  Dx_Bin = rep(filter(good_metaf, dx == "cancer")[, "Dx_Bin"], 2), 
  chemo = rep(filter(good_metaf, dx == "cancer")[, "chemo_received"], 2), 
  rads = rep(filter(good_metaf, dx == "cancer")[, "radiation_received"], 2), 
  EDRN = rep(filter(good_metaf, dx == "cancer")[, "EDRN"], 2))

write.csv(probability_data_table, 
          "data/process/tables/crc_reduced_follow_up_probability_summary.csv", row.names = F)




