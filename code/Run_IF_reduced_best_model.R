### Reduced IF model finalized model generation
### get the best mtry and build final model with all data on reduced data set
## Marc Sze

#Load needed libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "caret","scales", "wesanderson", "randomForest", "pROC"))

# Read in necessary data frames
test_data <- read.csv("results/tables/reduced_IF_test_tune_data.csv", header = TRUE)
split_data_results <- read.csv("results/tables/reduced_IF_ROC_model_summary.csv", header = TRUE, stringsAsFactors = F)


# Get best mtry to use
mtry_table <- table(split_data_results$best_mtry)
maximized_mtry <- as.numeric(names(mtry_table[mtry_table == max(mtry_table)]))

# Create Random Forest model

set.seed(3457)
full_model <- randomForest(lesion ~ ., data = test_data, mtry = maximized_mtry, ntree = 2000)

full_model_probs_predictions <- full_model$votes

full_model_roc <- roc(test_data$lesion ~ full_model_probs_predictions[, "Yes"])

# Write out needed data
test_data_roc <- cbind(sensitivities = full_model_roc$sensitivities, 
                       specificities = full_model_roc$specificities, 
                       run = rep("full_roc", length(full_model_roc$sensitivities)))

write.csv(test_data_roc, 
          "results/tables/reduced_IF_test_data_roc.csv", row.names = F)

#load in metadata to get IDs to select shared file by
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", 
                       stringsAsFactors = F, header = T)

#create vector in correct order

samples <- c(good_metaf$initial, good_metaf$followUp)

#read in shared file and keep only the samples that are needed
shared <- read.delim('data/process/final.0.03.subsample.shared', 
                     header=T, sep='\t') %>% select(-label, -numOtus)
rownames(shared) <- shared$Group
shared <- shared[as.character(samples), ]

# one sample has no follow up fit so need to remove that
samplesToRemove <- filter(good_metaf, is.na(fit_followUp)) %>% 
  select(initial, followUp)

# Update good_metaf
updated_metaf <- filter(
  good_metaf, initial != samplesToRemove[, "initial"])


# Remove the sample
shared <- filter(
  shared, Group != samplesToRemove[, "initial"], 
  Group != samplesToRemove[, "followUp"])

# Keep only OTUs in test data
shared <- select(shared, Group, one_of(colnames(test_data)))

# Create test data set for follow up samples
updated_metaf$followup_call[updated_metaf$Disease_Free == "n"] <- "Yes"
updated_metaf$followup_call[updated_metaf$Disease_Free != "n"] <- "No"
test_follow_up_data <- cbind(lesion = c(rep("Yes", length(updated_metaf$lesion)), updated_metaf$followup_call), 
                             fit_result = c(updated_metaf$fit_result, updated_metaf$fit_followUp), 
                             select(shared, -Group))


# Make predictions on samples with follow up
initial_predictions <- predict(full_model, 
                               newdata = test_follow_up_data[1:(length(rownames(test_follow_up_data))/2), ], 
                               type='prob')

followup_predictions <- predict(full_model, 
                                newdata = test_follow_up_data[(length(rownames(
                                  updated_metaf))+1):length(rownames(test_follow_up_data)), ], 
                                type='prob')


# Create data table needed for figure 4

probability_data_table <- cbind(
  No = c(initial_predictions[, "No"], followup_predictions[, "No"]), 
  Yes = c(initial_predictions[, "Yes"], followup_predictions[, "Yes"]), 
  sampleType = c(rep("initial", length(rownames(initial_predictions))), 
                 rep("followup", length(rownames(followup_predictions)))), 
  disease_free = rep(updated_metaf$Disease_Free, 2), 
  diagnosis = rep(updated_metaf$Diagnosis, 2), 
  Dx_Bin = rep(updated_metaf$Dx_Bin, 2), 
  chemo = rep(updated_metaf$chemo_received, 2), 
  rads = rep(updated_metaf$radiation_received, 2), 
  EDRN = rep(updated_metaf$EDRN, 2))

write.csv(probability_data_table, 
          "results/tables/reduced_IF_follow_up_probability_summary.csv", row.names = F)

