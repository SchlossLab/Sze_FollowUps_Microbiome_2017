### Create actual model based on training and testing
### Use the full data set now to create model
### Workflow based on Jenna Wiens suggested workflow
## Marc Sze

#Load needed libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "caret","scales", "wesanderson", "randomForest", "pROC"))


# Read in necessary data frames
load("exploratory/rocs.RData")
test_data <- read.csv("results/tables/full_test_data.csv", header = TRUE, row.names = 1)
split_data_results <- read.csv("results/tables/ROC_model_summary.csv", header = TRUE, stringsAsFactors = F)
test_data_roc <- read.csv("results/tables/test_data_roc.csv", header = TRUE, stringsAsFactors = F)
auc_data_table <- read.csv("results/tables/auc_summary.csv", header = TRUE, row.names = 1)

# Get best mtry to use
mtry_table <- table(split_data_results$best_mtry)
maximized_mtry <- as.numeric(names(mtry_table[mtry_table == max(mtry_table)]))

# Create Random Forest model

set.seed(3457)
full_model <- randomForest(lesion ~ ., data = test_data, mtry = maximized_mtry, ntree = 2000)

full_model_probs_predictions <- full_model$votes

full_model_roc <- roc(test_data$lesion ~ full_model_probs_predictions[, "Yes"])

# Add to needed data tables and write them back out

auc_data_table <- rbind(auc_data_table, full = c(full_model_roc$auc, NA, NA, NA))
write.csv(auc_data_table, "results/tables/auc_summary.csv")

test_data_roc <- rbind(test_data_roc, 
                       cbind(sensitivities = full_model_roc$sensitivities, 
                             specificities = full_model_roc$specificities, 
                             run = rep("full_roc", length(full_model_roc$sensitivities))))

write.csv(test_data_roc, 
          "results/tables/test_data_roc.csv", row.names = F)

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


# Load in test data set
follow_up_data <- read.csv("results/tables/follow_up_prediction_table.csv")

test_follow_up_data <- cbind(select(follow_up_data, -contains("Otu0")), 
                   select(shared, -Group))


# Make predictions on samples with follow up
initial_predictions <- predict(full_model, 
  newdata = test_follow_up_data[1:(length(rownames(test_follow_up_data))/2), ], 
  type='prob')

followup_predictions <- predict(full_model, 
  newdata = test_follow_up_data[((length(rownames(
    follow_up_data))/2)+1):length(rownames(test_follow_up_data)), ], 
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
          "results/tables/follow_up_probability_summary.csv", row.names = F)

# Create data table for significance ROC testing














