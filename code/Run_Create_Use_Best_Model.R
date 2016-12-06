### Create actual model based on training and testing
### Use the full data set now to create model
### Workflow based on Jenna Wiens suggested workflow
## Marc Sze

#Load needed libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "caret","scales", "wesanderson", "randomForest", "pROC"))


# Read in necessary data frames

test_data <- read.csv("results/tables/full_test_data.csv", header = TRUE)
split_data_results <- read.csv("results/tables/ROC_model_summary.csv", header = TRUE, stringsAsFactors = F)
test_data_roc <- read.csv("results/tables/test_data_roc.csv", header = TRUE, stringsAsFactors = F)
auc_data_table <- read.csv("results/tables/auc_summary.csv", header = TRUE, row.names = 1)

# Get best mtry to use

maximized_mtry <- as.numeric(names(test[test == max(test)]))

# Create Random Forest model

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
