### Create model using the important variables only
### Does this model perform similarily to the full model
## Marc Sze


#Load needed libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "scales", "wesanderson", "ggplot2", "caret"))

#Load needed data
test_data <- read.csv("results/tables/full_test_data.csv", header = TRUE, row.names = 1)
lesion_imp_vars <- read.csv("results/tables/rf_wCV_imp_vars_summary.csv", header = T, stringsAsFactors = F)


#Create data table with only reduced features (impvars only)
vars_to_keep <- lesion_imp_vars$Variable
test_data_imps <- select(test_data, lesion, one_of(vars_to_keep))


#################################################################################
#                                                                               #
#                                                                               #
#               Data Splitting to create train set                              #
#                                                                               #
#################################################################################

# Split data evenly 
set.seed(3457)
eighty_twenty_splits <- createDataPartition(test_data_imps$lesion, 
                                            p = 0.8, list = FALSE, times = 100)

#Create Overall specifications for model tuning
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated twenty times
  repeats = 20, 
  p = 0.8, 
  classProbs = TRUE, 
  summaryFunction = twoClassSummary)


save.image("exploratory/RF_reduced_features_model_setup.RData")


# Create Random Forest model
#set.seed(3457)
#full_model <- randomForest(lesion ~ ., data = test_data_imps, mtry = maximized_mtry, ntree = 2000)
#full_model_probs_predictions <- full_model$votes
#full_model_roc <- roc(test_data$lesion ~ full_model_probs_predictions[, "Yes"])











