### Create model using the important variables only
### Does this model perform similarily to the full model
## Marc Sze


#Load needed libraries
source('code/functions.R')

loadLibs(c("dplyr", "scales", "caret"))

#Load needed data
test_data <- read.csv("results/tables/full_test_data.csv", header = TRUE, row.names = 1)
common_imp_vars <- read.csv("results/tables/pvalue_IF_lesion_common_imp_vars.csv", 
                            header = T, stringsAsFactors = F, row.names = 1)


#Create data table with only reduced features (impvars only)
vars_to_keep <- common_imp_vars$otu
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


save.image("exploratory/RF_common_features_model_setup.RData")


