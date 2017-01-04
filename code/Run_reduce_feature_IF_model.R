### Create IF model using the important variables only
### Does this model perform similarily to the full model
## Marc Sze


#Load needed libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("randomForest", "dplyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "caret", "doMC"))

#Load needed data
test_data <- read.csv("results/tables/IF_test_tune_data.csv", header = TRUE)
lesion_imp_vars <- read.csv("results/tables/IF_rf_wCV_imp_vars_summary.csv", header = T, stringsAsFactors = F)

#Create data table with only reduced features (impvars only)
vars_to_keep <- lesion_imp_vars$Variable
test_data_imps <- select(test_data, lesion, one_of(vars_to_keep))

#################################################################################
#                                                                               #
#                                                                               #
#               Data Splitting to create train set                              #
#                                                                               #
#################################################################################

# Creat a data frame to store rows to be sampled for the training set
eighty_twenty_splits <- as.data.frame(matrix(nrow=104, ncol = 100))

# Make sure that the same individual is only in the train or test respectively
for(i in 1:length(colnames(eighty_twenty_splits))){
  
  test <- sample_n(test_data_imps[1:66, ], size = 52)
  eighty_twenty_splits[i] <- c(as.numeric(rownames(test)), 
                               as.numeric(rownames(test))+67)
  colnames(eighty_twenty_splits)[i] <- paste("split_", i, sep="")
}


#################################################################################
#                                                                               #
#                                                                               #
#               Model Training and Parameter Tuning                             #
#                                                                               #
#################################################################################



#Create Overall specifications for model tuning
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated twenty times
  repeats = 20, 
  p = 0.8, 
  classProbs = TRUE, 
  summaryFunction = twoClassSummary)


registerDoMC(cores = 8)

#Set up lists to store the data
test_tune_list <- list()
test_predictions <- list()

for(i in 1:length(colnames(eighty_twenty_splits))){
  
  #Get test data
  train_test_data <- test_data_imps[eighty_twenty_splits[, i], ]
  
  #Train the model
  set.seed(3457)
  test_tune_list[[paste("data_split", i, sep = "")]] <- 
    train(lesion ~ ., data = train_test_data, 
          method = "rf", 
          ntree = 2000, 
          trControl = fitControl, 
          metric = "ROC", 
          na.action = na.omit, 
          verbose = FALSE)
  
  
  test_test_data <- test_data[-eighty_twenty_splits[, i], ]
  
  test_predictions[[paste("data_split", i, sep = "")]] <- 
    predict(test_tune_list[[i]], test_test_data)
}


# Save image with data and relevant parameters
save.image("exploratory/IF_reduced_RF_model_Imp_OTU.RData")



