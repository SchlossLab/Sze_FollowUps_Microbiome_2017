### Build the best lesion model possible
### Try XG-Boost, RF, Logit (GLM), C5.0, SVM
### Find the best based on Jenna Wiens suggestions on test and training
## Marc Sze

#Load needed libraries
source('code/functions.R')

loadLibs(c("dplyr", "caret","scales", "doMC"))

load("exploratory/srn_RF_model_setup.RData")

# Set i variable

i = 99

#################################################################################
#                                                                               #
#                                                                               #
#               Model Training and Parameter Tuning                             #
#                                                                               #
#################################################################################

# Call number of processors to use
registerDoMC(cores = 4)

#Set up lists to store the data
test_tune_list <- list()
test_predictions <- list()

#Get test data
train_test_data <- test_data[eighty_twenty_splits[, i], ]
  
#Train the model
train_name <- paste("data_split", i, sep = "")

set.seed(3457)
test_tune_list[[paste("data_split", i, sep = "")]] <- assign(train_name, 
  train(lesion ~ ., data = train_test_data, 
        method = "rf", 
        ntree = 2000, 
        trControl = fitControl, 
        metric = "ROC", 
        verbose = FALSE))
   
test_test_data <- test_data[-eighty_twenty_splits[, i], ]
  
test_predictions[[paste("data_split", i, sep = "")]] <- 
predict(test_tune_list[[paste("data_split", i, sep = "")]], 
  test_test_data)


# Save image with data and relevant parameters
save.image(paste("exploratory/srn_RF_model_", i, ".RData", sep=""))
