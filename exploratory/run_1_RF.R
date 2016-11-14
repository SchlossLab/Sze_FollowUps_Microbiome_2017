### Build the best lesion model possible
### Try XG-Boost, RF, Logit (GLM), C5.0, SVM
### Find the best based on Jenna Wiens suggestions on test and training
## Marc Sze

#Load needed libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "caret","scales", "wesanderson", "doMC"))

load("exploratory/RF_model_setup.RData")

# Set i variable

i = 1

#################################################################################
#                                                                               #
#                                                                               #
#               Model Training and Parameter Tuning                             #
#                                                                               #
#################################################################################

# Create the tune grid to be used
tunegrid <- expand.grid(
  .mtry=c(1:441), .ntree=c(250, 500, 1000, 1500, 2000))

# Call number of processors to use
registerDoMC(cores = 10)

#Set up lists to store the data
test_tune_list <- list()
test_predictions <- list()

#Get test data
train_test_data <- test_data[eighty_twenty_splits[, i], ]
  
#Train the model
train_name <- paste("data_split", i, sep = "")

set.seed(3457)
assign(train_name, 
  train(lesion ~ ., data = train_test_data, 
        method = customRF, 
        tuneGrid = tunegrid, 
        trControl = fitControl, 
        metric = "ROC", 
        verbose = FALSE))
  
  
test_test_data <- test_data[-eighty_twenty_splits[, i], ]
  
test_predictions[[paste("data_split", i, sep = "")]] <- 
predict(test_tune_list[[paste("data_split", i, sep = "")]], 
  test_test_data)


# Save image with data and relevant parameters
save.image(paste("exploratory/RF_model_", i, ".RData", sep=""))
