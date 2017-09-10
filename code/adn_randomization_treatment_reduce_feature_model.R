### Create IF model using the important variables only
### Does this model perform similarily to the full model
## Marc Sze


#Load needed libraries
source('code/functions.R')

loadLibs(c("randomForest", "dplyr", "scales", "caret"))

#Load needed data
test_data <- read.csv("data/process/tables/adn_randomization_treatment_test_tune_data.csv", header = TRUE)
lesion_imp_vars <- read.csv("data/process/tables/adn_randomization_treatment_imp_vars_summary.csv", header = T, stringsAsFactors = F) %>% 
  rename(variable = Variable)
mda_summary <- read.csv("data/process/tables/adn_randomization_treatment_top_vars_MDA_Summary.csv", header = T, stringsAsFactors = F)


#Merge MDA and rank counts together to get top 10.
combined_ranks <- inner_join(lesion_imp_vars, mda_summary, by = "variable") %>% 
  arrange(desc(total_appearance), desc(mean_MDA))

#Create data table with only reduced features (impvars only)
vars_to_keep <- slice(combined_ranks, 1:10)[, "variable"]
test_data_imps <- test_data %>% select(lesion, one_of(vars_to_keep$variable))


#################################################################################
#                                                                               #
#                                                                               #
#               Model Training and Parameter Tuning                             #
#                                                                               #
#################################################################################



#Create Overall specifications for model tuning
fitControl <- trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  ## repeated twenty times
  p = 0.8, 
  classProbs = TRUE, 
  savePredictions = TRUE, 
  summaryFunction = twoClassSummary)

set.seed(3457)

#Set up lists to store the data
test_tune_list <- list()
test_predictions <- list()

for(i in 1:100){
  
  #Get test data
  train_test_data <- test_data_imps
  
  #Train the model
  test_tune_list[[paste("data_split", i, sep = "")]] <- 
    train(lesion ~ ., data = train_test_data, 
          method = "rf", 
          ntree = 100, 
          trControl = fitControl, 
          metric = "ROC", 
          na.action = na.omit, 
          verbose = FALSE)
  
  
  test_predictions[[paste("data_split", i, sep = "")]] <- 
    predict(test_tune_list[[i]], test_data_imps)
}


# Save image with data and relevant parameters
save.image("exploratory/adn_randomization_treatment_reduced_RF_model_Imp_OTU.RData")



