### Try to identify the important OTUs 
### Specifically try and identify which OTUs are involved
### with the decrease in positive probability between intial and follow up
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("randomForest", "dplyr", "ggplot2", "reshape2", 
  "gridExtra", "scales", "wesanderson", "caret", "doMC"))

# Read in data tables
good_metaf <- read.csv(
  "results/tables/mod_metadata/good_metaf_final.csv", 
  stringsAsFactors = F, header = T) %>% 
  mutate(lesion_follow = ifelse(Disease_Free == "n", 1, 0))

shared <- read.delim('data/process/final.0.03.subsample.shared', 
  header=T, sep='\t') %>% select(-label, -numOtus)


#################################################################################
#                                                                               #
#                                                                               #
#               Data Clean up for better and faster modeling                    #
#                                                                               #
#################################################################################

#Create an insync version to merge with filtered_shared
good_metaf_select <- data_frame(
  sample = c(good_metaf$initial, good_metaf$followUp)) %>% 
  mutate(lesion = c(good_metaf$lesion, good_metaf$lesion_follow), 
         fit_result = c(good_metaf$fit_result, good_metaf$fit_followUp))

test_data <- inner_join(good_metaf_select, shared, 
  by = c("sample" = "Group"))

rm(shared)

# Filter and use only specific data

test_data <- test_data %>% 
  select(sample, lesion, fit_result, contains("Otu0")) %>% 
  mutate(lesion = factor(lesion))

#Filter out rows that are not all complete

test_data <- test_data[complete.cases(test_data), ]

# Need to create dummy variables for those that are factor classes
#temp_dummy <- dummyVars(~lesion,  data = test_data)
#useable_dummy_vars <- as.data.frame(predict(temp_dummy, 
#  newdata = test_data))

# Combine data frame back together and covert samples to rownames
#test_data <- useable_dummy_vars %>% 
#mutate(sample = (test_data$sample)) %>% 
#  inner_join(test_data, by = "sample") %>% 
#  select(-lesion)

rownames(test_data) <- test_data$sample
test_data <- select(test_data, -sample)

# Remove those with near zero variance 
# Similar to removal of singletons or removal based on percentage in sample
# Using near zero variance should cover the above two so could use only this instead

nzv <- nearZeroVar(test_data)
# By default, nearZeroVar will return the positions of the variables that are flagged to be problematic
# using the saveMetrics argument will output a table with more in-depth info
# Group stays in this case (because it is numberic) 
# but maybe better in future to transfer this to row names if not a numberic label

test_data <- test_data[, -nzv]

# Find samples that are very highly correlated with each other and only keep one of the two
# Cutoff in manual is above r value of 0.75
# Do not have any OTUs that are correlated at 0.999
# Have some perfectly negatively correlated so will have to remove those
# similar to nzv provides vector of columns that need to be removed

#test_matrix <- cor(test_data, method = "spearman")

#summary(test_matrix[upper.tri(test_matrix)])

#highly_correlated <- findCorrelation(test_matrix, cutoff = 0.75)

#test_data <- test_data[, -highly_correlated]

#rm(test_matrix, temp_dummy, nzv, highly_correlated)

# Look for linear dependencies
# Looks for sets of linear combinations and removes columns
# based on whether or not the linear combination still exits 
# with or without the variable present

#test_for_combos <- findLinearCombos(test_data) # Found some 
# This is different from the initial samples so could be due to the limited samples
# Do not remove and keep in this instance


# Creates factor variables
#test_data <- test_data %>% 
#  mutate(Gender.m = factor(Gender.m, 
#    levels = c(0, 1), 
#    labels = c("No", "Yes")), 
#  Hx_Prev.1 = factor(Hx_Prev.1, 
#    levels = c(0, 1), 
#    labels = c("No", "Yes")), 
#  Hx_Fam_CRC.1 = factor(Hx_Fam_CRC.1, 
#    levels = c(0, 1), 
#    labels = c("No", "Yes")), 
#  White.1 = factor(White.1, 
#    levels = c(0, 1), 
#    labels = c("No", "Yes"))) 

# Add the lesion variable to the test_data
test_data <- test_data[-67, ]

test_data$lesion <- factor(test_data$lesion, 
  levels = c(0, 1), labels = c("No", "Yes"))

#write out table for future use

write.csv(test_data, "results/tables/follow_up_prediction_table.csv", 
  row.names = F)

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
  
  test <- sample_n(test_data[1:66, ], size = 52)
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

### Create a custom RF list to train and tune both mtry and ntree

#Create the initial list with baseline parameters
customRF <- list(type = "Classification", library = "randomForest", 
  loop = NULL)
# Create parameters for algorithm to be used
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), 
                                  class = rep("numeric", 2), 
                                  label = c("mtry", "ntree"))
# Function to look for the tunegrid data frame
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
# Function used to execute the random forest algorithm
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
# Function used to generate the prediction on the seperated data
customRF$predict <- function(modelFit, newdata, preProc = NULL, 
  submodels = NULL)
  predict(modelFit, newdata)
# Function used to generate the probabilities from the seperated data
customRF$prob <- function(modelFit, newdata, preProc = NULL, 
  submodels = NULL)
  predict(modelFit, newdata, type = "prob")
# Used to order by mtry
customRF$sort <- function(x) x[order(x[,1]),]
# Used to pull out variables optimized
customRF$levels <- function(x) x$classes


# Create the tune grid to be used
tunegrid <- expand.grid(.mtry=c(1:441), 
  .ntree=c(250, 500, 1000, 1500, 2000))


#Create Overall specifications for model tuning
# number controls fold of cross validation
# Repeats control the number of times to run it

fitControl <- trainControl(## 5-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
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
  train_test_data <- test_data[eighty_twenty_splits[, i], ]
  
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
save.image("exploratory/RF_model_Imp_OTU.RData")





