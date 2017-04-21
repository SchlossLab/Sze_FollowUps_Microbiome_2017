### CRC treatment model
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("randomForest", "dplyr", "ggplot2", "reshape2", 
  "gridExtra", "scales", "wesanderson", "caret", "doMC"))

# Read in data tables
good_metaf <- read.csv(
  "data/process/mod_metadata/good_metaf_final.csv", 
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
good_metaf_select <- data_frame(sample = c(good_metaf$initial, good_metaf$followUp)) %>% 
  mutate(lesion = c(rep(0, length(good_metaf$EDRN)), rep(1, length(good_metaf$EDRN))), 
         Dx_Bin = rep(good_metaf$Dx_Bin, 2))

test_data <- inner_join(good_metaf_select, shared, 
  by = c("sample" = "Group"))

rm(shared)

# Filter and use only specific data

test_data <- test_data %>% 
  filter(Dx_Bin == "cancer") %>% 
  select(sample, lesion, contains("Otu"))

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
test_data$lesion <- factor(test_data$lesion, 
  levels = c(0, 1), labels = c("No", "Yes"))

#write out table for future use

write.csv(test_data, "data/process/tables/crc_treatment_model.csv", 
  row.names = F)


#################################################################################
#                                                                               #
#                                                                               #
#               Model Training and Parameter Tuning                             #
#                                                                               #
#################################################################################

#Create Overall specifications for model tuning
# number controls fold of cross validation
# Repeats control the number of times to run it

fitControl <- trainControl(## 5-fold CV
  method = "cv",
  number = 10,
  p = 0.8, 
  classProbs = TRUE, 
  summaryFunction = twoClassSummary, 
  savePredictions = "final")

registerDoMC(cores = 8)

#Set up lists to store the data
test_tune_list <- list()
test_predictions <- list()

for(i in 1:100){
  
  #Train the model
  set.seed(3457)
  test_tune_list[[paste("run_", i, sep = "")]] <- 
    train(lesion ~ ., data = test_data, 
          method = "rf", 
          ntree = 2000, 
          trControl = fitControl, 
          metric = "ROC", 
          na.action = na.omit, 
          verbose = FALSE)
  
  
  test_predictions[[paste("data_split", i, sep = "")]] <- 
  predict(test_tune_list[[i]], test_data)
}


# Save image with data and relevant parameters
save.image("exploratory/crc_treatment_model.RData")





