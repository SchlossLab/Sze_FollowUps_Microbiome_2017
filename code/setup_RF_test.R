### Build the best lesion model possible
### Try XG-Boost, RF, Logit (GLM), C5.0, SVM
### Find the best based on Jenna Wiens suggestions on test and training
## Marc Sze

#Load needed libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "caret","scales", "wesanderson", "doMC"))


# Read in necessary data frames

shared <- read.delim('data/process/final.0.03.subsample.shared', 
  header=T, sep='\t') %>% select(Group, contains("Otu0"))

metaI <- read.csv("results/tables/mod_metadata/metaI_final.csv", 
  stringsAsFactors = F, header = T) 

#################################################################################
#                                                                               #
#                                                                               #
#               Data Clean up for better and faster modeling                    #
#                                                                               #
#################################################################################

# Remove follow up samples and join metadata with microbiome data
test_data <- inner_join(metaI, shared, by = c("sample" = "Group"))
rm(shared)

# Filter and use only specific data

test_data <- test_data %>% 
  select(sample, fit_result, Hx_Prev, Hx_Fam_CRC, White, 
    BMI, Age, Gender, contains("Otu0")) %>% 
  mutate(Gender = factor(Gender)) %>% 
  mutate(Hx_Prev = factor(Hx_Prev), 
    Hx_Fam_CRC = factor(Hx_Fam_CRC), White = factor(White))

#Filter out rows that are not all complete

test_data <- test_data[complete.cases(test_data), ]

#Reduce MetaI to match the test_data
metaI <- filter(
  metaI, as.character(sample) %in% as.character(test_data$sample))

# Need to create dummy variables for those that are factor classes
temp_dummy <- dummyVars(~Gender + Hx_Prev + Hx_Fam_CRC + White, 
  data = test_data)
useable_dummy_vars <- as.data.frame(predict(temp_dummy, 
  newdata = test_data))

# Combine data frame back together and covert samples to rownames
test_data <- useable_dummy_vars %>% 
mutate(sample = (test_data$sample)) %>% 
  inner_join(test_data, by = "sample") %>% 
  select(-Gender, -Hx_Prev, -Hx_Fam_CRC, -White)

rownames(test_data) <- test_data$sample
test_data <- select(test_data, -sample)

# Remove those not in at least 1% of all samples
#microbiome_data <- microbiome_data[, colSums(microbiome_data > 1)/length(rownames(microbiome_data)) > 0.01]

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

test_matrix <- cor(test_data, method = "spearman")

summary(test_matrix[upper.tri(test_matrix)])

highly_correlated <- findCorrelation(test_matrix, cutoff = 0.75)

test_data <- test_data[, -highly_correlated]

rm(test_matrix, temp_dummy, nzv, highly_correlated)

# Look for linear dependencies
    # Looks for sets of linear combinations and removes columns
    # based on whether or not the linear combination still exits 
    # with or without the variable present

test_for_combos <- findLinearCombos(test_data) # Found none so don't need to worry about this

# Creates factor variables
test_data <- test_data %>% 
  mutate(
    Gender.m = factor(Gender.m, 
      levels = c(0, 1), labels = c("No", "Yes")), 
    Hx_Prev.1 = factor(Hx_Prev.1, 
      levels = c(0, 1), labels = c("No", "Yes")), 
    Hx_Fam_CRC.1 = factor(Hx_Fam_CRC.1, 
      levels = c(0, 1), labels = c("No", "Yes")), 
    White.1 = factor(White.1, 
      levels = c(0, 1), labels = c("No", "Yes"))) 

# Add the lesion variable to the test_data
test_data <- cbind(lesion = factor(metaI$lesion, 
  levels = c(0, 1), labels = c("No", "Yes")), test_data)

# Write data table for future use
write.csv(test_data, "results/tables/full_test_data.csv")


#################################################################################
#                                                                               #
#                                                                               #
#               Data Splitting to create train set                              #
#                                                                               #
#################################################################################

# Split data evenly 
set.seed(3457)
eighty_twenty_splits <- createDataPartition(test_data$lesion, 
  p = 0.8, list = FALSE, times = 100)


#################################################################################
#                                                                               #
#                                                                               #
#               Model Training and Parameter Tuning                             #
#                                                                               #
#################################################################################

#Get test data
train_test_data <- test_data[eighty_twenty_splits[, 1], ]

### Create a custom RF list to train and tune both mtry and ntree

#Create the initial list with baseline parameters
customRF <- list(type = "Classification", 
  library = "randomForest", loop = NULL)
# Create parameters for algorithm to be used
customRF$parameters <- data.frame(
  parameter = c("mtry", "ntree"), 
  class = rep("numeric", 2), 
  label = c("mtry", "ntree"))
# Function to look for the tunegrid data frame
customRF$grid <- function(
  x, y, len = NULL, search = "grid") {}
# Function used to execute the random forest algorithm
customRF$fit <- function(
  x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
# Function used to generate the prediction on the seperated data
customRF$predict <- function(
  modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
# Function used to generate the probabilities from the seperated data
customRF$prob <- function(
  modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
# Used to order by mtry
customRF$sort <- function(x) x[order(x[,1]),]
# Used to pull out variables optimized
customRF$levels <- function(x) x$classes


# Create the tune grid to be used
tunegrid <- expand.grid(
  .mtry=c(1:441), .ntree=c(250, 500, 1000, 1500, 2000))

#Create Overall specifications for model tuning
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated twenty times
  repeats = 20, 
  p = 0.8, 
  classProbs = TRUE, 
  summaryFunction = twoClassSummary)

# Save image with data and relevant parameters
save.image("exploratory/RF_model_setup.RData")
