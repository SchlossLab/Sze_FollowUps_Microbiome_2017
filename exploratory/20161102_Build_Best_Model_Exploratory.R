### Build the best lesion model possible
### Try XG-Boost, RF, Logit (GLM), C5.0, SVM
### Find the best based on Jenna Wiens suggestions on test and training
## Marc Sze

#Load needed libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "caret","scales", "wesanderson"))


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

# Center, Scale, Yeo Johnson Transform, and pca for better modeling
woPCA <- preProcess(test_data, 
  method = c("center", "scale", "YeoJohnson"))  

wPCA <- preProcess(test_data, 
  method = c("center", "scale", "YeoJohnson", "pca"), thresh = 0.99)  

woPCA_test_data <- predict(woPCA, test_data)

wPCA_test_data <- predict(wPCA, test_data)


#################################################################################
#                                                                               #
#                                                                               #
#               Data Splitting to create train set                              #
#                                                                               #
#################################################################################

# Split data evenly 
set.seed(3457)
eighty_twenty_splits <- createDataPartition(woPCA$lesion, 
  p = 0.8, list = FALSE, times = 100)


#################################################################################
#                                                                               #
#                                                                               #
#               Model Training and Parameter Tuning                             #
#                                                                               #
#################################################################################

#Get test data
train_test_data <- woPCA[eighty_twenty_splits[, 1], ]

#Create Overall specifications for model tuning
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10, 
  p = 0.8, 
  classProbs = TRUE, 
  summaryFunction = twoClassSummary)

#Train the model
set.seed(3457)
test_tune <- train(lesion ~ ., data = train_test_data, 
                 method = "rf", 
                 trControl = fitControl,
                 metric = "Sens", 
                 verbose = FALSE)



save.image("exploratory/test_model_RF.RData")
