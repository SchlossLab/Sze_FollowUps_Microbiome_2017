### Build the best lesion model possible
### Find the best based on Jenna Wiens suggestions on test and training
## Marc Sze

#Load needed libraries
source('code/functions.R')

loadLibs(c("dplyr", "caret","scales", "doMC"))

# Read in necessary data frames

shared <- read.delim('data/process/final.0.03.subsample.shared', 
  header=T, sep='\t') %>% select(Group, contains("Otu0"))

metaI <- read.csv("data/process/mod_metadata/metaI_final.csv", 
  stringsAsFactors = F, header = T) 

good_metaf <- read.csv('data/process/mod_metadata/good_metaf_final.csv', 
                       header = T, stringsAsFactors = F) %>% select(initial)

metaI <- filter(metaI, !(sample %in% good_metaf$initial))

#################################################################################
#                                                                               #
#                                                                               #
#               Data Clean up for better and faster modeling                    #
#                                                                               #
#################################################################################

# Remove follow up samples and join metadata with microbiome data
test_data <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(dx == "cancer" | dx == "normal")
rm(shared)

# Filter and use only specific data

test_data <- test_data %>% 
  select(sample, contains("Otu0"))

#Filter out rows that are not all complete

test_data <- test_data[complete.cases(test_data), ]

#Reduce MetaI to match the test_data
metaI <- filter(metaI, dx == "cancer" | dx == "normal")

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

# Add the lesion variable to the test_data
test_data <- cbind(lesion = factor(metaI$lesion, 
  levels = c(0, 1), labels = c("No", "Yes")), test_data)

# Write data table for future use
write.csv(test_data, "data/process/tables/crc_full_test_data.csv")


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
save.image("exploratory/crc_RF_model_setup.RData")
