### Create actual model based on training and testing
### Use the full data set now to create model
### Workflow based on Jenna Wiens suggested workflow
## Marc Sze

#Load needed libraries
source('code/functions.R')

loadLibs(c("dplyr", "caret","scales", "wesanderson", "randomForest", "pROC"))


# Read in necessary data frames
load("exploratory/rocs.RData")
test_data <- read.csv("data/process/tables/full_test_data.csv", header = TRUE, row.names = 1)
split_data_results <- read.csv("data/process/tables/ROC_model_summary.csv", header = TRUE, stringsAsFactors = F)
test_data_roc <- read.csv("data/process/tables/test_data_roc.csv", header = TRUE, stringsAsFactors = F)
auc_data_table <- read.csv("data/process/tables/auc_summary.csv", header = TRUE, row.names = 1)

# Get best mtry to use
mtry_table <- table(split_data_results$best_mtry)
maximized_mtry <- as.numeric(names(mtry_table[mtry_table == max(mtry_table)]))

# Create Random Forest model

set.seed(3457)
full_model <- randomForest(lesion ~ ., data = test_data, mtry = maximized_mtry, ntree = 2000)

full_model_probs_predictions <- full_model$votes

full_model_roc <- roc(test_data$lesion ~ full_model_probs_predictions[, "Yes"])

# Add to needed data tables and write them back out

auc_data_table <- rbind(auc_data_table, full = c(full_model_roc$auc, NA, NA, NA))
write.csv(auc_data_table, "data/process/tables/auc_summary.csv")
write.csv(cbind(lesion = as.character(test_data$lesion), full_model_probs_predictions), 
          "data/process/tables/lesion_test_data_probs_summary.csv")

test_data_roc <- rbind(test_data_roc, 
                       cbind(sensitivities = full_model_roc$sensitivities, 
                             specificities = full_model_roc$specificities, 
                             run = rep("full_roc", length(full_model_roc$sensitivities))))

write.csv(test_data_roc, 
          "data/process/tables/test_data_roc.csv", row.names = F)

#load in metadata to get IDs to select shared file by
good_metaf <- read.csv("data/process/mod_metadata/good_metaf_final.csv", 
                       stringsAsFactors = F, header = T)

#create vector in correct order

samples <- c(good_metaf$initial, good_metaf$followUp)

#read in shared file and keep only the samples that are needed
shared <- read.delim('data/process/final.0.03.subsample.shared', 
                     header=T, sep='\t') %>% select(-label, -numOtus)
rownames(shared) <- shared$Group
shared <- shared[as.character(samples), ]

# Keep only OTUs in test data
OTUs_to_keep <- colnames(select(test_data, -lesion))
shared <- select(shared, Group, one_of(OTUs_to_keep))


# Load in test data set
follow_up_data <- read.csv("data/process/tables/follow_up_prediction_table.csv")

test_follow_up_data <- cbind(select(follow_up_data, -contains("Otu0")), 
                   select(shared, -Group))


# Make predictions on samples with follow up
initial_predictions <- predict(full_model, 
  newdata = test_follow_up_data[1:(length(rownames(test_follow_up_data))/2), ], 
  type='prob')

followup_predictions <- predict(full_model, 
  newdata = test_follow_up_data[((length(rownames(
    follow_up_data))/2)+1):length(rownames(test_follow_up_data)), ], 
  type='prob')


# Create data table needed for figure 4

probability_data_table <- cbind(
  No = c(initial_predictions[, "No"], followup_predictions[, "No"]), 
  Yes = c(initial_predictions[, "Yes"], followup_predictions[, "Yes"]), 
  sampleType = c(rep("initial", length(rownames(initial_predictions))), 
                 rep("followup", length(rownames(followup_predictions)))), 
  disease_free = rep(good_metaf$Disease_Free, 2), 
  diagnosis = rep(good_metaf$Diagnosis, 2), 
  Dx_Bin = rep(good_metaf$Dx_Bin, 2), 
  chemo = rep(good_metaf$chemo_received, 2), 
  rads = rep(good_metaf$radiation_received, 2), 
  EDRN = rep(good_metaf$EDRN, 2))

write.csv(probability_data_table, 
          "data/process/tables/follow_up_probability_summary.csv", row.names = F)

# Create data table for significance ROC testing

roc_data_list[["full_roc"]] <- full_model_roc

p_values_summary <- matrix(nrow = 4, ncol = 3, dimnames = list(
  nrow = c("best", "middle", "worse", "full"), 
  ncol = c("middle", "worse", "full")))

BH_pvalues_summary <- matrix(nrow = 4, ncol = 3, dimnames = list(
  nrow = c("best", "middle", "worse", "full"), 
  ncol = c("middle", "worse", "full")))

for(i in 1:length(roc_data_list)){
  
  for(j in i:length(roc_data_list)){
    
    if(i != j){
      
      p_values_summary[i, j-1] <- roc.test(roc_data_list[[i]], roc_data_list[[j]])$p.value
      BH_pvalues_summary[i, j-1] <- p.adjust(p_values_summary[i, j-1], method = "BH", n = 6)
      
    }
  }
}

write.csv(BH_pvalues_summary, "data/process/tables/roc_pvalue_summary.csv")
write.csv(p_values_summary, "data/process/tables/roc_nonBH_pvalue_summary.csv")







