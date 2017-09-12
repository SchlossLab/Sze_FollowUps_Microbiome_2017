### Combine models from the Training runs
### Specifically identify mean and SD of model 
### and pull out most important OTUs used within it
## Marc Sze


### Note:
## Random Forest: from the R package: “For each tree, the prediction accuracy 
## on the out-of-bag portion of the data is recorded. Then the same is done after permuting 
## each predictor variable. The difference between the two accuracies are then averaged over all trees, 
## and normalized by the standard error. For regression, the MSE is computed on the out-of-bag data for each 
## tree, and then the same computed after permuting a variable. The differences are averaged and normalized 
## by the standard error. If the standard error is equal to 0 for a variable, the division is not done.”

## Cutoff of 0.5 (default was used for this) for RF model

###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "reshape2", "scales", "caret", "pROC"))

# Set up relevent environment variables
imp_vars_list <- list()
run_info_list <- list()
run_predictions <- list()
best_tune <- list()
test_data_list <- list()
pred_roc <- list()
stored_sens_sepc <- list()

n <- 100

best_model_data <- as.data.frame(matrix(nrow = 100, ncol = 8))

for(i in 1:n){
  
  load(paste("exploratory/crc_RF_model_", i, ".RData", sep=""))
  
  if(i == 1){
    write.csv(eighty_twenty_splits, "data/process/tables/crc_test_data_splits.csv", 
      row.names = F)
    write.csv(test_data, "data/process/tables/crc_test_tune_data.csv", 
      row.names = F)
  }
  
  rm(list = setdiff(ls(), 
                    c("test_tune_list", "test_predictions", "best_tune", "test_test_data",  
                      "best_model_data", "imp_vars_list", "run_info_list", "test_data_list", 
                      "run_predictions", "n", "i", "probs_predictions", "all_runs_list", "probs_test", 
                      "pred_roc", "stored_sens_sepc")))
  
  best_tune[paste("run_", i, sep = "")] <- test_tune_list[[paste(
    "data_split", i, sep = "")]]$bestTune
  
  run_info_list[[paste("run_", i, sep = "")]] <- 
    test_tune_list[[paste("data_split", i, sep = "")]]$results 
  
  ### Need to store the test_data for future use
  test_data_list[[paste("run_", i, sep="")]] <- test_test_data
  
  
  imp_vars_list[[paste("run_", i, sep = "")]] <- 
    varImp(test_tune_list[[paste("data_split", i, sep = "")]], 
           scale = FALSE)$importance %>% 
    mutate(Variable = rownames(.), run = i)
  
  
  run_predictions[[paste("run_", i, sep = "")]] <- test_predictions[[paste(
    "data_split", i, sep = "")]]
  
  
  pred_roc[[paste("run_", i, sep = "")]] <- 
    roc(test_data_list[[paste("run_", i, sep="")]]$lesion ~ 
          factor(run_predictions[[paste("run_", i, sep = "")]], ordered = TRUE))
  
  stored_sens_sepc[[paste("run_", i, sep = "")]] <- 
    cbind(sens = pred_roc[[paste("run_", i, sep = "")]]$sensitivities, 
          spec = pred_roc[[paste("run_", i, sep = "")]]$specificities)
  
  best_model_data[i, ] <- c(filter(run_info_list[[i]], 
                                   mtry == best_tune[[i]])[1, ], pred_roc[[paste("run_", i, sep = "")]]$auc)
  
  colnames(best_model_data) <- c("mtry", "ROC", "Sens", "Spec", "ROCSD", 
                                 "SensSD", "SpecSD", "test_auc")
  
  
  rownames(best_model_data)[i] <- paste("run_", i, sep = "")
  
  rm(test_tune_list, test_test_data)
  
}


# Write out ROC summary table

write.csv(
  mutate(best_model_data, run = rownames(best_model_data)), 
  "data/process/tables/crc_ROC_model_summary.csv", row.names = F)


# Get Ranges of 100 10-fold 20 times CV data (worse, best)
best_run <- as.numeric(strsplit((
  mutate(best_model_data, run = rownames(best_model_data)) %>% 
    filter(test_auc == max(test_auc)) %>% 
    select(run))[1,], "_")[[1]][2])

worse_run <- as.numeric(strsplit((mutate(best_model_data, 
                                         run = rownames(best_model_data)) %>% 
                                    filter(test_auc == min(test_auc)) %>% 
                                    select(run))[1,], "_")[[1]][2])


roc_data_list <- list(
  best_roc = pred_roc[[best_run]], 
  worse_roc = pred_roc[[worse_run]])

# Get sensitivity and specificity for test data (best and worse)
test_roc_data <- cbind(
  sensitivities = c(roc_data_list[["best_roc"]]$sensitivities, 
                    roc_data_list[["worse_roc"]]$sensitivities), 
  specificities = c(roc_data_list[["best_roc"]]$specificities, 
                    roc_data_list[["worse_roc"]]$specificities), 
  run = c(
    rep("best_roc", length(roc_data_list[["best_roc"]]$sensitivities)), 
    rep("worse_roc", length(roc_data_list[["worse_roc"]]$sensitivities))))

write.csv(test_roc_data, 
          "data/process/tables/crc_test_data_roc.csv", row.names = F)


# Aggregate all the MDA values into a single table
raw_mda_otu_data <- bind_rows(imp_vars_list)


# Generate median and IQR
summary_otu_mda <- raw_mda_otu_data %>% select(-run) %>% 
  group_by(Variable) %>% 
  summarise_all(
    funs(median_mda = median, 
         iqr25 = quantile(.)["25%"], 
         iqr75 = quantile(.)["75%"])) %>% as.data.frame() %>% 
  arrange(desc(median_mda))


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "reshape2", "scales", "caret", "pROC"))


# Read in and generate taxa labels
tax <- read.delim('data/process/final.taxonomy', sep='\t', header=T, row.names=1)

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
otu <- rownames(tax_df)

tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))

low_tax_ID <- gsub("_", " ", 
                   gsub("2", "", 
                        gsub("_unclassified", "", createTaxaLabeller(tax_df))))

test <- data_frame(otu = otu, tax_ID = low_tax_ID)

summary_otu_mda <- summary_otu_mda %>% inner_join(test, by = c("Variable" = "otu"))


# Write out the raw information for the importantance by MDA to a table
write.csv(raw_mda_otu_data, 
          "data/process/tables/crc_raw_mda_values.csv", row.names = F)


# Write out the median and IQR for the MDA of the most important variables

write.csv(summary_otu_mda, 
          "data/process/tables/crc_MDA_Summary.csv", row.names = F)












