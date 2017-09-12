###Aggregate all the data from this IF model
### Get the best and worst models and correct mtry for new model
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "ggplot2", "reshape2", "gridExtra", "scales", 
           "wesanderson", "caret", "pROC"))


# Load needed data
load("exploratory/crc_treatment_model.RData")

# Set up relevent environment variables
imp_vars_list <- list()
run_info_list <- list()
run_predictions <- list()
best_tune <- list()
probs_predictions <- list()
n <- 100

best_model_data <- as.data.frame(matrix(nrow = 100, ncol = 6))


for(i in 1:n){
  
  if(i == 1){

    write.csv(test_data, "data/process/tables/crc_treatment_test_tune_data.csv", 
              row.names = F)
  }
  
  probs_predictions[[paste("run_", i, sep = "")]] <- 
    test_tune_list[[paste("run_", i, sep = "")]]$finalModel$votes %>% as.data.frame
  
  best_tune[paste("run_", i, sep = "")] <- test_tune_list[[paste(
    "run_", i, sep = "")]]$bestTune
  
  run_info_list[[paste("run_", i, sep = "")]] <- 
    test_tune_list[[paste("run_", i, sep = "")]]$results 
  
  imp_vars_list[[paste("run_", i, sep = "")]] <- 
    varImp(test_tune_list[[paste("run_", i, sep = "")]], 
           scale = FALSE)$importance %>% 
    mutate(Variable = rownames(.), run = i)
  
  
  run_predictions[[paste("run_", i, sep = "")]] <- probs_predictions[[paste(
    "data_split", i, sep = "")]]
  
  best_model_data[i, ] <- filter(run_info_list[[i]], 
                                 mtry == best_tune[[i]]) %>% select(-mtry)
  colnames(best_model_data) <- colnames(select(run_info_list[[i]], -mtry))
  rownames(best_model_data)[i] <- paste("run_", i, sep = "")
  
}


write.csv(
  mutate(best_model_data, run = rownames(best_model_data), 
         best_mtry = t(as.data.frame.list(best_tune))), 
  "data/process/tables/crc_treatment_ROC_model_summary.csv", row.names = F)

# Get Ranges of 100 10-fold 20 times CV data (worse, best)
best_run <- as.numeric(strsplit((
  mutate(best_model_data, run = rownames(best_model_data), 
         one_minus_ROC = ifelse(ROC < 0.5, invisible(1-ROC), invisible(ROC))) %>% 
    filter(one_minus_ROC == max(one_minus_ROC)) %>% 
    select(run))[1,], "_")[[1]][2])

worse_run <- as.numeric(strsplit((mutate(best_model_data, 
                                         run = rownames(best_model_data), 
                                         one_minus_ROC = ifelse(ROC < 0.5, invisible(1-ROC), invisible(ROC))) %>% 
                                    filter(one_minus_ROC == min(one_minus_ROC)) %>% 
                                    select(run))[1,], "_")[[1]][2])


roc_data_list <- list(
  best_roc = roc(
    test_data$lesion ~ 
      probs_predictions[[best_run]][, "Yes"]), 
  worse_roc = roc(test_data$lesion ~ 
                    probs_predictions[[worse_run]][, "Yes"]))

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
          "data/process/tables/crc_treatment_test_data_roc.csv", row.names = F)

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


# Write out the important variables to a table
write.csv(raw_mda_otu_data, 
          "data/process/tables/crc_treatment_raw_mda_values.csv", row.names = F)


write.csv(summary_otu_mda, 
          "data/process/tables/crc_treatment_MDA_Summary.csv", row.names = F)









