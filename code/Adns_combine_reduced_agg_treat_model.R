###Aggregate all the data from this IF model
### Get the best and worst models and correct mtry for new model
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "scales", "caret", "pROC", "tidyr"))


# Load needed data
load("exploratory/adn_treatment_reduced_RF_model_Imp_OTU.RData")

#Read in needed data
tax <- read.delim('data/process/final.taxonomy', sep='\t', header=T, row.names=1, stringsAsFactors = F)
# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)


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

    write.csv(test_data, "data/process/tables/red_adn_treatment_test_tune_data.csv", 
              row.names = F)
  }
  
  probs_predictions[[paste("run_", i, sep = "")]] <- 
    test_tune_list[[paste("data_split", i, sep = "")]]$finalModel$votes %>% as.data.frame
  
  best_tune[paste("run_", i, sep = "")] <- test_tune_list[[paste(
    "data_split", i, sep = "")]]$bestTune
  
  run_info_list[[paste("run_", i, sep = "")]] <- 
    test_tune_list[[paste("data_split", i, sep = "")]]$results 
  
  imp_vars_list[[paste("run_", i, sep = "")]] <- 
    varImp(test_tune_list[[paste("data_split", i, sep = "")]], 
           scale = FALSE)$importance %>% 
    mutate(Variable = rownames(.)) %>% arrange(desc(Overall))
  
  
  run_predictions[[paste("run_", i, sep = "")]] <- test_predictions[[paste(
    "data_split", i, sep = "")]]
  
  best_model_data[i, ] <- filter(run_info_list[[i]], 
                                 mtry == best_tune[[i]]) %>% select(-mtry)
  colnames(best_model_data) <- colnames(select(run_info_list[[i]], -mtry))
  rownames(best_model_data)[i] <- paste("run_", i, sep = "")
  
}


write.csv(
  mutate(best_model_data, run = rownames(best_model_data), 
         best_mtry = t(as.data.frame.list(best_tune))), 
  "data/process/tables/reduced_adn_treatment_ROC_model_summary.csv", row.names = F)

# Get Ranges of 100 10-fold 20 times CV data (worse, best)
best_run <- as.numeric(strsplit((
  mutate(best_model_data, run = rownames(best_model_data)) %>% 
    filter(ROC == max(best_model_data$ROC)) %>% 
    select(run))[1,], "_")[[1]][2])

worse_run <- as.numeric(strsplit((mutate(best_model_data, 
                                         run = rownames(best_model_data)) %>% 
                                    filter(ROC == min(best_model_data$ROC)) %>% 
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
          "data/process/tables/reduced_adn_treatment_test_data_roc.csv", row.names = F)


# Collect the mean and SD for the MDA of the most important variables

top_vars_MDA <- lapply(imp_vars_list, function(x) 
  x[order(x[, "Variable"]), ])

top_vars_MDA_by_run <- as.data.frame(matrix(nrow = length(top_vars_MDA$run_1$Variable), 
                                            ncol = length(imp_vars_list), 
                                            dimnames = list(
                                              nrow = top_vars_MDA[["run_1"]]$Variable, 
                                              ncol = paste("run_", seq(1:100), sep = ""))))

for(i in 1:length(top_vars_MDA_by_run)){
  
  tempData <- top_vars_MDA[[i]]
  rownames(tempData) <- tempData$Variable
  top_vars_MDA_by_run[, i] <- tempData[rownames(top_vars_MDA_by_run), "Overall"]
  rm(tempData)
}


# "1" pulls the value of mean or sd from the data frame
MDA_vars_summary <- as.data.frame(t(top_vars_MDA_by_run)) %>% summarise_all(funs(medianMDA = median, medianIQR = IQR)) %>% 
  gather(key = tempName) %>% separate(tempName, c("otu", "measurement")) %>% 
  spread(key = measurement, value = value) %>% rename(median_MDA = medianMDA, median_IQR = medianIQR) %>% 
  mutate(tax_id = as.character(tax_df[otu, "Genus"])) %>% 
  arrange(desc(median_MDA)) %>% 
  mutate(rank = seq(1:length(rownames(top_vars_MDA_by_run))))


write.csv(MDA_vars_summary, 
          "data/process/tables/reduced_adn_treatment_top_vars_MDA_Summary.csv", row.names = F)

adn_model_top_vars_MDA_full_data <- mutate(top_vars_MDA_by_run, otu = rownames(top_vars_MDA_by_run)) %>% 
  slice(match(MDA_vars_summary$otu, otu)) %>% 
  gather(key = variable, value = value, -otu)
  

write.csv(adn_model_top_vars_MDA_full_data, 
          "data/process/tables/reduced_adn_treatment_top_vars_MDA_full_data.csv", row.names = F)







