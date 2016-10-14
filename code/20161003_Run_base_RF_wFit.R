### Run the first Random Forest with Fit
### Base analysis without metadata or OTU groupings
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))


# Read in data tables
metaI <- read.csv("results/tables/mod_metadata/metaI_final.csv", stringsAsFactors = F, header = T)
metaF <- read.csv("results/tables/mod_metadata/metaF_final.csv", stringsAsFactors = F, header = T)
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)
shared <- read.delim('data/process/followUps.final.an.unique_list.0.03.subsample.0.03.filter.shared', 
                     header=T, sep='\t')

# createList with all data tables stored as a list

train <- list(cancer = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(cancer, fit_result, contains("Otu0")) %>% na.omit() %>% mutate(cancer = factor(cancer)), 
              lesion = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(lesion, fit_result, contains("Otu0")) %>% na.omit() %>% mutate(lesion = factor(lesion)), 
              SRNlesion = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(SRNlesion, fit_result, contains("Otu0")) %>% na.omit() %>% mutate(SRNlesion = factor(SRNlesion)), 
              threeway = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(threeway, fit_result, contains("Otu0")) %>% na.omit() %>% mutate(threeway = factor(threeway)))


# Set up initial list to store the run data
orig_RF_run <- list(cancer = c(), lesion = c(), SRNlesion = c(), threeway = c())
orig_rf_opt <- list(cancer = c(), lesion = c(), SRNlesion = c(), threeway = c())
orig_probs <- list(cancer = c(), lesion = c(), SRNlesion = c(), threeway = c())
orig_roc <- list(cancer = c(), lesion = c(), SRNlesion = c(), threeway = c())

# Run the actual AUCRF for each different model
for(i in 1:length(train)){
  
  set.seed(050416)
  orig_RF_run[[i]] <- AUCRF(as.formula(paste(colnames(train[[i]])[1], " ~ .", sep = "")), 
                           data=train[[i]], pdel=0.05, ntree=500, ranking='MDA')
  orig_rf_opt[[i]] <- orig_RF_run[[i]]$RFopt
  orig_probs[[i]] <- predict(orig_rf_opt[[i]], type = 'prob')[, 2]
  orig_roc[[i]] <- roc(train[[i]][, 1] ~ orig_probs[[i]])

}  


#Set up Boruta picked data list

impfactor_Data_List <- list(cancer = c(), lesion = c(), SRNlesion = c(), threeway = c())
confirmed_vars <- list(cancer = c(), lesion = c(), SRNlesion = c(), threeway = c())

# Set up the lists to store the selected run data
selected_RF_run <- list(cancer = c(), lesion = c(), SRNlesion = c(), threeway = c())
selected_rf_opt <- list(cancer = c(), lesion = c(), SRNlesion = c(), threeway = c())
selected_probs <- list(cancer = c(), lesion = c(), SRNlesion = c(), threeway = c())
selected_roc <- list(cancer = c(), lesion = c(), SRNlesion = c(), threeway = c())

for(i in 1:length(train)){
  # Limit the variables looked at to those from the optimal model
  selected_train <- select(train[[i]], one_of(colnames(train[[i]])[1], rownames(orig_rf_opt[[i]]$importance)))
  # Pick important variables of the model
  set.seed(050416)
  impfactor_Data_List[[i]] <- Boruta(as.formula(paste(colnames(selected_train)[1], " ~ .", sep = "")), 
                                     data = selected_train, mcAdj = TRUE, maxRuns = 1000)
  # Get confirmed important variables
  confirmed_vars[[i]] <- as.data.frame(impfactor_Data_List[[i]]['finalDecision'])  %>% 
    mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)
  
  # write table to results to be used by other analysis components
  write.csv(confirmed_vars[[i]], paste("results/tables/", colnames(train[[i]])[1], "_confirmed_vars.csv", sep = ""), 
            row.names = F)
  # Select specific variables from optimum RF
  test_set <- select(train[[i]], one_of(colnames(train[[i]])[1], confirmed_vars[[i]][, 'otus']))
  # Use the selected data set in AUCRF now
  set.seed(050416)
  selected_RF_run[[i]] <- AUCRF(as.formula(paste(colnames(train[[i]])[1], " ~ .", sep = "")), 
                                data=test_set, pdel=0.05, ntree=500, ranking='MDA')
  selected_rf_opt[[i]] <- selected_RF_run[[i]]$RFopt
  selected_probs[[i]] <- predict(selected_rf_opt[[i]], type = 'prob')[, 2]
  selected_roc[[i]] <- roc(test_set[, 1] ~ selected_probs[[i]])
  # Remove remaining unuseful variables
  rm(test_set, i)
}

### Graph the ROC curves for each of the different models and test for difference

# Created needed vectors and lists
rocNameList <- list(
  SRNlesion_train_roc = orig_roc[["SRNlesion"]], SRNlesion_selected_train_roc = selected_roc[["SRNlesion"]], 
                    lesion_train_roc = orig_roc[["lesion"]], lesion_selected_train_roc = selected_roc[["lesion"]])

variableList <- c("sensitivities", "specificities")
modelList <- c("SRNlesionALL", "SRNlesionSELECT", "lesionALL", "lesionSELECT")

sens_specif_table <- makeSensSpecTable(rocNameList, variableList, modelList)
write.csv(sens_specif_table, "results/tables/ROCCurve_sens_spec.csv")

# Obtain the pvalue statistics as well as the bonferroni corrected values

corr_pvalue_ROC_table <- getROCPvalue(rocNameList, modelList, 4, multi = T)
write.csv(corr_pvalue_ROC_table, "results/tables/ROCCurve_corrected_pvalues.csv")

# Save generated data so it does not need to be saved as multiple tables or csv files
save.image("exploratory/RFwFit.RData")

# Create the graph
ggplot(sens_specif_table, aes(sensitivities, specificities)) + 
  geom_line(aes(group = model, color = model), size = 1.5) + scale_x_continuous(trans = "reverse") + 
  theme_bw() + xlab("Sensitivity") + ylab("Specificity") + geom_abline(intercept = 1, linetype = 2) + 
  theme(axis.title = element_text(face = "bold"))

ggsave(file = "results/figures/ROCCurve_withFit.tiff", width=8, height = 8, dpi = 300)

# Write out .csv for lesion variables of optimum RFmodel
write.csv(rownames(orig_rf_opt$lesion$importance), "results/tables/lesion_RFOpt_Imp_Vars.csv", row.names = F)

