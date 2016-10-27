### Random Forest analysis
### With and without Fit train and test.
### Extrapolation of Important Vairables in Optimum Random Forest
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("pROC","randomForest","AUCRF", "dplyr", "tidyr", "ggplot2", 
  "reshape2", "gridExtra", "scales", "wesanderson"))


# Read in data tables
metaI <- read.csv("results/tables/mod_metadata/metaI_final.csv", 
  stringsAsFactors = F, header = T)
metaF <- read.csv("results/tables/mod_metadata/metaF_final.csv", 
  stringsAsFactors = F, header = T)
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", 
  stringsAsFactors = F, header = T)
shared <- read.delim('data/process/final.0.03.subsample.shared', 
                     header=T, sep='\t')

# createList with all data tables stored as a list

train <- list(
  # Add training set with fit
  wFit = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(lesion, fit_result, contains("Otu0")) %>% 
                na.omit() %>% mutate(lesion = factor(lesion)), 
  # Add training set without fit            
  woFit = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(lesion, contains("Otu0")) %>% na.omit() %>% 
                mutate(lesion = factor(lesion)))


# Set up initial lists to store the most important components of the training AUCRF run
train_RF_run <- list(wfit = c(), wofit = c())
train_rf_opt <- list(wfit = c(), wofit = c())
train_probs <- list(wfit = c(), wofit = c())
train_rocs <- list(wfit = c(), wofit = c())

# Run the actual AUCRF for each different model
for(i in 1:length(train)){
  
  set.seed(050416)
  train_RF_run[[i]] <- AUCRF(as.formula(
    paste(colnames(train[[i]])[1], " ~ .", sep = "")), 
  data=train[[i]], pdel=0.05, ntree=500, ranking='MDA')

  train_rf_opt[[i]] <- train_RF_run[[i]]$RFopt
  train_probs[[i]] <- predict(train_rf_opt[[i]], type = 'prob')[, 2]
  train_rocs[[i]] <- roc(train[[i]][, 1] ~ train_probs[[i]])
  
}  

# Get threshold used by each of the different models to make the call and write table for use later

cutoffs <- as.data.frame.list(lapply(train_rocs, 
  function(y) coords(y, x='best', ret='thr'))) 
#take the threshold point used for the best sensitivity and specificty

write.csv(cutoffs, "results/tables/aucrf_model_cutoffs.csv")


# Create data frames to be used for initial and follow up samples
initials <- list(
  wfit = inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
    select(lesion, fit_result, contains("Otu0")), 
  
  wofit = inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
    select(lesion, contains("Otu0")))


followups <- list(
  wfit = inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
    select(lesionf, fit_followUp, contains("Otu0")) %>% 
    rename(lesion = lesionf) %>% rename(fit_result = fit_followUp), 
  
  wofit = inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
    select(lesionf, contains("Otu0")) %>% rename(lesion = lesionf))
  

# Get the prediction tables for each group using initial and follow up data with and without Fit

initial_predictions <- list(wfit = c(), wofit = c())
followup_predictions <- list(wfit = c(), wofit = c())

for(i in 1:length(followups)){
  
  initial_predictions[[i]] <- predict(train_rf_opt[[i]], 
    newdata = initials[[i]], type='prob')
  
  followup_predictions[[i]] <- predict(train_rf_opt[[i]], 
    newdata = followups[[i]], type='prob')
}


#Create Test set ROCs for AUCs 
#The 2 represents the specificity

test_rocs <- list(
  wfit = roc(c(initials$wfit$lesion, followups$wfit$lesion) ~ 
               c(initial_predictions$wfit[, 2], 
                followup_predictions$wfit[, 2])), 
  
  wofit = roc(c(initials$wofit$lesion, followups$wofit$lesion) ~ 
                c(initial_predictions$wofit[, 2], 
                  followup_predictions$wofit[, 2])))


# Creat needed vectors and lists to create a sensitivity and specificity table
rocNameList <- list(
  train_wfit = train_rocs[["wfit"]], train_wofit = train_rocs[["wofit"]], 
  test_wfit = test_rocs[["wfit"]], test_wofit = test_rocs[["wofit"]])

variableList <- c("sensitivities", "specificities")
modelList <- c("train_wfit", "train_wofit", "test_wfit", "test_wofit")

# Create the sensitivity and specifiicty table
sens_specif_table <- makeSensSpecTable(rocNameList, variableList, modelList)

# Obtain the pvalue statistics
corr_pvalue_ROC_table <- getROCPvalue(rocNameList, modelList, 4, multi = T)

# Make a table with train and test AUCs with cutoffs used
models_summary <- cbind(
  model = modelList, 
  AUC = c(train_rocs$wfit$auc, train_rocs$wofit$auc, 
    test_rocs$wfit$auc, test_rocs$wofit$auc), 
  cutoffs = rep(c(cutoffs$wfit, cutoffs$wofit), 2))

# Make a table with important variables for with and without fit models for later use

rf_imp_vars_summary <- rbind(
  cbind(imp_variable = train_RF_run$wfit$Xopt, 
    model = rep("wfit", length(train_RF_run$wfit$Xopt))), 
  cbind(imp_variable = train_RF_run$wofit$Xopt, 
    model = rep("wofit", length(train_RF_run$wofit$Xopt))))

### Graph the ROC curves for each of the different models and test for difference

roc_curve_graph <- ggplot(sens_specif_table, 
  aes(sensitivities, specificities)) + 
  geom_line(aes(group = model, color = model), size = 1.5) + 
  scale_x_continuous(trans = "reverse") + 
  theme_bw() + xlab("Sensitivity") + ylab("Specificity") + 
  geom_abline(intercept = 1, linetype = 2) + 
  scale_color_discrete(name = "Lesion Model", 
                       breaks = c("train_wfit", "train_wofit", "test_wfit", "test_wofit"), 
                       labels = c(paste("with Fit Train Set\n(AUC = ", 
                                        round(train_rocs$wfit$auc, 
                                          digits = 3), ")", sep = ""), 
                                  paste("without Fit Train Set\n(AUC = ", 
                                        round(train_rocs$wofit$auc, 
                                          digits = 3), ")", sep=""), 
                                  paste("with Fit Test Set\n(AUC = ", 
                                        round(test_rocs$wfit$auc, 
                                          digits = 3), ")", sep=""), 
                                  paste("without Fit Test Set\n(AUC = ", 
                                        round(test_rocs$wofit$auc, 
                                          digits = 3), ")", sep=""))) + 
  theme(axis.title = element_text(face = "bold"), 
    legend.position = c(0.75, 0.4), 
        legend.title = element_text(face="bold"))


# Save graph image as a pdf
ggsave(file = "results/figures/Figure3.pdf", roc_curve_graph, 
       width=8, height = 8, dpi = 300)


# Write out relevant tables to be used in manuscript
write.csv(models_summary, "results/tables/models_summary.csv")
write.csv(corr_pvalue_ROC_table, 
  "results/tables/roc_corr_pvalue_summary.csv")
write.csv(rf_imp_vars_summary, 
  "results/tables/rf_imp_vars_summary.csv", row.names = F)





