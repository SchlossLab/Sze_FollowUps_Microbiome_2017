### Prediction of Follow Ups Analysis
### How do the adenomas and SRN specifically compare against each other
### P-value tables for results of SRN versus Adenoma in with and without FIT models
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))

# Read in data tables
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", 
                       stringsAsFactors = F, header = T)
rf_prediction_summary <- read.csv(
  'results/tables/rf_prediction_summary.csv', header = T)
aucrf_model_cutoff <- read.csv('results/tables/aucrf_model_cutoffs.csv', 
                               header = T, row.names = 1)


# Compare differences between SRN and adenoma positive probability decrease

# Add an extra column for definitions to the rf summary table
rf_prediction_summary$Dx_Bin <- rep(good_metaf$Dx_Bin, 4)
rf_prediction_summary$EDRN <- rep(good_metaf$EDRN, 4)

# Create variable vectors to cycle through during the for loop
model_used <- c("wfit", "wofit")
samples <- c("initial", "followup")

srn_compare_pvalue_summary <- as.data.frame(matrix(
  nrow = 2, ncol = 4, dimnames = list(
    rows = c("wfit", "wofit"), 
    cols = c("prob_init", "prob_follow", "change", "proportion"))))

# Difference in probabilities between SRN and Adenoma Initials and Follow Ups
for(i in 1:length(model_used)){
  
  for(j in 1:length(samples)){
    
    srn_compare_pvalue_summary[i, j] <- wilcox.test(
      filter(rf_prediction_summary, 
        model == paste(model_used[i]), 
        sample_type == paste(samples[j]), 
        Dx_Bin == "adv_adenoma")[, "postive_probability"], 
      filter(rf_prediction_summary, 
        model == paste(model_used[i]), 
        sample_type == paste(samples[j]), 
        Dx_Bin == "adenoma")[, "postive_probability"])$p.value
  }
}

#Difference in change in probabilities between SRN and Adenoma

for(i in 1:length(model_used)){
  
  srn_compare_pvalue_summary[i, "change"] <- wilcox.test(
    filter(rf_prediction_summary, 
      model == paste(model_used[i]), 
      sample_type == "initial", 
      Dx_Bin == "adv_adenoma")[, "postive_probability"] - 
    filter(rf_prediction_summary, 
      model == paste(model_used[i]), 
      sample_type == "followup", 
      Dx_Bin == "adv_adenoma")[, "postive_probability"], 
    filter(rf_prediction_summary, 
      model == paste(model_used[i]), 
      sample_type == "initial", 
      Dx_Bin == "adenoma")[, "postive_probability"] - 
    filter(rf_prediction_summary, 
      model == paste(model_used[i]), 
      sample_type == "followup", 
      Dx_Bin == "adenoma")[, "postive_probability"])$p.value
}


#Difference in proportion above and below cutoff between SRN and Adenoma
for(i in 1:length(model_used)){
  
  srn_compare_pvalue_summary[i, "proportion"] <- 
    makeANDfish_2by2(
      select(rf_prediction_summary, -diagnosis) %>% 
      rename(diagnosis = Dx_Bin), c("initial", "followup"), 
      c("Yes", "No"), aucrf_model_cutoff, 
      model = FALSE, model_sample_type = NULL, 
      model_type = paste(model_used[i]), remove_sample = "cancer")
}


# Write out updated rf_prediction table
# Write out SRN comparison P-value summary

write.csv(rf_prediction_summary, 
          "results/tables/rf_prediction_summary.csv", row.names = FALSE)
write.csv(srn_compare_pvalue_summary, 
          "results/tables/adn_vs_srn_pvalue_summary.csv")



