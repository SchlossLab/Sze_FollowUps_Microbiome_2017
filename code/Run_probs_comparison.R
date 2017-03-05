### Prediction of Follow Ups Analysis
### How do the with and without FIT models do in correctly calling initial and follow up samples
### P-value table of comparison and Figure 4 graph
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "caret"))

# Load lesion model data
follow_up_probability <- read.csv(
  "results/tables/follow_up_probability_summary.csv", 
  header = T, stringsAsFactors = F)

red_follow_up_probability <- read.csv("results/tables/reduced_follow_up_probability_summary.csv", 
                                      stringsAsFactors = F, header = T)

# Load IF model data
IF_follow_up_probability <- read.csv(
  "results/tables/IF_follow_up_probability_summary.csv", 
  header = T, stringsAsFactors = F)

IF_red_follow_up_probability <- read.csv("results/tables/reduced_IF_follow_up_probability_summary.csv", 
                                      stringsAsFactors = F, header = T)


# Read in meta data tables
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", 
                       stringsAsFactors = F, header = T)


# create tables to hold wilcoxson paired tests pvalues with BH correction
lesion_wilcox_pvalue_summary <- getProb_PairedWilcox(follow_up_probability)
lesion_red_wilcox_pvalue_summary <- getProb_PairedWilcox(red_follow_up_probability)
IF_wilcox_pvalue_summary <- getProb_PairedWilcox(IF_follow_up_probability)
IF_red_wilcox_pvalue_summary <- getProb_PairedWilcox(IF_red_follow_up_probability)

all_wilcox_summary <- rbind(lesion_wilcox_pvalue_summary, 
                    lesion_red_wilcox_pvalue_summary, 
                    IF_wilcox_pvalue_summary, 
                    IF_red_wilcox_pvalue_summary)

all_wilcox_summary <- as.data.frame(all_wilcox_summary) %>% 
  mutate(model_type = c(rep("lesion", 4), rep("red_lesion", 4), rep("IF", 4), rep("red_IF", 4))) %>% 
  mutate(comparison = rep(c("lesion", "all_adenoma", "carcinoma_only", "SRN_only"), 4))


# Create temporary list to store used data
tempList <- list(
  lesion = follow_up_probability, 
  red_lesion = red_follow_up_probability, 
  IF = IF_follow_up_probability, 
  red_IF = IF_red_follow_up_probability
)

# Get model summary information

model_summary_info <- as.data.frame(c())
models <- c("lesion", "red_lesion", "IF", "red_IF")

for(i in 1:length(tempList)){
  
  model_summary_info <- 
    rbind(model_summary_info, 
          cbind(info = rownames(get_confusion_data(tempList[[i]], good_metaf)), 
                rbind(
                  get_confusion_data(tempList[[i]], good_metaf), 
                  get_confusion_data(filter(tempList[[i]], Dx_Bin != "cancer"), 
                                     filter(good_metaf, Dx_Bin != "cancer")), 
                  get_confusion_data(filter(tempList[[i]], Dx_Bin == "cancer"), 
                                     filter(good_metaf, Dx_Bin == "cancer"))), 
                samples_tested = c(rep("all", 18), rep("adn", 18), rep("crc", 18)), 
                model_type = rep(models[i], 54)
  ))
}


# create a new matrix to store the data
confusion_counts_summary <- c()

for(i in 1:length(tempList)){
  
  confusion_counts_summary <- rbind(
    confusion_counts_summary, 
    make_confusionTable(tempList[[i]], good_metaf), 
    make_confusionTable(tempList[[i]], good_metaf, n = 67, m = 132))

}

# Add final column to confusion counts table on model type
confusion_counts_summary <- as.data.frame(confusion_counts_summary) %>% 
  mutate(model_type = c(rep("lesion", 4), rep("red_lesion", 4), rep("IF", 4), rep("red_IF", 4)))

#Write out data tables for other use
write.csv(all_wilcox_summary, 
          "results/tables/all_models_wilcox_paired_pvalue_summary.csv", row.names = F)

write.csv(model_summary_info, 
          "results/tables/all_models_summary_info.csv", row.names = F)

write.csv(confusion_counts_summary, "results/tables/all_models_confusion_summary.csv", row.names = F)



