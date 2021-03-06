### Prediction of Follow Ups Analysis
### How do the with and without FIT models do in correctly calling initial and follow up samples
### P-value table of comparison and Figure 4 graph
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "caret"))

# Load adenoma model data
adn_follow_up_probability <- read.csv("data/process/tables/adn_follow_up_probability_summary.csv", 
                                      stringsAsFactors = F, header = T)

# Load srn model data
srn_follow_up_probability <- read.csv("data/process/tables/srn_follow_up_probability_summary.csv", 
                                      stringsAsFactors = F, header = T)

# Load crc model data
crc_follow_up_probability <- read.csv("data/process/tables/crc_follow_up_probability_summary.csv", 
                                      stringsAsFactors = F, header = T)


# Read in meta data tables
good_metaf <- read.csv("data/process/mod_metadata/metaF_final.csv", 
                       stringsAsFactors = F, header = T)


# create tables to hold wilcoxson paired tests pvalues with BH correction
adn_wilcox_pvalue_summary <- getProb_PairedWilcox(
  adn_follow_up_probability, 
  rown = "adenoma", 
  not_group = c("cancer", "adv_adenoma"), 
  extra_specifics = "adenoma")

srn_wilcox_pvalue_summary <- getProb_PairedWilcox(
  srn_follow_up_probability, 
  rown = "adv_adenoma", 
  not_group = c("cancer", "adenoma"), 
  extra_specifics = "adv_adenoma")

crc_wilcox_pvalue_summary <- getProb_PairedWilcox(
  crc_follow_up_probability, 
  rown = "carcinoma", 
  not_group = "adenoma", 
  extra_specifics = "cancer")

all_wilcox_summary <- rbind(adn_wilcox_pvalue_summary, 
                            srn_wilcox_pvalue_summary, 
                            crc_wilcox_pvalue_summary) %>% 
  mutate(BH_correction = p.adjust(Pvalue, method = "BH"))

all_wilcox_summary <- as.data.frame(all_wilcox_summary) %>% 
  mutate(model_type = c("adn_m", "srn_m", "crc_m")) %>% 
  mutate(comparison = c("adenoma", "SRN", "carcinoma"))


# Create temporary list to store used data
tempList <- list(
  red_adn = adn_follow_up_probability, 
  red_srn = srn_follow_up_probability, 
  red_crc = crc_follow_up_probability
)

# Get model summary information

model_summary_info <- as.data.frame(c())

model_summary_info <- 
  cbind(rbind(
                get_confusion_data(tempList[["red_adn"]], good_metaf, to_filter1 = "cancer", 
                                   to_filter2 = "adv_adenoma"), 
                get_confusion_data(tempList[["red_srn"]], good_metaf, to_filter1 = "cancer", 
                                   to_filter2 = "adenoma"), 
                get_confusion_data(tempList[["red_crc"]], good_metaf, to_filter1 = "adenoma", 
                                   to_filter2 = "adv_adenoma", column = "Dx_Bin", DF_data = "Yes")), 
              samples_tested = c(rep("adn", 18), rep("SRN", 18), rep("crc", 18)), 
              model_type = c(rep("red_adn", 18), rep("red_SRN", 18), rep("red_crc", 18)))




#Write out data tables for other use
write.csv(all_wilcox_summary, 
          "data/process/tables/all_crc_srn_adn_models_wilcox_paired_pvalue_summary.csv", row.names = F)

write.csv(model_summary_info, 
          "data/process/tables/all_crc_srn_adn_models_summary_info.csv")



