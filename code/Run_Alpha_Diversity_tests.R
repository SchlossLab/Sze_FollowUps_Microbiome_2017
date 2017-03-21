### Analyze Alpha diversity
### Test if there are differences in intial and follow ups based on alpha diversity
### Marc Sze

# Load in needed functions and libraries

source('code/functions.R')

loadLibs("dplyr")

# Load in needed data
good_metaf <- read.csv("data/process/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)
alpha_summary <- read.delim("data/process/final.groups.ave-std.summary", stringsAsFactors = F)

# Create data set to be tested on
alpha_data <- alpha_summary[match(c(good_metaf$initial, good_metaf$followUp), alpha_summary$group), ]

# Add custom columns from metaf and select only ones to be used for testing
alpha_data <- mutate(alpha_data, 
                     sampleType = 
                       ifelse(alpha_data$group %in% good_metaf$initial, "initial", "followups"), 
                     diagnosis = 
                       ifelse(alpha_data$group %in% filter(good_metaf, Diagnosis != "adenoma")[, "initial"] | 
                         alpha_data$group %in% filter(good_metaf, Diagnosis != "adenoma")[, "followUp"], 
                         "adenocarcinoma", "adenoma"), 
                     Dx_Bin = rep(good_metaf$Dx_Bin, 2)) %>% 
              select(group, sobs, shannon, shannoneven, sampleType, diagnosis, Dx_Bin)

# Select out specific columns for significance testing with a paired wilcoxson test
alpha_table_summary <- rbind(
  #Adenoma follow ups
  get_alpha_pvalues(filter(alpha_data, Dx_Bin == "adenoma")), 
  #SRN follow ups
  get_alpha_pvalues(filter(alpha_data, Dx_Bin == "adv_adenoma")), 
  #CRC follow ups
  get_alpha_pvalues(filter(alpha_data, Dx_Bin == "cancer"))) %>% 
  # Adjust for multiple comparisons
  mutate(BH_adj_pvalue = p.adjust(pvalue, method = "BH"))

rownames(alpha_table_summary) <- c("adn_sobs", "adn_shannon", "adn_evenness", 
                                   "srn_sobs", "srn_shannon", "srn_evenness", 
                                   "crc_sobs", "crc_shannon", "crc_evenness")

write.csv(alpha_table_summary, "data/process/tables/alpha_table_summary.csv")


