### Analyze Alpha diversity
### Test if there are differences in intial and follow ups based on alpha diversity
### Marc Sze

# Load in needed functions and libraries

source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr"))


# Load in needed data
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)
alpha_summary <- read.delim("data/process/final.groups.ave-std.summary", stringsAsFactors = F)

# Create data set to be tested on
alpha_data <- alpha_summary[match(c(good_metaf$initial, good_metaf$followUp), alpha_summary$group), ]

# Add custom columns from metaf
alpha_data$sampleType[alpha_data$group %in% good_metaf$initial] <- "initial"
alpha_data$sampleType[alpha_data$group %in% good_metaf$followUp] <- "followups"

alpha_data$diagnosis[alpha_data$group %in% filter(good_metaf, Diagnosis != "adenoma")[, "initial"]] <- "adenomocarcinoma"
alpha_data$diagnosis[alpha_data$group %in% filter(good_metaf, Diagnosis != "adenoma")[, "followUp"]] <- "adenomaocarcinoma"
alpha_data$diagnosis[alpha_data$group %in% filter(good_metaf, Diagnosis == "adenoma")[, "initial"]] <- "adenoma"
alpha_data$diagnosis[alpha_data$group %in% filter(good_metaf, Diagnosis == "adenoma")[, "followUp"]] <- "adenoma"

# Select out specific columns for significance testing with a paired wilcoxson test
alpha_data <- select(alpha_data, group, sobs, shannon, shannoneven, sampleType, diagnosis)


All_followups <- get_alpha_pvalues(alpha_data, rows_names = c("all_sobs", "all_shannon", "all_evenness"))

Adn_followups <- get_alpha_pvalues(filter(alpha_data, diagnosis == "adenoma"), rows_names = c("adn_sobs", "adn_shannon", "adn_evenness"))

CRC_followups <- get_alpha_pvalues(filter(alpha_data, diagnosis != "adenoma"), rows_names = c("crc_sobs", "crc_shannon", "crc_evenness"))


alpha_table_summary <- rbind(All_followups, Adn_followups, CRC_followups)

write.csv(alpha_table_summary, "results/tables/alpha_table_summary.csv")


