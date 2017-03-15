### Run to check probabilities of chemo and rads effects
### Could these treatments influence findings of surgical removal
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr"))

# Read in data tables

data_list <- list(
  red_lesion_probs = read.csv("data/process/tables/reduced_follow_up_probability_summary.csv", 
                              header = T, stringsAsFactors = F), 
  red_IF_probs = read.csv("data/process/tables/reduced_IF_follow_up_probability_summary.csv", 
                          header = T, stringsAsFactors = F)
)

difference_table_treatment <- read.csv("data/process/tables/difference_table.csv", 
                                       header = T, stringsAsFactors = F)

good_metaf <- read.csv("data/process/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)
alpha_summary <- read.delim("data/process/final.groups.ave-std.summary", stringsAsFactors = F)

# Create data set more amenable to changes
alpha_data <- alpha_summary[match(c(good_metaf$initial, good_metaf$followUp), alpha_summary$group), ]

# Add custom columns from metaf and select only ones to be used for testing
alpha_data <- mutate(alpha_data, 
                     sampleType = 
                       ifelse(alpha_data$group %in% good_metaf$initial, "initial", "followups"), 
                     diagnosis = 
                       ifelse(alpha_data$group %in% filter(good_metaf, Diagnosis != "adenoma")[, "initial"] | 
                                alpha_data$group %in% filter(good_metaf, Diagnosis != "adenoma")[, "followUp"], 
                              "adenocarcinoma", "adenoma"), 
                     chemo = rep(good_metaf$chemo_received, 2), 
                     rads = rep(good_metaf$radiation_received, 2)) %>% 
  select(group, sobs, shannon, shannoneven, sampleType, diagnosis, chemo, rads)

# Create common data table to generate mean, sd, and testing on
test <- (filter(alpha_data, sampleType == "initial")[, c("sobs", "shannon", "shannoneven")] - 
  filter(alpha_data, sampleType == "followups")[, c("sobs", "shannon", "shannoneven")]) %>% 
  mutate(distance = difference_table_treatment$distance, 
         red_lesion = (filter(data_list[["red_lesion_probs"]], sampleType == "initial")[, "Yes"] - 
           filter(data_list[["red_lesion_probs"]], sampleType == "followup")[, "Yes"]), 
         red_IF = (filter(data_list[["red_IF_probs"]], sampleType == "initial")[, "Yes"] - 
                         filter(data_list[["red_IF_probs"]], sampleType == "followup")[, "Yes"]), 
         chemo = good_metaf$chemo_received, 
         rads = good_metaf$radiation_received, 
         Dx_Bin = good_metaf$Dx_Bin)

# Create variables for actual analysis to be automated
tests <- c("sobs", "sobs", "shannon", "shannon", "shannoneven", "shannoneven", 
           "distance", "distance", "red_lesion", "red_lesion", "red_IF", "red_IF")

mean_table <- as.data.frame(matrix(nrow = 12, ncol = 6, dimnames = list(
  rown = c("chemo_sobs", "rads_sobs", "chemo_shannon", "rads_shannon", "chem_even", "rads_even", 
           "chemo_dist", "rads_dist", "chemo_les", "rad_les", "chemo_IF", "rad_IF"), 
  coln = c("chemo_rad_mean", "removal_only_mean", "cr_sd", "ro_sd", "pvalue", "bh"))))

# Run automated tests
for(i in 1:length(tests)){
  if(i %% 2 != 0){
    
    mean_table[i, "chemo_rad_mean"] <- mean(filter(test, chemo == "yes")[, tests[i]])
    mean_table[i, "removal_only_mean"] <- mean(filter(test, chemo == "no")[, tests[i]])
    mean_table[i, "cr_sd"] <- sd(filter(test, chemo == "yes")[, tests[i]])
    mean_table[i, "ro_sd"] <- sd(filter(test, chemo == "no")[, tests[i]])
    
    mean_table[i, "pvalue"] <- wilcox.test(filter(test, chemo == "yes")[, tests[i]], 
                                          filter(test, chemo == "no")[, tests[i]])$p.value
    
  } else{
    
    mean_table[i, "chemo_rad_mean"] <- mean(filter(test, rads == "yes")[, tests[i]])
    mean_table[i, "removal_only_mean"] <- mean(filter(test, rads == "no")[, tests[i]])
    mean_table[i, "cr_sd"] <- sd(filter(test, rads == "yes")[, tests[i]])
    mean_table[i, "ro_sd"] <- sd(filter(test, rads == "no")[, tests[i]])
    
    mean_table[i, "pvalue"] <- wilcox.test(filter(test, chemo == "yes")[, tests[i]], 
                                           filter(test, chemo == "no")[, tests[i]])$p.value
  }

}

# Implement P-value correction
mean_table$bh <- p.adjust(mean_table$pvalue, method = "BH")

# Test whether radiation is significantly different then chemo in IF proability reduction
mean_table <- rbind(mean_table, 
                    chemo_v_rads_IF = c(NA, NA, NA, NA, 
                                        wilcox.test(filter(test, chemo == "yes" & rads == "no")[, "red_IF"], 
                                                  filter(test, rads == "yes")[, "red_IF"])$p.value, NA))

# Write out final table
write.csv(mean_table, "data/process/tables/probs_chemo_rad_pvalue_summary.csv")
write.csv(test, "data/process/tables/chemo_rad_summary.csv")


