### Run to check probabilities of chemo and rads effects
### Could these treatments influence findings of surgical removal
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr"))

# Read in data tables

data_list <- list(
  red_adn_probs = read.csv("data/process/tables/adn_reduced_follow_up_probability_summary.csv", 
                              header = T, stringsAsFactors = F), 
  red_srn_probs = read.csv("data/process/tables/srn_reduced_follow_up_probability_summary.csv", 
                           header = T, stringsAsFactors = F), 
  red_crc_probs = read.csv("data/process/tables/crc_reduced_follow_up_probability_summary.csv", 
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
                     rads = rep(good_metaf$radiation_received, 2), 
                     Dx_Bin = rep(good_metaf$Dx_Bin, 2)) %>% 
  select(group, sobs, shannon, shannoneven, sampleType, diagnosis, chemo, rads, Dx_Bin)

# Create common data table to generate mean, sd, and testing on for adn and crc
adn_test <- (filter(alpha_data, sampleType == "initial" & Dx_Bin != "cancer")[, c("sobs", "shannon", "shannoneven")] - 
  filter(alpha_data, sampleType == "followups" & Dx_Bin != "cancer")[, c("sobs", "shannon", "shannoneven")]) %>% 
  mutate(distance = filter(difference_table_treatment, Dx_Bin != "cancer")[, "distance"], 
         red_adn = (filter(data_list[["red_adn_probs"]], sampleType == "initial" & Dx_Bin != "cancer")[, "Yes"] - 
           filter(data_list[["red_adn_probs"]], sampleType == "followup" & Dx_Bin != "cancer")[, "Yes"]), 
         chemo = filter(good_metaf, Dx_Bin != "cancer")[, "chemo_received"], 
         rads = filter(good_metaf, Dx_Bin != "cancer")[, "radiation_received"], 
         surgery = filter(good_metaf, Dx_Bin !="cancer")[, "Surgery"], 
         Dx_Bin = filter(good_metaf, Dx_Bin != "cancer")[, "Dx_Bin"])


crc_test <- (filter(alpha_data, sampleType == "initial" & Dx_Bin == "cancer")[, c("sobs", "shannon", "shannoneven")] - 
               filter(alpha_data, sampleType == "followups" & Dx_Bin == "cancer")[, c("sobs", "shannon", "shannoneven")]) %>% 
  mutate(distance = filter(difference_table_treatment, Dx_Bin == "cancer")[, "distance"], 
         red_crc = (filter(data_list[["red_crc_probs"]], sampleType == "initial")[, "Yes"] - 
                      filter(data_list[["red_crc_probs"]], sampleType == "followup")[, "Yes"]), 
         chemo = filter(good_metaf, Dx_Bin == "cancer")[, "chemo_received"], 
         rads = filter(good_metaf, Dx_Bin == "cancer")[, "radiation_received"], 
         surgery = filter(good_metaf, Dx_Bin =="cancer")[, "Surgery"], 
         Dx_Bin = filter(good_metaf, Dx_Bin == "cancer")[, "Dx_Bin"])


# Create variables for actual analysis to be automated
adn_tests <- c("sobs", "shannon", "shannoneven", "distance", "red_adn")

crc_tests <- c("sobs", "sobs", "shannon", "shannon", "shannoneven", "shannoneven", 
               "distance", "distance", "red_crc", "red_crc")

adn_mean_table <- as.data.frame(matrix(nrow = 5, ncol = 6, dimnames = list(
  rown = c("sobs", "shannon", "even", "dist", "red_adn"), 
  coln = c("surg_mean", "no_surg_mean", "surg_sd", "no_surg_sd", "pvalue", "bh"))))

crc_mean_table <- as.data.frame(matrix(nrow = 10, ncol = 6, dimnames = list(
  rown = c("chemo_sobs", "rads_sobs", "chemo_shannon", "rads_shannon", "chem_even", "rads_even", 
           "chemo_dist", "rads_dist", "chemo_red_crc", "rad_red_crc"), 
  coln = c("chemo_rad_mean", "removal_only_mean", "cr_sd", "ro_sd", "pvalue", "bh"))))

# Run automated tests for adn
for(i in 1:length(adn_tests)){
  
  adn_mean_table[i, "surg_mean"] <- mean(filter(adn_test, surgery == "Y")[, adn_tests[i]])
  adn_mean_table[i, "no_surg_mean"] <- mean(filter(adn_test, surgery == "N")[, adn_tests[i]])
  adn_mean_table[i, "surg_sd"] <- sd(filter(adn_test, surgery == "Y")[, adn_tests[i]])
  adn_mean_table[i, "no_surg_sd"] <- sd(filter(adn_test, surgery == "N")[, adn_tests[i]])
  
  adn_mean_table[i, "pvalue"] <- wilcox.test(filter(adn_test, surgery == "Y")[, adn_tests[i]], 
                                         filter(adn_test, surgery == "N")[, adn_tests[i]])$p.value
}

# Implement P-value correction
adn_mean_table$bh <- p.adjust(adn_mean_table$pvalue, method = "BH")


# Run automated tests for crc
for(i in 1:length(crc_tests)){
  if(i %% 2 != 0){
    
    crc_mean_table[i, "chemo_rad_mean"] <- mean(filter(crc_test, chemo == "yes")[, crc_tests[i]])
    crc_mean_table[i, "removal_only_mean"] <- mean(filter(crc_test, chemo == "no")[, crc_tests[i]])
    crc_mean_table[i, "cr_sd"] <- sd(filter(crc_test, chemo == "yes")[, crc_tests[i]])
    crc_mean_table[i, "ro_sd"] <- sd(filter(crc_test, chemo == "no")[, crc_tests[i]])
    
    crc_mean_table[i, "pvalue"] <- wilcox.test(filter(crc_test, chemo == "yes")[, crc_tests[i]], 
                                          filter(crc_test, chemo == "no")[, crc_tests[i]])$p.value
    
  } else{
    
    crc_mean_table[i, "chemo_rad_mean"] <- mean(filter(crc_test, rads == "yes")[, crc_tests[i]])
    crc_mean_table[i, "removal_only_mean"] <- mean(filter(crc_test, rads == "no")[, crc_tests[i]])
    crc_mean_table[i, "cr_sd"] <- sd(filter(crc_test, rads == "yes")[, crc_tests[i]])
    crc_mean_table[i, "ro_sd"] <- sd(filter(crc_test, rads == "no")[, crc_tests[i]])
    
    crc_mean_table[i, "pvalue"] <- wilcox.test(filter(crc_test, chemo == "yes")[, crc_tests[i]], 
                                           filter(crc_test, chemo == "no")[, crc_tests[i]])$p.value
  }

}

# Implement P-value correction
crc_mean_table$bh <- p.adjust(crc_mean_table$pvalue, method = "BH")

# Test whether radiation is significantly different then chemo in IF proability reduction
crc_mean_table <- rbind(crc_mean_table, 
                    chemo_v_rads_red = c(NA, NA, NA, NA, 
                                        wilcox.test(filter(crc_test, chemo == "yes" & rads == "no")[, "red_crc"], 
                                                  filter(crc_test, rads == "yes")[, "red_crc"])$p.value, NA))


# Run fisher exact test on surgery proportions for adn and srn
test_data <- matrix(nrow = 2, ncol = 2, dimnames = list(rown = c("adn", "srn"), coln = c("surg_Y", "surg_N")))

test_data[, "surg_Y"] <- c(length(rownames(filter(good_metaf, Dx_Bin == "adenoma" & Surgery == "Y"))), 
                           length(rownames(filter(good_metaf, Dx_Bin == "adv_adenoma" & Surgery == "Y"))))

test_data[, "surg_N"] <- c(length(rownames(filter(good_metaf, Dx_Bin == "adenoma" & Surgery == "N"))), 
                           length(rownames(filter(good_metaf, Dx_Bin == "adv_adenoma" & Surgery == "N"))))

adn_mean_table <- rbind(adn_mean_table, 
                        prop_surg_adn_srn = c(NA, NA, NA, NA, 
                                              fisher.test(test_data)$p.value, NA))



# Write out final table
write.csv(crc_mean_table, "data/process/tables/crc_probs_chemo_rad_pvalue_summary.csv")
write.csv(adn_mean_table, "data/process/tables/adn_combined_probs_surgery_pvalue_summary.csv")
write.csv(crc_test, "data/process/tables/crc_chemo_rad_summary.csv")
write.csv(adn_test, "data/process/tables/adn_combined_surgery_summary.csv")


