### Run to check probabilities of chemo and rads effects
### Could these treatments influence findings of surgical removal
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr"))

# Read in data tables

data_list <- list(
  lesion_probs = read.csv("results/tables/follow_up_probability_summary.csv", header = T, stringsAsFactors = F), 
  red_lesion_probs = read.csv("results/tables/reduced_follow_up_probability_summary.csv", 
                              header = T, stringsAsFactors = F), 
  IF_probs = read.csv("results/tables/IF_follow_up_probability_summary.csv", 
                      header = T, stringsAsFactors = F), 
  red_IF_probs = read.csv("results/tables/reduced_IF_follow_up_probability_summary.csv", 
                          header = T, stringsAsFactors = F)
)


# Test for differences in positive probability changes
pvalue_table <- as.data.frame(matrix(nrow = 16, ncol = 2, dimnames = list(
  rown = c("chemo_les", "chemo_red_les", "chemo_IF", "chemo_red_IF", 
           "rad_les", "rad_red_les", "rad_IF", "rad_red_IF", "chemo_b", "rad_b", 
           "chemo_sobs", "chemo_shannon", "chemo_shannoneven", 
           "rads_sobs", "rads_shannon", "rads_shannoneven"), coln = c("pvalue", "bh"))))

for(i in 1:length(data_list)){
  
  pvalue_table[i, "pvalue"] <- wilcox.test(
    filter(data_list[[i]], sampleType == "initial", chemo == "yes")[, "Yes"] - 
                filter(data_list[[i]], sampleType == "followup", chemo == "yes")[, "Yes"], 
              filter(data_list[[i]], sampleType == "initial", chemo == "no")[, "Yes"] - 
                filter(data_list[[i]], sampleType == "followup", chemo == "no")[, "Yes"])$p.value
  
  pvalue_table[i+4, "pvalue"] <- wilcox.test(
    filter(data_list[[i]], sampleType == "initial", rads == "yes")[, "Yes"] - 
                filter(data_list[[i]], sampleType == "followup", rads == "yes")[, "Yes"], 
              filter(data_list[[i]], sampleType == "initial", rads == "no")[, "Yes"] - 
                filter(data_list[[i]], sampleType == "followup", rads == "no")[, "Yes"])$p.value
}


# Test thetayc value differences between treatment type

difference_table_treatment <- read.csv("results/tables/difference_table.csv", 
                                       header = T, stringsAsFactors = F) %>% 
  filter(!is.na(fit_followUp)) %>% 
  mutate(chemo = data_list[["lesion_probs"]][1:66, "chemo"], 
         rads = data_list[["lesion_probs"]][1:66, "rads"])

pvalue_table["chemo_b", "pvalue"] <- wilcox.test(
  filter(difference_table_treatment, chemo == "yes")[, "distance"], 
  filter(difference_table_treatment, chemo == "no")[, "distance"])$p.value

pvalue_table["rad_b", "pvalue"] <- wilcox.test(
  filter(difference_table_treatment, rads == "yes")[, "distance"], 
  filter(difference_table_treatment, rads == "no")[, "distance"])$p.value


# Test alpha diversity metrics by treatment type
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)
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
                     chemo = rep(good_metaf$chemo_received, 2), 
                     rads = rep(good_metaf$radiation_received, 2)) %>% 
  select(group, sobs, shannon, shannoneven, sampleType, diagnosis, chemo, rads)


alpha_to_test <- c("sobs", "shannon", "shannoneven")

for(i in 1:length(alpha_to_test)){
 
  pvalue_table[10+i, "pvalue"] <- wilcox.test(
    filter(alpha_data, sampleType == "initial", chemo == "yes")[, alpha_to_test[i]] - 
      filter(alpha_data, sampleType == "followups", chemo == "yes")[, alpha_to_test[i]], 
    filter(alpha_data, sampleType == "initial", chemo == "no")[, alpha_to_test[i]] - 
      filter(alpha_data, sampleType == "followups", chemo == "no")[, alpha_to_test[i]])$p.value
  
  pvalue_table[13+i, "pvalue"] <- wilcox.test(
    filter(alpha_data, sampleType == "initial", rads == "yes")[, alpha_to_test[i]] - 
      filter(alpha_data, sampleType == "followups", rads == "yes")[, alpha_to_test[i]], 
    filter(alpha_data, sampleType == "initial", rads == "no")[, alpha_to_test[i]] - 
      filter(alpha_data, sampleType == "followups", rads == "no")[, alpha_to_test[i]])$p.value
  
}

# Implement P-value correction
pvalue_table$bh <- p.adjust(pvalue_table$pvalue, method = "BH")

write.csv(pvalue_table, "results/tables/probs_chemo_rad_pvalue_summary.csv")


