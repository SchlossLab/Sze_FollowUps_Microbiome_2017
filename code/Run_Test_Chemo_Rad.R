### Run to check probabilities of chemo and rads effects
### Could these treatments influence findings of surgical removal
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

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


pvalue_table <- as.data.frame(matrix(nrow = 8, ncol = 2, dimnames = list(
  rown = c("chemo_les", "chemo_red_les", "chemo_IF", "chemo_red_IF", 
           "rad_les", "rad_red_les", "rad_IF", "rad_red_IF"), coln = c("pvalue", "bh"))))

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


pvalue_table$bh <- p.adjust(pvalue_table$pvalue, method = "BH")

write.csv(pvalue_table, "results/tables/probs_chemo_rad_pvalue_summary.csv")


