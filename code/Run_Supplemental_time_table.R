### Used to create figure S1
### Investigate if days has an impact on the overall outcome
## Marc Sze

# Load needed packages and functions
source('code/functions.R')

loadLibs("dplyr")

# Load needed data
thetaCompTotal <- read.dist('data/process/final.thetayc.0.03.lt.ave.dist')
metaI <- read.csv("data/process/mod_metadata/metaI_final.csv", 
  stringsAsFactors = F, header = T)
metaF <- read.csv("data/process/mod_metadata/metaF_final.csv", 
  stringsAsFactors = F, header = T)
good_metaf <- read.csv("data/process/mod_metadata/good_metaf_final.csv", 
  stringsAsFactors = F, header = T)

#Generate distance table to be used
difference_table_treatment <- pickDistanceValues(
  as.character(metaF$initial), as.character(metaF$followUp), 
  thetaCompTotal, metaF, 
  c("fit_result", "fit_followUp", "Dx_Bin", "dx", "time", "EDRN"))

# Create a summary table 
summary_statistics <- bind_rows(
  c(filter(difference_table_treatment, Dx_Bin == "adenoma") %>% 
      select(distance, time) %>% summarise_each(funs(mean, sd))), 
  c(filter(difference_table_treatment, Dx_Bin == "adv_adenoma") %>% 
      select(distance, time) %>% summarise_each(funs(mean, sd))), 
  c(filter(difference_table_treatment, Dx_Bin == "cancer") %>% 
      select(distance, time) %>% summarize_each(funs(mean, sd)))) 

colnames(summary_statistics) <- c(
  "mean_thetayc_change", "mean_days", "sd_thetayc_change", "sd_days")

rownames(summary_statistics) <- c("adenoma", "srn", "carcinoma")


pvalues_summary <- as.data.frame(
  matrix(nrow = 3, ncol = 3, 
         dimnames = list(rown = c(), coln = c("comparison", "pvalue", "bh")))) %>% 
  mutate(
    comparison = c("adnVsrn", "adnVcrc", "srnVcrc"), 
    pvalue = c(
    wilcox.test(
      filter(difference_table_treatment, Dx_Bin == "adenoma")[, "time"], 
      filter(difference_table_treatment, Dx_Bin == "adv_adenoma")[, "time"])$p.value, 
    wilcox.test(
      filter(difference_table_treatment, Dx_Bin == "adenoma")[, "time"], 
      filter(difference_table_treatment, Dx_Bin == "cancer")[, "time"])$p.value, 
    wilcox.test(
      filter(difference_table_treatment, Dx_Bin == "adv_adenoma")[, "time"], 
      filter(difference_table_treatment, Dx_Bin == "cancer")[, "time"])$p.value), 
    bh = p.adjust(pvalue, method = "BH"))


write.csv(pvalues_summary, "data/process/tables/time_pvalues.csv", row.names = F)
write.csv(summary_statistics, "data/process/tables/time_summary_data.csv")
write.csv(difference_table_treatment, "data/process/tables/time_datatable.csv", 
  row.names = F)

