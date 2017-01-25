### Used to create figure S1
### Investigate if days has an impact on the overall outcome
## Marc Sze

# Load needed packages and functions
source('code/functions.R')

loadLibs("dplyr")

# Load needed data
thetaCompTotal <- dissplit(
  'data/process/final.thetayc.0.03.lt.ave.dist',
  split=F, meta = F)
metaI <- read.csv("results/tables/mod_metadata/metaI_final.csv", 
  stringsAsFactors = F, header = T)
metaF <- read.csv("results/tables/mod_metadata/metaF_final.csv", 
  stringsAsFactors = F, header = T)
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", 
  stringsAsFactors = F, header = T)

#Generate distance table to be used
difference_table_treatment <- pickDistanceValues(
  as.character(metaF$initial), as.character(metaF$followUp), 
  thetaCompTotal, metaF, 
  c("fit_result", "fit_followUp", "Dx_Bin", "dx", "time", "EDRN"))

# Create a summary table 
summary_statistics <- as.data.frame(rbind(
  c("Adenoma", filter(difference_table_treatment, dx == "adenoma") %>% 
      select(distance, time) %>% summarise_each(funs(mean, sd))), 
  c("Cancer", filter(difference_table_treatment, dx == "cancer") %>% 
      select(distance, time) %>% summarize_each(funs(mean, sd)))))

colnames(summary_statistics) <- c(
  "dx", "mean_thetayc_change", "mean_days", "sd_thetayc_change", "sd_days")

pvalues_summary <- rbind(
  c("thetayc_change", wilcox.test(distance ~ dx, 
    data = difference_table_treatment)$p.value), 
  c("time", wilcox.test(time ~ dx, 
    data = difference_table_treatment)$p.value))

colnames(pvalues_summary) <- c("test_values", "uncorrected_p_value")

write.csv(pvalues_summary, "results/tables/time_pvalues.csv", row.names = F)
write.csv(difference_table_treatment, "results/tables/time_datatable.csv", 
  row.names = F)

