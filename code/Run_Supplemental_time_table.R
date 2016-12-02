### Used to create figure S1
### Investigate if days has an impact on the overall outcome
## Marc Sze


source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "vegan"))


# Load needed data
thetaCompTotal <- dissplit('data/process/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist',split=F, meta = F)
metaI <- read.csv("results/tables/mod_metadata/metaI_final.csv", stringsAsFactors = F, header = T)
metaF <- read.csv("results/tables/mod_metadata/metaF_final.csv", stringsAsFactors = F, header = T)
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)


# Load needed data
thetaCompTotal <- dissplit('data/process/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist',split=F, meta = F)
metaI <- read.csv("results/tables/mod_metadata/metaI_final.csv", stringsAsFactors = F, header = T)
metaF <- read.csv("results/tables/mod_metadata/metaF_final.csv", stringsAsFactors = F, header = T)
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)

#Generate distance table to be used
difference_table_treatment <- pickDistanceValues(
  as.character(metaF$initial), as.character(metaF$followUp), thetaCompTotal, metaF, 
  c("fit_result", "fit_followUp", "Dx_Bin", "dx", "time", "EDRN"))

# Plot the change in distance versus time
time_graph <- ggplot(difference_table_treatment, aes(x = time, y = distance)) + 
  geom_point(aes(color = dx), size = 4) + theme_bw() + 
  scale_color_manual(name = "Lesion Type", values = wes_palette("GrandBudapest"), 
                     breaks = c("adenoma", "cancer"), labels = c("Adenoma", "Cancer")) + 
  coord_cartesian(ylim = c(0, 1)) + ylab("Thetayc Distance") + xlab("Time (Days)") + 
  theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"))

# Create a summary table 
summary_statistics <- rbind(
  c("Adenoma", mean(filter(difference_table_treatment, dx == "adenoma")[, "distance"]), 
    sd(filter(difference_table_treatment, dx == "adenoma")[, "distance"]), 
    mean(filter(difference_table_treatment, dx == "adenoma")[, "time"]), 
    sd(filter(difference_table_treatment, dx == "adenoma")[, "time"])), 
  c("Cancer", mean(filter(difference_table_treatment, dx == "cancer")[, "distance"]), 
    sd(filter(difference_table_treatment, dx == "cancer")[, "distance"]), 
    mean(filter(difference_table_treatment, dx == "cancer")[, "time"]), 
    sd(filter(difference_table_treatment, dx == "cancer")[, "time"])))

colnames(summary_statistics) <- c("dx", "mean_thetayc_change", "sd_thetayc_change", "mean_days", "sd_days")

pvalues_summary <- rbind(c("thetayc_change", wilcox.test(distance ~ dx, data = difference_table_treatment)$p.value), 
                         c("time", wilcox.test(time ~ dx, data = difference_table_treatment)$p.value))

colnames(pvalues_summary) <- c("test_values", "uncorrected_p_value")

write.csv(pvalues_summary, "results/tables/time_pvalues.csv", row.names = F)
write.csv(difference_table_treatment, "results/tables/time_datatable.csv", row.names = F)

