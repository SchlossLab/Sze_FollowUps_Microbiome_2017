### Create Figure 1
### Investigate differences between initial and follow up for thetayc and Fit

# Load needed functions and libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "vegan", "knitr"))

# Load needed data
thetaCompTotal <- dissplit('data/process/final.thetayc.0.03.lt.ave.dist',
  split=F, meta = F)
metaF <- read.csv("results/tables/mod_metadata/metaF_final.csv", 
  stringsAsFactors = F, header = T)

# Crate distance table with only initial and follow ups
difference_table_treatment <- pickDistanceValues(
  as.character(metaF$initial), as.character(metaF$followUp), 
  thetaCompTotal, metaF, c("fit_result", "fit_followUp", "Dx_Bin", "dx")) %>% 
  mutate(fit_difference = fit_result - fit_followUp)

# Get the pvalue for thetayc and fit difference between initial and follow up
pValueList <- list(thetaDiff = wilcox.test(distance ~ dx, 
  data = difference_table_treatment)$p.value, 
                   fitDiff = wilcox.test(fit_difference ~ dx, 
                    data = difference_table_treatment)$p.value)

# Create table to store mean, sd, and pvalue from test
change_theta_fit_summary <- as.data.frame(matrix(nrow = 2, ncol = 7, 
  dimnames = list(rows = c("thetayc", "fit"), 
    cols = c(
      "adn_mean", "adn_SD", "adn_n", "crc_mean", "crc_SD", "crc_n", "Pvalue"))))

# Add the thetayc row of table
change_theta_fit_summary['thetayc', ] <- 
c(mean(filter(difference_table_treatment, dx == "adenoma")[, "distance"]), 
    sd(filter(difference_table_treatment, dx == "adenoma")[, "distance"]), 
    length(filter(difference_table_treatment, dx == "adenoma")[, "distance"]), 
    mean(filter(difference_table_treatment, dx != "adenoma")[, "distance"]), 
    sd(filter(difference_table_treatment, dx != "adenoma")[, "distance"]), 
    length(filter(difference_table_treatment, dx != "adenoma")[, "distance"]), 
    pValueList$thetaDiff)

# Add the Fit row of table
change_theta_fit_summary['fit', ] <- c(
  -mean(filter(difference_table_treatment, dx == "adenoma")[, "fit_difference"], 
      na.rm = TRUE), 
   sd(filter(difference_table_treatment, dx == "adenoma")[, "fit_difference"], 
      na.rm = TRUE), 
   length(
    filter(difference_table_treatment, dx == "adenoma")[, "fit_difference"]), 
   -mean(
    filter(difference_table_treatment, dx != "adenoma")[, "fit_difference"]), 
   sd(
    filter(difference_table_treatment, dx != "adenoma")[, "fit_difference"]), 
   length(
    filter(difference_table_treatment, dx != "adenoma")[, "fit_difference"]), 
   pValueList$fitDiff)

# Save table for later use
write.csv(change_theta_fit_summary, 
  "results/tables/change_theta_fit_summary.csv")


#Difference between initial and follow up for thetayc and fit broken down by adenoma and cancer
diff_adn_v_crc <- grid.arrange(
  # Difference from bacterial community structure between adenoma and cancer
  ggplot(difference_table_treatment, 
    aes(factor(dx, levels = c("adenoma", "cancer")), distance, group = 1)) + 
    geom_jitter(aes(color=Dx_Bin), width = 0.3, size = 4, alpha = 0.5) + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + 
    scale_color_manual(name = "Lesion Type", 
      values = wes_palette("GrandBudapest"), 
                       breaks = c("Adenoma", "adv Adenoma", "Cancer"), 
                       labels = c("Adenoma", "SRN", "Cancer")) + 
    coord_cartesian(ylim = c(0, 1)) + ylab("Thetayc Distance") + 
    xlab("") + theme_bw() + ggtitle("A") + 
    theme(axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
          title = element_text(face="bold", hjust = 0)) + 
    annotate("text", label = paste("P-value = ", 
      round(pValueList[["thetaDiff"]], digits = 2)), x = 1.5, y = 0.7), 
  
  # Difference from fit between adenoma and cancer
  ggplot(difference_table_treatment, 
    aes(factor(dx, levels = c("adenoma", "cancer")), 
      fit_difference*-1, group = 1)) + 
    geom_jitter(aes(color=Dx_Bin), width = 0.3, size = 4, alpha = 0.5) + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + 
    scale_color_manual(name = "Lesion Type", 
      values = wes_palette("GrandBudapest"), 
                       breaks = c("Adenoma", "adv Adenoma", "Cancer"), 
                       labels = c("Adenoma", "SRN", "Cancer")) + 
    ylab("Change in Fit from Follow up to Initial") + 
    xlab("") + theme_bw() + ggtitle("B") + 
    theme(axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
          title = element_text(face="bold", hjust = 0)) + 
    annotate("text", label = paste("P-value = ", 
      format(round(pValueList[["fitDiff"]], digits = 7))), x = 1.5, y = 100)
)

ggsave(file = "results/figures/Figure1.pdf", diff_adn_v_crc, 
  width=8, height = 8, dpi = 300)
