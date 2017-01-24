### Create Figure 1
### Investigate differences between initial and follow up for thetayc and Fit

# Load needed functions and libraries
source('code/functions.R')

loadLibs(c("dplyr", "ggplot2", "reshape2", 
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

# Create table to store mean, sd, and pvalue from test
change_theta_fit_summary <- as.data.frame(matrix(nrow = 2, ncol = 7, 
  dimnames = list(rows = c("thetayc", "fit"), 
    cols = c(
      "adn_mean", "adn_SD", "adn_n", "crc_mean", "crc_SD", "crc_n", "Pvalue"))))

# Add the thetayc row of table
change_theta_fit_summary['thetayc', ] <- c(
  #Adds the summary stats for adenoma
  filter(difference_table_treatment, dx == "adenoma") %>% 
    summarise(adn_mean = mean(distance), adn_SD = sd(distance), adn_n = length(distance)),  
  #Adds the summary stats for CRC
  filter(difference_table_treatment, dx != "adenoma") %>% 
    summarise(crc_mean = mean(distance), crc_SD = sd(distance), crc_n = length(distance)), 
  #Adds the pvalue comparison between adenoma and CRC
  Pvalue = wilcox.test(distance ~ dx, 
                       data = difference_table_treatment)$p.value) 

# Add the Fit row of table
change_theta_fit_summary['fit', ] <- c(
  #Adds the summary stats for adenoma
  filter(difference_table_treatment, dx == "adenoma") %>% 
    summarise(adn_mean = -mean(fit_difference, na.rm = TRUE), adn_SD = sd(fit_difference, na.rm = TRUE), 
              adn_n = length(fit_difference)), 
  #Adds the summary stats for CRC
  filter(difference_table_treatment, dx != "adenoma") %>% 
    summarise(adn_mean = -mean(fit_difference, na.rm = TRUE), adn_SD = sd(fit_difference, na.rm = TRUE), 
              adn_n = length(fit_difference)), 
  #Adds the pvalue comparison between adenoma and CRC
  Pvalue = wilcox.test(fit_difference ~ dx, 
                       data = difference_table_treatment)$p.value)

# Save table for later use
write.csv(change_theta_fit_summary, 
  "results/tables/change_theta_fit_summary.csv")

#Write out difference table
write.csv(difference_table_treatment, "results/tables/difference_table.csv", row.names = F)
