### Investigate differences between initial and follow up for thetayc and Fit
### Compare whether changes are different between the two groups
## Marc Sze

# Load needed functions and libraries
source('code/functions.R')

loadLibs(c("dplyr", "aod"))

# Load needed data
thetaCompTotal <- read.dist('data/process/final.thetayc.0.03.lt.ave.dist')

good_metaf <- read.csv("data/process/mod_metadata/good_metaF_final.csv", stringsAsFactors = F, header = T)

# Crate distance table with only initial and follow ups
difference_table_treatment <- pickDistanceValues(
  as.character(good_metaf$initial), as.character(good_metaf$followUp), 
  thetaCompTotal, good_metaf, c("Dx_Bin", "dx", "chemo_received", "radiation_received", "Surgery")) %>% 
  rename(chemo = chemo_received, rads = radiation_received, surgery = Surgery)

# Create table to store mean, sd, and pvalue from test
change_theta_fit_summary <- as.data.frame(matrix(nrow = 2, ncol = 12, 
  dimnames = list(rows = c("thetayc", "fit"), 
    cols = c(
      "adn_mean", "adn_SD", "adn_n", "srn_mean", "srn_SD", "srn_n", 
      "crc_mean", "crc_SD", "crc_n", "Pvalue_adnVsrn", "Pvalue_adnVcrc", "Pvalue_srnVcrc"))))

# Add the thetayc row of table
change_theta_fit_summary['thetayc', ] <- c(
  #Adds the summary stats for adenoma
  filter(difference_table_treatment, Dx_Bin == "adenoma") %>% 
    summarise(adn_mean = mean(distance), adn_SD = sd(distance), adn_n = length(distance)),  
  #Adds the summary stats for SRN
  filter(difference_table_treatment, Dx_Bin == "adv_adenoma") %>% 
    summarise(srn_mean = mean(distance), srn_SD = sd(distance), srn_n = length(distance)), 
  #Adds the summary stats for CRC
  filter(difference_table_treatment, Dx_Bin == "cancer") %>% 
    summarise(crc_mean = mean(distance), crc_SD = sd(distance), crc_n = length(distance)), 
  #Adds the pvalue comparisons
  Pvalue_adnVsrn = wilcox.test(
    filter(difference_table_treatment, Dx_Bin == "adenoma")[, "distance"], 
    filter(difference_table_treatment, Dx_Bin == "adv_adenoma")[, "distance"])$p.value, 
  Pvalue_adnVcrc = wilcox.test(
    filter(difference_table_treatment, Dx_Bin == "adenoma")[, "distance"], 
    filter(difference_table_treatment, Dx_Bin == "cancer")[, "distance"])$p.value, 
  Pvalue_srnVcrc = wilcox.test(
    filter(difference_table_treatment, Dx_Bin == "adv_adenoma")[, "distance"], 
    filter(difference_table_treatment, Dx_Bin == "cancer")[, "distance"])$p.value) 


# Save table for later use
write.csv(change_theta_fit_summary, 
  "data/process/tables/change_theta_fit_summary.csv")

#Write out difference table
write.csv(difference_table_treatment, "data/process/tables/difference_table.csv", row.names = F)

