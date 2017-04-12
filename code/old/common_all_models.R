# Rough comparison of common OTUs
### Most concerned with commonalities between models (treatment and normal)
### Marc Sze


source('code/functions.R')

loadLibs(c("tidyr", "dplyr", "scales", "knitr", "rmarkdown"))

adn_treat <- read.csv("data/process/tables/adn_treatment_imp_vars_summary.csv")
srn_treat <- read.csv("data/process/tables/srn_treatment_imp_vars_summary.csv")
crc_treat <- read.csv("data/process/tables/crc_treatment_imp_vars_summary.csv")
common_norm <- read.csv('data/process/tables/pvalue_adn_srn_crc_common_imp_vars.csv', 
                                                       header = T, row.names = 1, stringsAsFactors = F)



adn_otus <- as.character(adn_treat$Variable[1:round(length(adn_treat$Variable)*0.10)])
srn_otus <- as.character(srn_treat$Variable[1:round(length(srn_treat$Variable)*0.10)])
crc_otus <- as.character(crc_treat$Variable[1:round(length(crc_treat$Variable)*0.10)])


length(rownames(filter(common_norm, otu %in% adn_otus) %>% select(otu)))
length(rownames(filter(common_norm, otu %in% srn_otus) %>% select(otu)))
length(rownames(filter(common_norm, otu %in% crc_otus) %>% select(otu)))


