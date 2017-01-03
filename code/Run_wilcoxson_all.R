### All OTU Paired Wilcoxson Test
### Are there differences in any OTU in initial and follow up
### Uses the full lesion model data
## Marc Sze


#Load needed libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "scales", "wesanderson", "ggplot2"))

#Load needed data sets
lesion_data <- read.csv("results/tables/full_test_data.csv", header = T, row.names = 1)
shared <- read.delim('data/process/final.0.03.subsample.shared', 
                     header=T, sep='\t') %>% select(Group, contains("Otu0"))

metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", header = T, stringsAsFactors = F)

#set row names for shared file
rownames(shared) <- shared$Group

#get OTUs to keep
OTUs_to_keep <- colnames(select(lesion_data, contains("Otu0")))

#get initial and follow up sample IDs (ALL, Adenoma Only, and CRC Only)
ALL_ids_to_get <- c(metaf$initial, metaf$followUp)
adn_ids_to_get <- c(filter(metaf, Diagnosis == "adenoma")[, "initial"], 
                    filter(metaf, Diagnosis == "adenoma")[, "followUp"])

crc_ids_to_get <- c(filter(metaf, Diagnosis != "adenoma")[, "initial"], 
                    filter(metaf, Diagnosis != "adenoma")[, "followUp"])

#subset the original shared file into an amenable testing formats
shared_to_test_ALL <- shared[as.character(ALL_ids_to_get), OTUs_to_keep]










