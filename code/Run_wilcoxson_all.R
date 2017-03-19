### All OTU Paired Wilcoxson Test
### Are there differences in any OTU in initial and follow up
### Uses the full lesion model data
## Marc Sze


#Load needed libraries
source('code/functions.R')

loadLibs(c("dplyr", "scales", "wesanderson", "caret"))

# Read in necessary data frames

shared <- read.delim('data/process/final.0.03.subsample.shared', 
                     header=T, sep='\t') %>% select(Group, contains("Otu0"))

metaI <- read.csv("data/process/mod_metadata/metaI_final.csv", 
                  stringsAsFactors = F, header = T) 

good_metaf <- read.csv('data/process/mod_metadata/good_metaf_final.csv', 
                       header = T, stringsAsFactors = F) %>% select(initial)

metaI <- filter(metaI, !(sample %in% good_metaf$initial))

#################################################################################
#                                                                               #
#                                                                               #
#               Data Clean up for better and faster modeling                    #
#                                                                               #
#################################################################################

# Remove follow up samples and join metadata with microbiome data
full_data <- inner_join(metaI, shared, by = c("sample" = "Group"))

# Filter and use only specific data

full_data <- full_data %>% 
  select(sample, contains("Otu0"))

#Filter out rows that are not all complete

full_data <- full_data[complete.cases(full_data), ]

#Reduce MetaI to match the full_data
metaI <- filter(
  metaI, as.character(sample) %in% as.character(full_data$sample))

rownames(full_data) <- full_data$sample
full_data <- select(full_data, -sample)

# Remove those with near zero variance 
# Similar to removal of singletons or removal based on percentage in sample
# Using near zero variance should cover the above two so could use only this instead

nzv <- nearZeroVar(full_data)
# By default, nearZeroVar will return the positions of the variables that are flagged to be problematic
# using the saveMetrics argument will output a table with more in-depth info
# Group stays in this case (because it is numberic) 
# but maybe better in future to transfer this to row names if not a numberic label

full_data <- full_data[, -nzv]

# Add the lesion variable to the full_data
full_data <- cbind(lesion = factor(metaI$lesion, 
                                   levels = c(0, 1), labels = c("No", "Yes")), full_data)

#Load needed data sets
metaf <- read.csv("data/process/mod_metadata/good_metaf_final.csv", header = T, stringsAsFactors = F)

#set row names for shared file
rownames(shared) <- shared$Group

#get OTUs to keep
OTUs_to_keep <- colnames(select(full_data, contains("Otu0")))

#get initial and follow up sample IDs (Adenoma, SRN, and CRC Only)
adn_ids_to_get <- c(filter(metaf, Dx_Bin == "adenoma")[, "initial"], 
                    filter(metaf, Dx_Bin == "adenoma")[, "followUp"])

srn_ids_to_get <- c(filter(metaf, Dx_Bin == "adv_adenoma")[, "initial"], 
                    filter(metaf, Dx_Bin == "adv_adenoma")[, "followUp"])


crc_ids_to_get <- c(filter(metaf, Dx_Bin == "cancer")[, "initial"], 
                    filter(metaf, Dx_Bin == "cancer")[, "followUp"])


#subset the original shared file into an amenable testing formats
shared_to_test <- list(adn = shared[as.character(adn_ids_to_get), OTUs_to_keep], 
                       srn = shared[as.character(srn_ids_to_get), OTUs_to_keep], 
                       crc = shared[as.character(crc_ids_to_get), OTUs_to_keep])

#Run and get pvalues for each of the data sets and multiple comparison correct
pvalue_results <- list(adn = NA, srn = NA, crc = NA)
bh_corrected <- list(adn = NA, srn = NA, crc = NA)

for(i in 1:length(pvalue_results)){
  
  pvalue_results[[i]] <- apply(shared_to_test[[i]], 2, function(x) 
    wilcox.test(x[1:(length(x)/2)], x[(length(x)/2+1):length(x)], paired = TRUE)$p.value)
  
  pvalue_results[[i]] <- pvalue_results[[i]][order(pvalue_results[[i]])]
  
  bh_corrected[[i]] <- p.adjust(pvalue_results[[i]], method = "BH")
}


#Create data table with pvalue results
final_data <- rbind(
  cbind(pvalue_results[["adn"]], 
        bh_corrected[["adn"]], 
        analysis = rep("adn", length(pvalue_results[["adn"]]))), 
  cbind(pvalue_results[["srn"]], 
        bh_corrected[["srn"]], 
        analysis = rep("srn", length(pvalue_results[["adn"]]))), 
  cbind(pvalue_results[["crc"]], 
        bh_corrected[["crc"]], 
        analysis = rep("crc", length(pvalue_results[["crc"]]))))

colnames(final_data) <- c("p_value", "BH_corrected", "analysis")

#Write out data table for future use
write.csv(final_data, "data/process/tables/OTU_paired_wilcoxson_test.csv", row.names = F)








