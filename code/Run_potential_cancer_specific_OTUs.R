### Pull out 4 specific OTUs of Interest
### Final figure on Potential Difference between Adenoma and Carcinoma
### P. micra, P. stomatis, P. assacharolytica, and F. nucleatum
## Marc Sze

# Load needed functions
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr"))

# Load needed data tables
tax <- read.delim('data/process/final.taxonomy', sep='\t', 
                  header=T, row.names=1)

# Convert taxa table to a data frame with columns for each taxonomic division and shared OTU assignment
tax_df <- data.frame(do.call('rbind', 
                             strsplit(as.character(tax$Taxonomy), ';'))) %>% 
  select(Domain = X1, Phyla = X2, Order = X3, Class = X4, Family = X5, Genus = X6) %>% 
  mutate(otu = rownames(tax))

# Remove the (100) from the columns and remove the unused tax data table
tax_df <- as.data.frame(apply(tax_df, 2, 
                              function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)
# Load in metad data, create Disease_Free variable, and remove values with NA
good_metaf <- read.csv(
  "data/process/mod_metadata/good_metaf_final.csv", 
  stringsAsFactors = F, header = T) %>% 
  mutate(lesion_follow = ifelse(Disease_Free == "n", 1, 0))

# Gather by lowest level classification (genus)
genera_data <- as.data.frame(get_tax_level_shared('data/process/final.shared', 'data/process/final.taxonomy', 6))

# Generate relative abundance
total_seqs <- rowSums(genera_data)

genera_data <- cbind(Group = rownames(genera_data), 
              as.data.frame(apply(genera_data, 2, 
                                  function(x) x/total_seqs)))

# Aggregate data and reshape for graphing
temp_data <- genera_data %>% filter(Group %in% as.character(c(good_metaf$initial, good_metaf$followUp))) %>% 
  slice(match(as.character(c(good_metaf$initial, good_metaf$followUp)), Group)) %>% 
  mutate(sampleType = c(rep("initial", length(good_metaf$initial)), 
                        rep("followup", length(good_metaf$followUp))), 
         Dx_Bin = rep(good_metaf$Dx_Bin, 2))

crc_select_data <- cbind(
  EDRN = rep(good_metaf$EDRN, 2), 
  Group = temp_data$Group, 
  sampleType = temp_data$sampleType, 
  Dx_Bin = temp_data$Dx_Bin, 
  Disease_free = rep(good_metaf$Disease_Free, 2), 
  select(temp_data, Fusobacterium, Parvimonas, Peptostreptococcus, Porphyromonas)) %>% 
  gather(key = Genus, 
         value = rel.abund, Fusobacterium, Parvimonas, Peptostreptococcus, Porphyromonas)


# Write out table for future use
write.csv(crc_select_data, "data/process/tables/adn_crc_maybe_diff.csv", row.names = F)

# Run statistics testing and BH correction
test_data <- temp_data %>% select(Fusobacterium, Parvimonas, Peptostreptococcus, Porphyromonas) %>% 
  mutate(Dx_Bin = rep(good_metaf$Dx_Bin, 2), 
         sampleType = temp_data$sampleType, 
         Disease_free = rep(good_metaf$Disease_Free, 2))

good_counts_init <- c("Fusobacterium", "Parvimonas", "Peptostreptococcus", "Porphyromonas")

pvalue_summary <- cbind(
  lesion_pvalue = apply(select(test_data, one_of(good_counts_init)), 2, 
                        function(x){
                          wilcox.test(x[test_data$sampleType == "initial"], 
                                      x[test_data$sampleType == "followup"], 
                                      paired = TRUE)$p.value}),
  
  crc_pvalue = apply(select(test_data, one_of(good_counts_init)), 2, 
                     function(x){
                       wilcox.test(x[test_data$sampleType == "initial" & test_data$Dx_Bin == "cancer"], 
                                   x[test_data$sampleType == "followup" & test_data$Dx_Bin == "cancer"], 
                                   paired = TRUE)$p.value}), 
  adn_pvalue = apply(select(test_data, one_of(good_counts_init)), 2, 
                     function(x){
                       wilcox.test(x[test_data$sampleType == "initial" & test_data$Dx_Bin != "cancer"], 
                                   x[test_data$sampleType == "followup" & test_data$Dx_Bin != "cancer"], 
                                   paired = TRUE)$p.value}))

adjusted_pvalues <- p.adjust(c(pvalue_summary[, "lesion_pvalue"], 
                               pvalue_summary[, "crc_pvalue"], 
                               pvalue_summary[, "adn_pvalue"]), method = "BH")

pvalue_summary <- cbind(
  pvalue_summary, lesion_BH = adjusted_pvalues[1:length(good_counts_init)], 
  crc_BH = adjusted_pvalues[(length(good_counts_init) + 1):(length(good_counts_init)*2)], 
  adn_BH = adjusted_pvalues[(length(good_counts_init)*2 + 1):(length(good_counts_init)*3)])

rownames(pvalue_summary) <- c("porp", "fn", "parv", "pept")

# Write out the pvalue table for future use
write.csv(pvalue_summary, "data/process/tables/adn_crc_maybe_pvalue_summary.csv")


