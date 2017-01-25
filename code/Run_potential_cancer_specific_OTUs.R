### Pull out 4 specific OTUs of Interest
### Final figure on Potential Difference between Adenoma and Carcinoma
### P. micra, P. stomatis, P. assacharolytica, and F. nucleatum
## Marc Sze

# Load needed functions
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "reshape2"))

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
  "results/tables/mod_metadata/good_metaf_final.csv", 
  stringsAsFactors = F, header = T) %>% 
  mutate(lesion_follow = ifelse(Disease_Free == "n", 1, 0)) %>% 
  filter(!is.na(fit_followUp))

# Load in shared file and modify to have initial and follow up
shared <- read.delim('data/process/final.0.03.subsample.shared', 
                     header=T, sep='\t') %>% select(-label, -numOtus) %>% 
  filter(Group %in% as.character(c(good_metaf$initial, good_metaf$followUp))) %>% 
  slice(match(as.character(c(good_metaf$initial, good_metaf$followUp)), Group)) %>% 
  mutate(sampleType = c(rep("initial", length(good_metaf$initial)), 
                        rep("followup", length(good_metaf$followUp))), 
         Dx_Bin = rep(good_metaf$Dx_Bin, 2))

# Generate relative abundance
total_seqs <- rowSums(select(shared, -Group, -sampleType))[1]

shared <- cbind(Group = shared$Group, 
                sampleType = shared$sampleType, 
                Dx_Bin = shared$Dx_Bin, 
                as.data.frame(apply(select(shared, -Group, -sampleType, -Dx_Bin), 2, 
                                    function(x) x/total_seqs)))

# Identify which OTUs are part of the genera of intrest
IDs <- filter(tax_df, Genus == "Fusobacterium" | Genus == "Parvimonas" | 
           Genus == "Peptostreptococcus" | Genus == "Porphyromonas")

# Isolate the specific OTUs of interest  
shared_imp_init <- select(shared, Group, sampleType, Dx_Bin, one_of(as.character(IDs[, "otu"]))) %>% 
  filter(sampleType == "initial") %>% select(-Group, -sampleType, -Dx_Bin)

# Get number of positive counts by column
total_counts_init <- colSums(shared_imp_init != 0)
good_counts_init <- names(total_counts_init[total_counts_init > 10])

# Create data table for graphing
crc_select_data <- select(shared, Group, sampleType, one_of(good_counts_init)) %>% 
  mutate(EDRN = rep(good_metaf$EDRN, length(unique(shared$sampleType))), 
         Disease_Free = rep(good_metaf$Disease_Free, length(unique(shared$sampleType))), 
         Dx_Bin = rep(good_metaf$Dx_Bin, length(unique(shared$sampleType)))) %>% 
  melt(id.vars = c("Group", "sampleType", "EDRN", "Disease_Free", "Dx_Bin"), variable.name = "otu") %>% 
  mutate(tax_id = c(rep(as.character(tax_df[tax_df$otu == good_counts_init[1], "Genus"]), length(good_metaf$initial)*2), 
                    rep(as.character(tax_df[tax_df$otu == good_counts_init[2], "Genus"]), length(good_metaf$initial)*2), 
                    rep(as.character(tax_df[tax_df$otu == good_counts_init[3], "Genus"]), length(good_metaf$initial)*2), 
                    rep(as.character(tax_df[tax_df$otu == good_counts_init[4], "Genus"]), length(good_metaf$initial)*2)))
  
# Write out table for future use
write.csv(crc_select_data, "results/tables/adn_crc_maybe_diff.csv", row.names = F)

# Run statistics testing and BH correction
test_data <- select(shared, Group, sampleType, Dx_Bin, one_of(good_counts_init))

pvalue_summary <- cbind(
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

pvalue_summary <- cbind(pvalue_summary, 
  crc_BH = p.adjust(pvalue_summary[, "crc_pvalue"], method = "BH"),
  adn_BH = p.adjust(pvalue_summary[, "adn_pvalue"], method = "BH"))

rownames(pvalue_summary) <- c("porp", "fn", "parv", "pept")

# Write out the pvalue table for future use
write.csv(pvalue_summary, "results/tables/adn_crc_maybe_pvalue_summary.csv")








