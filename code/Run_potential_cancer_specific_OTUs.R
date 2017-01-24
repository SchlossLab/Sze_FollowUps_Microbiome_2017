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
# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', 
                             strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c(
  "Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, 
                              function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)

good_metaf <- read.csv(
  "results/tables/mod_metadata/good_metaf_final.csv", 
  stringsAsFactors = F, header = T) %>% 
  mutate(lesion_follow = ifelse(Disease_Free == "n", 1, 0)) %>% 
  filter(!is.na(fit_followUp))

shared <- read.delim('data/process/final.0.03.subsample.shared', 
                     header=T, sep='\t') %>% select(-label, -numOtus)

shared <- filter(shared, Group %in% as.character(c(good_metaf$initial, good_metaf$followUp))) %>% 
  slice(match(as.character(c(good_metaf$initial, good_metaf$followUp)), Group)) %>% 
  mutate(sampleType = c(rep("initial", length(good_metaf$initial)), 
                        rep("followup", length(good_metaf$followUp))))

# Generate relative abundance
total_seqs <- rowSums(select(shared, -Group, -sampleType))[1]

shared <- cbind(Group = shared$Group, 
                sampleType = shared$sampleType, 
                as.data.frame(apply(select(shared, -Group, -sampleType), 2, 
                                    function(x) x/total_seqs)))


# Identify which OTUs are part of the genera of intrest
IDs <- mutate(tax_df, otu = rownames(tax_df)) %>% 
  filter(Genus == "Fusobacterium" | Genus == "Parvimonas" | 
           Genus == "Peptostreptococcus" | Genus == "Porphyromonas")

# Isolate the specific OTUs of interest  
shared_imp_init <- select(shared, Group, sampleType, one_of(IDs[, "otu"])) %>% 
  filter(sampleType == "initial") %>% select(-Group, -sampleType)

# Get number of positive counts by column
total_counts_init <- colSums(shared_imp_init != 0)
good_counts_init <- names(total_counts_init[total_counts_init > 10])

# Create data table for graphing

crc_select_data <- select(shared, Group, sampleType, one_of(good_counts_init)) %>% 
  mutate(EDRN = rep(good_metaf$EDRN, length(unique(shared$sampleType))), 
         Disease_Free = rep(good_metaf$Disease_Free, length(unique(shared$sampleType))), 
         Dx_Bin = rep(good_metaf$Dx_Bin, length(unique(shared$sampleType)))) %>% 
  melt(id.vars = c("Group", "sampleType", "EDRN", "Disease_Free", "Dx_Bin"), variable.name = "otu") %>% 
  mutate(tax_id = c(rep(as.character(tax_df[good_counts_init[1], "Genus"]), length(good_metaf$EDRN)*2), 
                    rep(as.character(tax_df[good_counts_init[2], "Genus"]), length(good_metaf$EDRN)*2), 
                    rep(as.character(tax_df[good_counts_init[3], "Genus"]), length(good_metaf$EDRN)*2), 
                    rep(as.character(tax_df[good_counts_init[4], "Genus"]), length(good_metaf$EDRN)*2)))
  
# Write out table for future use

write.csv(crc_select_data, "results/tables/adn_crc_maybe_diff.csv", row.names = F)

# Run statistics testing
pvalue_summary <- matrix(nrow = 4, ncol = 4, dimnames = list(
  nrow = c("fn", "parv", "pept", "porp"), ncol = c("crc_pvalue", "crc_BH", "adn_pvalue", "adn_BH")))

for(i in 1:length(good_counts_init)){
  
  pvalue_summary[i, "crc_pvalue"] <- wilcox.test(
    filter(crc_select_data, otu == good_counts_init[i], Dx_Bin == "cancer", 
           sampleType == "initial")[, "value"], 
    filter(crc_select_data, otu == good_counts_init[i], Dx_Bin == "cancer", 
           sampleType == "followup")[, "value"], paired = TRUE)$p.value
  
  pvalue_summary[i, "crc_BH"] <- p.adjust(pvalue_summary[i, "crc_pvalue"], method = "BH", n = 4)
  
}


for(i in 1:length(good_counts_init)){
  
  pvalue_summary[i, "adn_pvalue"] <- wilcox.test(
    filter(crc_select_data, otu == good_counts_init[i], Dx_Bin != "cancer", 
           sampleType == "initial")[, "value"], 
    filter(crc_select_data, otu == good_counts_init[i], Dx_Bin != "cancer", 
           sampleType == "followup")[, "value"], paired = TRUE)$p.value
  
  pvalue_summary[i, "adn_BH"] <- p.adjust(pvalue_summary[i, "adn_pvalue"], method = "BH", n = 4)
  
}

# Write out the pvalue table for future use

write.csv(pvalue_summary, "results/tables/adn_crc_maybe_pvalue_summary.csv")








