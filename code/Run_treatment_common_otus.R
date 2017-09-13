### Create data on treatment affects on common OTUs from each model
### Could these treatments influence findings of surgical removal
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr"))

# Read in data tables
chemo_rad_stats_summary <- read.csv('data/process/tables/crc_probs_chemo_rad_pvalue_summary.csv', 
                              header = T, stringsAsFactors = F, row.names = 1)

good_metaf <- read.csv('data/process/mod_metadata/metaF_final.csv', header = T, stringsAsFactors = F)

common_vars_summary <- read.csv('data/process/tables/pvalue_adn_srn_crc_common_imp_vars.csv', 
                                header = T, row.names = 1, stringsAsFactors = F)

chemo_rad_summary <- read.csv('data/process/tables/crc_chemo_rad_summary.csv', 
                              header = T, stringsAsFactors = F, row.names = 1)

shared <- read.delim('data/process/final.shared', header = T, stringsAsFactors = F) %>% 
  select(-label, -numOtus) %>% mutate(Group = as.character(Group))

# Convert to relative abundance
total_seqs <- rowSums(select(shared, -Group))

shared <- cbind(Group = shared$Group, 
                as.data.frame(apply(select(shared, -Group), 2, 
                                    function(x) x/total_seqs)))

# Shrink shared file down to only the samples initial then follow ups
samples_to_keep <- as.character(c(good_metaf$initial, good_metaf$followUp))
shared <- filter(shared, Group %in% samples_to_keep) %>% 
  slice(match(samples_to_keep, Group))

# Select only OTUs that were common
crc_shared <- shared %>% select(one_of(common_vars_summary$Variable)) %>% 
  mutate(sampleType = c(rep("initial", length(good_metaf$initial)), rep("followups", length(good_metaf$followUp))), 
         Dx_Bin = rep(good_metaf$Dx_Bin, 2), 
         chemo = rep(good_metaf$chemo_received, 2), 
         rads = rep(good_metaf$radiation_received, 2)) %>% 
  filter(Dx_Bin == "cancer")

test_data <- ((filter(crc_shared, sampleType == "initial") %>% select(contains("Otu0"))) - 
  (filter(crc_shared, sampleType == "followups") %>% select(contains("Otu0")))) %>% 
  mutate(chemo = crc_shared$chemo[crc_shared$sampleType == "initial"], 
         rads = crc_shared$rads[crc_shared$sampleType == "initial"])

# Create lowest ID labels
common_lowest_IDs <- common_vars_summary$tax_ID
table_labels <- c()
for(i in 1:length(common_lowest_IDs)){
  
  table_labels <- c(table_labels, rep(common_lowest_IDs[i], length(good_metaf$EDRN)))
}


# Add a lowest taxa ID to test data table
graph_ids <- c()
otus <- unique(common_vars_summary$Variable)

for(i in 1:length(otus)){
  
  graph_ids <- c(graph_ids, rep(common_lowest_IDs[i], 
                                length(rownames(filter(good_metaf, Dx_Bin == "cancer")))))
}


# Create data table to be used for graphing
common_data <- test_data %>% gather(key = otu, value = change, one_of(common_vars_summary$Variable)) %>% 
  mutate(lowest_ID = graph_ids)


# Calculate P-values and write out to csv
pvalue_table <- as.data.frame(matrix(nrow = length(otus), 
                                     ncol = 4, dimnames = list(
                                       rown = otus, 
                                       coln = c("Pvalue_chemo", "bh_chemo", "Pvalue_rads", "bh_rads"))))

for(i in 1:length(otus)){
  
  pvalue_table[i, "Pvalue_chemo"] <- wilcox.test(filter(common_data, chemo == "yes" & otu == otus[i])[, "change"], 
              filter(common_data, chemo == "no" & otu == otus[i])[, "change"])$p.value
  
  pvalue_table[i, "Pvalue_rads"] <- wilcox.test(filter(common_data, rads == "yes" & otu == otus[i])[, "change"], 
                                                 filter(common_data, rads == "no" & otu == otus[i])[, "change"])$p.value
}

pvalue_table <- pvalue_table %>% mutate(bh_chemo = p.adjust(Pvalue_chemo, method = "BH"), 
                                        bh_rads = p.adjust(Pvalue_rads, method = "BH"), 
                                        lowest_ID = common_lowest_IDs, 
                                        otu = otus)

write.csv(pvalue_table, "data/process/tables/chemo_rads_treatment_pvalue_summary.csv", row.names = F)



# Calculate P-values for surgery within adenoma combined group

# Select only OTUs that were common
adn_shared <- shared %>% select(one_of(common_vars_summary$Variable)) %>% 
  mutate(sampleType = c(rep("initial", length(good_metaf$initial)), rep("followups", length(good_metaf$followUp))), 
         Dx_Bin = rep(good_metaf$Dx_Bin, 2), 
         Surgery = rep(good_metaf$Surgery, 2)) %>% 
  filter(Dx_Bin == "adenoma" | Dx_Bin == "adv_adenoma")


adn_test_data <- ((filter(adn_shared, sampleType == "initial") %>% select(contains("Otu0"))) - 
                (filter(adn_shared, sampleType == "followups") %>% select(contains("Otu0")))) %>% 
  mutate(surg = adn_shared$Surgery[adn_shared$sampleType == "initial"])

# Add a lowest taxa ID to test data table
adn_graph_ids <- c()

for(i in 1:length(otus)){
  
  adn_graph_ids <- c(adn_graph_ids, rep(common_lowest_IDs[i], 
                                length(rownames(filter(good_metaf, Dx_Bin != "cancer")))))
}


# Create data table to be used for graphing
adn_common_data <- adn_test_data %>% gather(key = otu, value = change, one_of(common_vars_summary$Variable)) %>% 
  mutate(lowest_ID = adn_graph_ids)


# Calculate P-values and write out to csv
adn_pvalue_table <- as.data.frame(matrix(nrow = length(otus), 
                                     ncol = 2, dimnames = list(
                                       rown = otus, 
                                       coln = c("Pvalue_surg", "bh_surg"))))

for(i in 1:length(otus)){
  
  adn_pvalue_table[i, "Pvalue_surg"] <- wilcox.test(filter(adn_common_data, surg == "Y" & 
                                                         otu == otus[i])[, "change"], 
                                                 filter(adn_common_data, surg == "N" & 
                                                          otu == otus[i])[, "change"])$p.value
}

adn_pvalue_table <- adn_pvalue_table %>% mutate(bh_surg = p.adjust(Pvalue_surg, method = "BH"), 
                                        lowest_ID = common_lowest_IDs, otu = otus)

write.csv(adn_pvalue_table, "data/process/tables/all_adn_surg_pvalue_summary.csv", row.names = F)

