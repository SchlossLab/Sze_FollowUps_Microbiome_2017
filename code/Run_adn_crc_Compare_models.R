### Try to identify the important OTUs part 2
### Specifically compare initial follow up model to lesion model
### find similar OTUs and compare outcomes with multiple comparison testing
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr","reshape2", "scales"))

# Read in data tables
adn_rf_otu_tax <- read.csv("data/process/tables/adn_rf_otu_tax.csv", 
                       header = T, stringsAsFactors = F) %>% rename(otu = X)

srn_rf_otu_tax <- read.csv("data/process/tables/srn_rf_otu_tax.csv", 
                           header = T, stringsAsFactors = F) %>% rename(otu = X)

crc_rf_otu_tax <- read.csv("data/process/tables/crc_rf_otu_tax.csv", 
                          header = T, stringsAsFactors = F) %>% rename(otu = X)

shared <- read.delim('data/process/final.0.03.subsample.shared', 
  header=T, sep='\t') %>% select(-label, -numOtus)

good_metaf <- read.csv(
  "data/process/mod_metadata/good_metaf_final.csv", 
  stringsAsFactors = F, header = T) %>% 
  mutate(lesion_follow = ifelse(Disease_Free == "n", 1, 0))

# Make Relative Abundance
subsample <- rowSums(select(shared, -Group))
shared <- as.data.frame(apply(select(shared, -Group), 2, 
  function(x) x/subsample)) %>% mutate(Group = shared$Group)

# Find common OTUs
adn_model_variables <- adn_rf_otu_tax$otu
srn_model_variables <- srn_rf_otu_tax$otu
crc_model_variables <- crc_rf_otu_tax$otu
adn_srn_common_vars <- adn_model_variables[adn_model_variables %in% srn_model_variables]
common_variables <- crc_model_variables[crc_model_variables %in% adn_srn_common_vars]

# Filter taxonomy for only common OTUs
common_taxa <- filter(adn_rf_otu_tax, otu %in% common_variables)

# Filter Shared for only follow up samples and commmon OTUs
test_data <- slice(shared, match(c(good_metaf$initial, good_metaf$followUp), Group)) %>% 
  select(one_of(as.character(common_variables))) %>% 
  mutate(lesion = c(good_metaf$lesion, good_metaf$lesion_follow), 
         EDRN = rep(good_metaf$EDRN, 2), 
         sampleType = c(rep("initial", length(rownames(good_metaf))), 
                        rep("followup", length(rownames(good_metaf)))))

# Use a paired wilcoxson 
otu_pvalue <- apply(select(test_data, -lesion, -EDRN, -sampleType), 2, 
      function(x){
        wilcox.test(
          x[test_data$sampleType == "initial"], x[test_data$sampleType == "followup"], 
          paired = TRUE)$p.value})

# Create P-value table for later use
pvalue_table <- cbind(
  otu = common_variables, 
  lowest_ID = gsub("_unclassified", "", createTaxaLabeller(common_taxa)), 
  Pvalue = otu_pvalue, 
  BH_corr = p.adjust(otu_pvalue, method = "BH"))

write.csv(pvalue_table, "data/process/tables/pvalue_adn_srn_crc_common_imp_vars.csv")
