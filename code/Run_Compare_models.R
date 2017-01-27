### Try to identify the important OTUs part 2
### Specifically compare initial follow up model to lesion model
### find similar OTUs and compare outcomes with multiple comparison testing
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr","reshape2", "scales"))

# Read in data tables
rf_otu_tax <- read.csv("results/tables/rf_otu_tax.csv", 
                       header = T, stringsAsFactors = F) %>% rename(otu = X)
if_rf_otu_tax <- read.csv("results/tables/if_rf_otu_tax.csv", 
                          header = T, stringsAsFactors = F) %>% rename(otu = X)

shared <- read.delim('data/process/final.0.03.subsample.shared', 
  header=T, sep='\t') %>% select(-label, -numOtus)

good_metaf <- read.csv(
  "results/tables/mod_metadata/good_metaf_final.csv", 
  stringsAsFactors = F, header = T) %>% 
  mutate(lesion_follow = ifelse(Disease_Free == "n", 1, 0)) %>% 
  filter(!is.na(fit_followUp))

# Make Relative Abundance
subsample <- rowSums(select(shared, -Group))
shared <- as.data.frame(apply(select(shared, -Group), 2, 
  function(x) x/subsample)) %>% mutate(Group = shared$Group)

# Find common OTUs
lesion_model_variables <- rf_otu_tax$otu
IF_model_variables <- if_rf_otu_tax$otu
common_variables <- lesion_model_variables[lesion_model_variables %in% IF_model_variables]

# Filter taxonomy for only common OTUs
common_taxa <- filter(if_rf_otu_tax, otu %in% common_variables)

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

write.csv(pvalue_table, "results/tables/pvalue_IF_lesion_common_imp_vars.csv")
