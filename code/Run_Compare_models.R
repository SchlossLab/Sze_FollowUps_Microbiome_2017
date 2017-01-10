### Try to identify the important OTUs part 2
### Specifically compare initial follow up model to lesion model
### find similar OTUs and compare outcomes with multiple comparison testing
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", 
  "gridExtra", "scales", "wesanderson", "caret", "randomForest"))

# Read in data tables
rf_otu_tax <- read.csv("results/tables/rf_otu_tax.csv", header = T, row.names = 1)
if_rf_otu_tax <- read.csv("results/tables/if_rf_otu_tax.csv", header = T, row.names = 1)

shared <- read.delim('data/process/final.0.03.subsample.shared', 
  header=T, sep='\t') %>% select(-label, -numOtus)
rownames(shared) <- shared$Group

good_metaf <- read.csv(
  "results/tables/mod_metadata/good_metaf_final.csv", 
  stringsAsFactors = F, header = T) %>% 
  mutate(lesion_follow = ifelse(Disease_Free == "n", 1, 0)) %>% 
  filter(!is.na(fit_followUp))

# Make Relative Abundance
subsample <- rowSums(select(shared, -Group))
shared <- as.data.frame(apply(select(shared, -Group), 2, 
  function(x) x/subsample))

# Find common OTUs
lesion_model_variables <- rownames(rf_otu_tax)
IF_model_variables <- rownames(if_rf_otu_tax)
common_variables <- lesion_model_variables[lesion_model_variables %in% IF_model_variables]

# Filter taxonomy for only common OTUs
common_taxa <- if_rf_otu_tax[common_variables, ]

# Filter Shared for only follow up samples and commmon OTUs
test_data <- shared[as.character(c(good_metaf$initial, good_metaf$followUp)), common_variables] %>% 
  mutate(lesion = c(good_metaf$lesion, good_metaf$lesion_follow), 
         EDRN = rep(good_metaf$EDRN, 2), 
         sampleType = c(rep("initial", length(rownames(good_metaf))), 
                        rep("followup", length(rownames(good_metaf)))))


# Use a paired wilcoxson 

otu_pvalue <- as.data.frame(matrix(nrow = length(rownames(common_taxa)), 
  ncol = 1, dimnames = list(
    nrow = rownames(common_taxa), ncol = "Pvalue")))

for(i in 1:length(rownames(common_taxa))){
  
  temp_init <- filter(test_data, sampleType == "initial") %>% select(contains(common_variables[i]))
  
  temp_follow <- filter(test_data, sampleType == "followup") %>% select(contains(common_variables[i]))
  
  otu_pvalue[i, "Pvalue"] <- wilcox.test(
      temp_init[, common_variables[i]], temp_follow[, common_variables[i]], 
      paired = TRUE)$p.value
}


# Create P-value table for later use
pvalue_table <- cbind(otu = rownames(common_taxa), lowest_ID = createTaxaLabeller(common_taxa), 
  Pvalue = otu_pvalue, BH_corr = p.adjust(otu_pvalue$Pvalue, 
    method = "BH"))

write.csv(pvalue_table, "results/tables/pvalue_IF_lesion_common_imp_vars.csv")

