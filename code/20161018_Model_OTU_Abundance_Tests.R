### Test amount increase/decrease of each model OTU
### Does the overall amount increase/decrease even if the proportion positive does not
## Marc Sze


# Load in needed functions and libraries

source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "vegan", "knitr"))

# Read in needed data

good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)
shared <- read.delim('data/process/followUps.final.an.unique_list.0.03.subsample.0.03.filter.shared', 
                     header=T, sep='\t')

tax <- read.delim('data/process/followUps.final.an.unique_list.0.03.cons.taxonomy', 
                  sep='\t', header=T, row.names=1)


wfit_RFOpt_vars <- read.csv("results/tables/lesion_RFOpt_Imp_Vars.csv", stringsAsFactors = F, header = T) %>% 
  rename(vars = x)

wofit_RFOpt_vars <- read.csv("results/tables/lesion_RFOpt_NOFIT_Imp_Vars.csv", stringsAsFactors = F, header = T) %>% 
  rename(vars = x)

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)
# This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components

# Create needed initial and follow up data tables for with fit model and test statistically
initial_lesion_wfit <- inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
  select(lesion, one_of(filter(wfit_RFOpt_vars, vars != "fit_result")[, "vars"]))

followups_lesion_wfit <- inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
  select(lesionf, one_of(filter(wfit_RFOpt_vars, vars != "fit_result")[, "vars"])) %>% 
           rename(lesion = lesionf)

adjusted_pvalues_abund <- get_abund_pvalues(initial_lesion_wfit, followups_lesion_wfit)


# create needed labels for Boruta picked important variables for each model
lesion_wfit_taxa <- tax_df[filter(adjusted_pvalues_abund, adj_pvalues < 0.1)[, 'otus'], ]
lesion_wfit_labs <- createTaxaLabeller(lesion_wfit_taxa)


# Create needed initial and follow up data tables for without fit model and test statistically
initial_lesion_wofit <- inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
  select(lesion, one_of(filter(wofit_RFOpt_vars, vars != "fit_result")[, "vars"]))

followups_lesion_wofit <- inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
  select(lesionf, one_of(filter(wofit_RFOpt_vars, vars != "fit_result")[, "vars"])) %>% 
  rename(lesion = lesionf)

wofit_adjusted_pvalues_abund <- get_abund_pvalues(initial_lesion_wofit, followups_lesion_wofit)

# create needed labels for Boruta picked important variables for each model
lesion_wofit_taxa <- tax_df[filter(wofit_adjusted_pvalues_abund, adj_pvalues < 0.1)[, 'otus'], ]
lesion_wofit_labs <- createTaxaLabeller(lesion_wofit_taxa)






