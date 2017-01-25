### Create specialized taxonomy for each of the OTU Tables
### TO be used in conjuction with original OTU pulled variables
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", 
           "wesanderson", "caret", "randomForest"))

#Read in needed data
tax <- read.delim('data/process/final.taxonomy', sep='\t', header=T, row.names=1)
# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)

if_model_imp_vars <- read.csv("results/tables/IF_rf_wCV_imp_vars_summary.csv", 
                           header = T, stringsAsFactors = F) %>% 
  filter(Variable != "fit_result")

model_imp_vars <- read.csv("results/tables/rf_wCV_imp_vars_summary.csv", 
                              header = T, stringsAsFactors = F) %>% 
  filter(Variable != "fit_result")

# Create taxonomy tables based on OTUs ID'd as important

tax_model <- tax_df[model_imp_vars$Variable, ]

if_tax_model <- tax_df[if_model_imp_vars$Variable, ]


# Write out the data
write.csv(tax_model, "results/tables/rf_otu_tax.csv")
write.csv(if_tax_model, "results/tables/if_rf_otu_tax.csv")






