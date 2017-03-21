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

adn_model_imp_vars <- read.csv("data/process/tables/adn_rf_wCV_imp_vars_summary.csv", 
                           header = T, stringsAsFactors = F) %>% 
  filter(Variable != "fit_result")

crc_model_imp_vars <- read.csv("data/process/tables/crc_rf_wCV_imp_vars_summary.csv", 
                              header = T, stringsAsFactors = F) %>% 
  filter(Variable != "fit_result")

# Create taxonomy tables based on OTUs ID'd as important

adn_tax_model <- tax_df[adn_model_imp_vars$Variable, ]

crc_tax_model <- tax_df[crc_model_imp_vars$Variable, ]


# Write out the data
write.csv(adn_tax_model, "data/process/tables/adn_rf_otu_tax.csv")
write.csv(crc_tax_model, "data/process/tables/crc_rf_otu_tax.csv")






