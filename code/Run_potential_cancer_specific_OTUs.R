### Pull out 4 specific OTUs of Interest
### Final figure on Potential Difference between Adenoma and Carcinoma
### P. micra, P. stomatis, P. assacharolytica, and F. nucleatum
## Marc Sze


# Load needed functions
source('code/functions.R')
source('code/graphFunctions.R')


# Load needed libraries
loadLibs(c("dplyr", "ggplot2", "reshape2", "gridExtra", "scales", "wesanderson", "knitr", "rmarkdown"))

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
rownames(shared) <- shared$Group

shared <- shared[as.character(c(good_metaf$initial, good_metaf$followUp)), ]

# Generate row sums
total_seqs <- rowSums(select(shared, -Group))[1]

shared <- apply(select(shared, -Group), 2, function(x) x/total_seqs)
shared <- as.data.frame(shared)

# Identify which OTUs are part of the genera of intrest

fusobacterium_IDs <- mutate(tax_df, otu = rownames(tax_df)) %>% filter(Genus == "Fusobacterium" )
parvimonas_IDs <- mutate(tax_df, otu = rownames(tax_df)) %>% filter(Genus == "Parvimonas" )
peptostrep_IDs <- mutate(tax_df, otu = rownames(tax_df)) %>% filter(Genus == "Peptostreptococcus" )
porpho_IDs <- mutate(tax_df, otu = rownames(tax_df)) %>% filter(Genus == "Porphyromonas" )

# Isolate the specific OTUs of interest  

shared_imp_init <- select(shared, one_of(c(fusobacterium_IDs[, "otu"], parvimonas_IDs[, "otu"], 
                                   peptostrep_IDs[, "otu"], porpho_IDs[, "otu"]))) %>% slice(1:66)

# Get number of positive counts by column
total_counts_init <- colSums(shared_imp_init != 0)
good_counts_init <- names(total_counts_init[total_counts_init > 10])

# Create data table for graphing

crc_select_data <- as.data.frame(cbind(value = c(shared[, good_counts_init[1]], 
                                                 shared[, good_counts_init[2]], 
                                                 shared[, good_counts_init[3]], 
                                                 shared[, good_counts_init[4]]))) %>% 
  mutate(otu = c(rep(good_counts_init[1], length(rownames(shared))), 
                 rep(good_counts_init[2], length(rownames(shared))), 
                 rep(good_counts_init[3], length(rownames(shared))), 
                 rep(good_counts_init[4], length(rownames(shared))))) %>% 
  mutate(tax_id = c(rep(as.character(tax_df[good_counts_init[1], "Genus"]), length(good_metaf$EDRN)*2), 
                    rep(as.character(tax_df[good_counts_init[2], "Genus"]), length(good_metaf$EDRN)*2), 
                    rep(as.character(tax_df[good_counts_init[3], "Genus"]), length(good_metaf$EDRN)*2), 
                    rep(as.character(tax_df[good_counts_init[4], "Genus"]), length(good_metaf$EDRN)*2))) %>% 
  mutate(EDRN = rep(good_metaf$EDRN, length(good_counts_init)*2)) %>% 
  mutate(Disease_Free = rep(good_metaf$Disease_Free, length(good_counts_init)*2)) %>% 
  mutate(Dx_Bin = rep(good_metaf$Dx_Bin, length(good_counts_init)*2)) %>% 
  mutate(sampleType = rep(c(rep("initial", length(good_metaf$EDRN)), 
                            rep("followup", length(good_metaf$EDRN))), length(good_counts_init)))

# Write out table for future use

write.csv(crc_select_data, "results/tables/adn_crc_maybe_diff.csv", row.names = F)






