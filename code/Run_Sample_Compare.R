### Select Sample Exploration 
### Investigate similarity and differences in individuals with increased positive probability on follow up.
### Focus specifically on carcinoma samples
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))

# Read in data tables
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", 
                       stringsAsFactors = F, header = T)
rf_prediction_summary <- read.csv(
  'results/tables/rf_prediction_summary.csv', header = T)
aucrf_model_cutoff <- read.csv('results/tables/aucrf_model_cutoffs.csv', 
                               header = T, row.names = 1)

shared <- read.delim('data/process/final.0.03.subsample.shared', header=T, sep='\t')

imp_vars <- read.csv('results/tables/rf_imp_vars_summary.csv', 
                     stringsAsFactors = F, header = T)

tax <- read.delim('data/process/final.taxonomy',
 sep='\t', stringsAsFactors = F, header=T, row.names=1)

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)
# This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components


# Create a data frame with only the higher than in follow up samples for wFit and woFit
increased_data <- list(wfit = c(), wofit = c())
models = c("wfit", "wofit")

for(i in 1:length(increased_data)){
  
  increased_data[[i]] <- filter(rf_prediction_summary, 
    model == paste(models[i]) & diagnosis != "adenoma")
  
  
  increased_data[[i]] <- increased_data[[i]][filter(
    increased_data[[i]], sample_type == "initial")[, "postive_probability"] < 
      filter(increased_data[[i]], sample_type == "followup")[, "postive_probability"], ]
  
  increased_data[[i]] <- filter(good_metaf, 
    EDRN %in% unique(increased_data[[i]][, "EDRN"])) %>% 
    select(initial, followUp) %>% melt() %>% 
    cbind(increased_data[[i]], .) %>% rename(sample = value) %>% 
    select(-variable)
  
}


### Try to select out the most relevant OTUs 

# Pare down the shared file and keep order of initial then follow up for wfit and wofit models
increased_shared<- list(wfit = c(), wofit = c())

for(i in 1:length(increased_shared)){
  
  increased_shared[[i]] <- 
    filter(shared, Group %in% rbind(good_metaf$initial, good_metaf$followUp)) %>% 
    select(Group, one_of(filter(imp_vars, imp_variable != "fit_result" & 
                                  model == paste(models[i]))[, 'imp_variable']))
  
  rownames(increased_shared[[i]]) <- increased_shared[[i]]$Group
  
  increased_shared[[i]] <- increased_shared[[i]][
    as.character(c(good_metaf$initial, good_metaf$followUp)), ]
}


# Get a delta change for each OTU between follow up and initial (follow - initial)

delta_shared <- list(wfit = c(), wofit = c())

for(i in 1:length(delta_shared)){
  
  delta_shared[[i]] <- select(
    increased_shared[[i]], -Group)[as.character(good_metaf$followUp), ] - 
    select(increased_shared[[i]], -Group)[as.character(good_metaf$initial), ]
  
  delta_shared[[i]] <- cbind(EDRN = good_metaf$EDRN, delta_shared[[i]])
}

rm(i, models, shared)

# Use a paired wilcoxson 

otu_pvalue_list <- list(wfit = c(), wofit = c())

for(i in 1:length(otu_pvalue_list)){
  
  otu_pvalue_list[[i]] <- sort(apply(select(delta_shared[[i]], -EDRN), 2, function(x) 
    wilcox.test(
      x[delta_shared[[i]]$EDRN %in% unique(increased_data[[i]][, "EDRN"])], 
      x[!(delta_shared[[i]]$EDRN %in% unique(increased_data[[i]][, "EDRN"]))])$p.value), 
    decreasing = FALSE)
}


# Correct for multiple comparisons
otu_bh_corr_list <- list(wfit = c(), wofit = c())

for(i in 1:length(otu_bh_corr_list)){
  
  otu_bh_corr_list[[i]] <- p.adjust(otu_pvalue_list[[i]], method = "BH")
}

    # None are below a 0.05 threshold
    # Use pvalue under 0.05 in uncorrected to see which OTUs might be important to this


# Create reference taxa list

imp_taxa_list <- list(wfit = createTaxaLabeller(
  tax_df[names(otu_pvalue_list[["wfit"]][otu_pvalue_list[[1]] < 0.05]), ]), 
  wofit = createTaxaLabeller(
    tax_df[names(otu_pvalue_list[["wofit"]][otu_pvalue_list[[1]] < 0.05]), ]))


# Create table with those under 0.05 threshold

wfit_in_probs_summary <- select(increased_shared[["wfit"]], 
  one_of(names(otu_pvalue_list[["wfit"]][otu_pvalue_list[["wfit"]] < 0.05]))) %>% 
  mutate(sampleType = c(rep("initial", length(rownames(good_metaf))), 
                        rep("followup", length(rownames(good_metaf))))) %>% 
  mutate(diagnosis = rep(good_metaf$Diagnosis, 2)) %>% 
  mutate(disease_free = rep(good_metaf$Disease_Free, 2)) %>% 
  mutate(EDRN = rep(good_metaf$EDRN, 2)) %>% 
  mutate(
    increased = ifelse(EDRN %in% unique(increased_data[[1]][, 'EDRN']), "Yes", "No")) %>% 
  melt(id = c("EDRN", "increased", "sampleType", "diagnosis", "disease_free")) %>% 
  mutate(taxa = imp_taxa_list[["wfit"]][variable])


wofit_in_probs_summary <- select(increased_shared[["wofit"]], 
                                one_of(names(
                                  otu_pvalue_list[["wofit"]][otu_pvalue_list[["wofit"]] < 0.05]))) %>% 
  mutate(sampleType = c(rep("initial", length(rownames(good_metaf))), 
                        rep("followup", length(rownames(good_metaf))))) %>% 
  mutate(diagnosis = rep(good_metaf$Diagnosis, 2)) %>% 
  mutate(disease_free = rep(good_metaf$Disease_Free, 2)) %>% 
  mutate(EDRN = rep(good_metaf$EDRN, 2)) %>% 
  mutate(
    increased = ifelse(EDRN %in% unique(increased_data[[1]][, 'EDRN']), "Yes", "No")) %>% 
  melt(id = c("EDRN", "increased", "sampleType", "diagnosis", "disease_free")) %>% 
  mutate(taxa = imp_taxa_list[["wofit"]][variable])

# Create data table with pvalues


wfit_inc_otu_pvalue_summary <- cbind(wilcox_pvalue = otu_pvalue_list[["wfit"]], 
                              benjamini_hochberg = otu_bh_corr_list[["wfit"]])

wofit_inc_otu_pvalue_summary <- cbind(wilcox_pvalue = otu_pvalue_list[["wofit"]], 
                              benjamini_hochberg = otu_bh_corr_list[["wofit"]])


# Write out relevant tables to be used in manuscript
write.csv(wfit_in_probs_summary, 
  "results/tables/wfit_in_probs_summary.csv", row.names = F)
write.csv(wofit_in_probs_summary, 
  "results/tables/wofit_in_probs_summary.csv", row.names = F)
write.csv(wfit_inc_otu_pvalue_summary, 
  "results/tables/wfit_inc_otu_pvalue_summary.csv")
write.csv(wofit_inc_otu_pvalue_summary, 
  "results/tables/wofit_inc_otu_pvalue_summary.csv")







