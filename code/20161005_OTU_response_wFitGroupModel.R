### Investigate Imp OTU variables at initial and follow up with fit as categorical
### Do the OTUs change that is reflective of loss of tumor or polyp?
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))


# Read in data and remove unneeded tables and lists
load("exploratory/RFwFitGroup.RData")

good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)


initial <- inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
  select(SRNlesion, lesion, fit_init_positive, contains("Otu0")) %>% rename(fit_positive = fit_init_positive)

followups <- inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
  select(lesionf, fit_follow_positive, contains("Otu0")) %>% rename(SRNlesion = lesionf) %>% 
  rename(fit_positive = fit_follow_positive)

followups <- cbind(good_metaf$lesionf, followups)
colnames(followups)[1] <- c("lesion")

rm(corr_pvalue_ROC_table, metaF, metaI, selected_train, sens_specif_table, impfactor_Data_List, modelList, 
   orig_probs, orig_rf_opt, orig_RF_run, orig_roc, rocNameList, selected_probs, selected_rf_opt, selected_RF_run, 
   selected_roc, fit_group_train, variableList)

tax <- read.delim('data/process/followUps.final.an.unique_list.0.03.cons.taxonomy', 
                  sep='\t', header=T, row.names=1)

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)
# This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components

# create needed labels for Boruta picked important variables for each model

SRNlesion_selected_taxa <- tax_df[filter(confirmed_vars[["SRNlesion"]], otus != "fit_positive")[, 'otus'], ]
SRNlesion_selected_labs <- createTaxaLabeller(SRNlesion_selected_taxa)

lesion_selected_taxa <- tax_df[filter(confirmed_vars[["lesion"]], otus != "fit_positive")[, 'otus'], ]
lesion_selected_labs <- createTaxaLabeller(lesion_selected_taxa)

cancer_positive <- as.data.frame(cbind(rep(as.character(good_metaf$Disease_Free), 18), rep(good_metaf$EDRN, 18))) %>% 
  rename(Disease_Free = V1, EDRN = V2)

# Graph of the SRN lesion model

SRN_model_impf_graphs <- follow_Abundance_FitGroup_Graph(SRNlesion_selected_taxa, SRNlesion_selected_labs, 
                                                    cancer_positive, initial, followups, select(good_metaf, -lesionf), 
                                                    0.08, 67, "SRNlesion", "Disease_Free",  "Cancer")

SRNplot <- grid.arrange(SRN_model_impf_graphs[['adenoma_OTUs']], SRN_model_impf_graphs[['cancer_OTUs']], 
                        SRN_model_impf_graphs[['adenoma_fit']], SRN_model_impf_graphs[['cancer_fit']])


ggsave(file = "results/figures/OTUselect_response_SRNFitGroupModel.tiff", SRNplot, width=10, height = 20, dpi = 300)


# Graph of the lesion model

lesion_model_impf_graphs <- follow_Abundance_FitGroup_Graph(lesion_selected_taxa, lesion_selected_labs, 
                                                       cancer_positive, initial, followups, select(good_metaf, -lesionf), 
                                                       0.08, 67, "lesion")

Lesionplot <- grid.arrange(lesion_model_impf_graphs[['adenoma_OTUs']], lesion_model_impf_graphs[['cancer_OTUs']], 
                           lesion_model_impf_graphs[['adenoma_fit']], lesion_model_impf_graphs[['cancer_fit']])

ggsave(file = "results/figures/OTUselect_response_LesionFitGroupModel.tiff", SRNplot, width=10, height = 20, dpi = 300)











