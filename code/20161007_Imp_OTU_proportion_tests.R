### Investigate Imp OTU variables at initial and follow up (Fit and noFit model)
### Measure proportions originally present versus no longer there
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))

# Read in data and remove unneeded tables and lists
load("exploratory/RFwFit.RData")
rm(corr_pvalue_ROC_table, data, metaF, metaI, selected_train, sens_specif_table, impfactor_Data_List, modelList, 
   orig_probs, orig_rf_opt, orig_RF_run, orig_roc, rocNameList, selected_probs, selected_rf_opt, selected_RF_run, 
   selected_roc, train, variableList, confirmed_vars)

# Convert taxa table to a data frame with columns for each taxonomic division
tax <- read.delim('data/process/followUps.final.an.unique_list.0.03.cons.taxonomy', 
                  sep='\t', header=T, row.names=1)

tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)
# This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components


# Read in data and remove unneeded tables and lists
load("exploratory/RFwoFit.RData")
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)

initial <- inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
  select(SRNlesion, lesion, fit_result, contains("Otu0"))

followups <- inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
  select(lesionf, fit_followUp, contains("Otu0")) %>% rename(SRNlesion = lesionf) %>% rename(fit_result = fit_followUp)

followups <- cbind(good_metaf$lesionf, followups)
colnames(followups)[1] <- c("lesion")

rm(corr_pvalue_ROC_table, metaF, metaI, selected_train, sens_specif_table, impfactor_Data_List, modelList, 
   orig_probs, orig_rf_opt, orig_RF_run, orig_roc, rocNameList, selected_probs, selected_rf_opt, selected_RF_run, 
   selected_roc, train, variableList, confirmed_vars)


RFopt_vars_wFit <- read.csv("results/tables/lesion_RFOpt_Imp_Vars.csv", stringsAsFactors = F, header = T) %>% rename(otus = x)
RFopt_vars_woFit <- read.csv("results/tables/lesion_RFOpt_NOFIT_Imp_Vars.csv", stringsAsFactors = F, header = T) %>% rename(otus = x)


# Read in data tables to be used with prediction 
withfit_model_data <- read.csv("results/tables/withFIT.models.datatable.csv", header = T, stringsAsFactors = F)
wofit_model_data <- read.csv("results/tables/noFIT.models.datatable.csv", header = T, stringsAsFactors = F)

withfit_cutoffs <- read.csv("results/tables/withFIT.cutoffs.csv", header = T, stringsAsFactors = F)
wofit_cutoffs <- read.csv("results/tables/noFIT.cutoffs.csv", header = T, stringsAsFactors = F)

# create needed labels for Boruta picked important variables for each model

lesion_selected_taxa_WF <- tax_df[filter(RFopt_vars_wFit, otus != "fit_result")[, 'otus'], ]
lesion_selected_labs_WF <- createTaxaLabeller(lesion_selected_taxa_WF)

lesion_selected_taxa_WoF <- tax_df[filter(RFopt_vars_woFit)[, 'otus'], ]
lesion_selected_labs_WoF <- createTaxaLabeller(lesion_selected_taxa_WoF)


# pull metadata of those that did not respond

test <- filter(wofit_model_data, L1 == "lesion_rf_opt", diagnosis != "adenoma", time_point == "followup", 
               positive > wofit_cutoffs$cutoff[3])

# Count how many are positive for the four Imp OTUs previously associated with cancer

shared_select_WF <- select(shared, Group, one_of(rownames(lesion_selected_taxa_WF))) %>% 
  filter(as.character(Group) %in% as.character(good_metaf$initial)) %>% rbind(
    select(shared, Group, one_of(rownames(lesion_selected_taxa_WF))) %>% 
      filter(as.character(Group) %in% as.character(good_metaf$followUp)))

init_shared_select_WF <- inner_join(shared_select_WF, good_metaf, by = c("Group" = "initial")) %>% 
  select(Group, contains("Otu0"), Diagnosis, Disease_Free, Surgery, chemo_received, radiation_received) %>% 
  filter(Diagnosis != "adenoma")

follow_shared_select_WF <- inner_join(shared_select_WF, good_metaf, by = c("Group" = "followUp")) %>% 
  select(Group, contains("Otu0"), Diagnosis, Disease_Free, Surgery, chemo_received, radiation_received) %>% 
  filter(Diagnosis != "adenoma")


shared_select_WoF <- select(shared, Group, one_of(rownames(lesion_selected_taxa_WoF))) %>% 
  filter(as.character(Group) %in% as.character(good_metaf$initial)) %>% rbind(
    select(shared, Group, one_of(rownames(lesion_selected_taxa_WoF))) %>% 
      filter(as.character(Group) %in% as.character(good_metaf$followUp)))

init_shared_select_WoF <- inner_join(shared_select_WoF, good_metaf, by = c("Group" = "initial")) %>% 
  select(Group, contains("Otu0"), Diagnosis, Disease_Free, Surgery, chemo_received, radiation_received) %>% 
  filter(Diagnosis != "adenoma")

follow_shared_select_WoF <- inner_join(shared_select_WoF, good_metaf, by = c("Group" = "followUp")) %>% 
  select(Group, contains("Otu0"), Diagnosis, Disease_Free, Surgery, chemo_received, radiation_received) %>% 
  filter(Diagnosis != "adenoma")


# Run fisher exact test on the sample sets

OTU_WF_pvalues <- c()
for(i in 1:length(rownames(lesion_selected_taxa_WF))){
  
  tempData <- create_two_by_two(init_shared_select_WF, follow_shared_select_WF, OTU = rownames(lesion_selected_taxa_WF)[i])
  OTU_WF_pvalues <- c(OTU_WF_pvalues, fisher.test(tempData)$p.value)
  names(OTU_WF_pvalues)[i] <- rownames(lesion_selected_taxa_WF)[i]
  
}

OTU_WF_pvalues_summary <- as.data.frame(cbind(
  OTU_WF_pvalues, p.adjust(OTU_WF_pvalues, method = "bonferroni"), 
  lesion_selected_labs_WF, rownames(lesion_selected_taxa_WF)))

colnames(OTU_WF_pvalues_summary) <- c("OTU_WF_pvalues", "OTU_WF_adjust_pvalues", "Tax_ID", "otus") 


OTU_WoF_pvalues <- c()
for(i in 1:length(rownames(lesion_selected_taxa_WoF))){
  
  tempData <- create_two_by_two(init_shared_select_WoF, follow_shared_select_WoF, OTU = rownames(lesion_selected_taxa_WoF)[i])
  OTU_WoF_pvalues <- c(OTU_WoF_pvalues, fisher.test(tempData)$p.value)
  names(OTU_WoF_pvalues)[i] <- rownames(lesion_selected_taxa_WoF)[i]
}

OTU_WoF_pvalues_summary <- as.data.frame(cbind(
  OTU_WoF_pvalues, p.adjust(OTU_WoF_pvalues, method = "bonferroni"), 
  lesion_selected_labs_WoF, rownames(lesion_selected_taxa_WoF)))

colnames(OTU_WoF_pvalues_summary) <- c("OTU_WoF_pvalues", "OTU_WoF_adjust_pvalues", "Tax_ID", "otus") 


# figure out if having any of the OTUs in more represented in those with chemo/rads versus none

chemo_OTU_WF_pvalues <- c()
for(i in 1:length(rownames(lesion_selected_taxa_WF))){
  
  tempData <- advanced_two_by_two(init_shared_select_WF, follow_shared_select_WF, OTU = rownames(lesion_selected_taxa_WF)[i], 
                                Treatment = "chemo_received")
  chemo_OTU_WF_pvalues <- c(chemo_OTU_WF_pvalues, fisher.test(tempData)$p.value)
  names(chemo_OTU_WF_pvalues)[i] <- rownames(lesion_selected_taxa_WF)[i]
  
}

chemo_OTU_WF_pvalues_summary <- as.data.frame(cbind(
  chemo_OTU_WF_pvalues, p.adjust(chemo_OTU_WF_pvalues, method = "bonferroni"), 
  lesion_selected_labs_WF, rownames(lesion_selected_taxa_WF)))

colnames(chemo_OTU_WF_pvalues_summary) <- c("chemo_pvalues", "chemo_adjust_pvalues", "Tax_ID", "otus") 


chemo_OTU_WoF_pvalues <- c()
for(i in 1:length(rownames(lesion_selected_taxa_WoF))){
  
  tempData <- advanced_two_by_two(init_shared_select_WoF, follow_shared_select_WoF, 
                                  OTU = rownames(lesion_selected_taxa_WoF)[i], Treatment = "chemo_received")
  chemo_OTU_WoF_pvalues <- c(chemo_OTU_WoF_pvalues, fisher.test(tempData)$p.value)
  names(chemo_OTU_WoF_pvalues)[i] <- rownames(lesion_selected_taxa_WoF)[i]
  
}

chemo_OTU_WoF_pvalues_summary <- as.data.frame(cbind(
  chemo_OTU_WoF_pvalues, p.adjust(chemo_OTU_WoF_pvalues, method = "bonferroni"), 
  lesion_selected_labs_WoF, rownames(lesion_selected_taxa_WoF)))

colnames(chemo_OTU_WoF_pvalues_summary) <- c("chemo_pvalues", "chemo_adjust_pvalues", "Tax_ID", "otus") 


# Do same thing but for only the big 4 OTU126, OTU566, OTU397, OTU205

shared_select_crcSp <- select(shared, Group, one_of(c("Otu000126", "Otu000566", "Otu000397", "Otu000205"))) %>% 
  filter(as.character(Group) %in% as.character(good_metaf$initial)) %>% rbind(
    select(shared, Group, one_of(c("Otu000126", "Otu000566", "Otu000397", "Otu000205"))) %>% 
      filter(as.character(Group) %in% as.character(good_metaf$followUp)))

init_shared_select_crcSp <- inner_join(shared_select_crcSp, good_metaf, by = c("Group" = "initial")) %>% 
  select(Group, contains("Otu0"), Diagnosis, Disease_Free, Surgery, chemo_received, radiation_received) %>% 
  filter(Diagnosis != "adenoma")

follow_shared_select_crcSp <- inner_join(shared_select_crcSp, good_metaf, by = c("Group" = "followUp")) %>% 
  select(Group, contains("Otu0"), Diagnosis, Disease_Free, Surgery, chemo_received, radiation_received) %>% 
  filter(Diagnosis != "adenoma")


# Run fisher exact test on the sample sets

OTU_crcSp_pvalues <- c()
crcOTUs <- c("Otu000126", "Otu000566", "Otu000397", "Otu000205")

for(i in 1:length(crcOTUs)){
  
  tempData <- create_two_by_two(init_shared_select_crcSp, follow_shared_select_crcSp, OTU = crcOTUs[i])
  OTU_crcSp_pvalues <- c(OTU_crcSp_pvalues, fisher.test(tempData)$p.value)
  names(OTU_crcSp_pvalues)[i] <- rownames(crcOTUs)[i]
  
}

OTU_crcSp_pvalues_summary <- as.data.frame(cbind(
  OTU_crcSp_pvalues, p.adjust(OTU_crcSp_pvalues, method = "bonferroni"), as.character(tax_df[crcOTUs, "Genus"]), crcOTUs))

colnames(OTU_crcSp_pvalues_summary) <- c("pvalues", "adjust_pvalues", "Tax_ID", "otus") 


# figure out if having any of the OTUs in more represented in those with chemo/rads versus none

chemo_OTU_crcSp_pvalues <- c()
for(i in 1:length(crcOTUs)){
  
  tempData <- advanced_two_by_two(init_shared_select_crcSp, follow_shared_select_crcSp, 
                                  OTU = crcOTUs[i], Treatment = "chemo_received")
  chemo_OTU_crcSp_pvalues <- c(chemo_OTU_crcSp_pvalues, fisher.test(tempData)$p.value)
  names(chemo_OTU_crcSp_pvalues)[i] <- crcOTUs[i]
  
}

chemo_OTU_crcSp_pvalues_summary <- as.data.frame(cbind(
  chemo_OTU_crcSp_pvalues, p.adjust(chemo_OTU_crcSp_pvalues, method = "bonferroni"), 
  as.character(tax_df[crcOTUs, "Genus"]), crcOTUs))

colnames(chemo_OTU_crcSp_pvalues_summary) <- c("chemo_pvalues", "chemo_adjust_pvalues", "Tax_ID", "otus") 

# Write out important summary pvalue tables

write.csv(OTU_WF_pvalues_summary, "results/tables/lesion_OTU_WF_Pvalue_summary.csv", row.names = F)
write.csv(OTU_WoF_pvalues_summary, "results/tables/lesion_OTU_WoF_Pvalue_summary.csv", row.names = F)

write.csv(chemo_OTU_WF_pvalues_summary, "results/tables/chemo_lesion_OTU_WF_Pvalue_summary.csv", row.names = F)
write.csv(chemo_OTU_WoF_pvalues_summary, "results/tables/chemo_lesion_OTU_WoF_Pvalue_summary.csv", row.names = F)

write.csv(OTU_crcSp_pvalues_summary, "results/tables/lesion_OTU_crcSp_Pvalue_summary.csv", row.names = F)
write.csv(chemo_OTU_crcSp_pvalues_summary, "results/tables/chemo_lesion_OTU_crcSp_Pvalue_summary.csv", row.names = F)













