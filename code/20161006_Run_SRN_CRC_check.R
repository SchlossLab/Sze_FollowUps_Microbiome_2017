### Assess if we gain anything from the 
### Base analysis without metadata or OTU groupings
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))


# Load in data tables needed
withFit_data <- read.csv("results/tables/withFIT.models.datatable.csv", header = T, stringsAsFactors = F)
woFit_data <- read.csv("results/tables/noFIT.models.datatable.csv", header = T, stringsAsFactors = F)
withFit_cutoffs <- read.csv("results/tables/withFIT.cutoffs.csv", header = T, stringsAsFactors = F)
woFit_cutoffs <- read.csv("results/tables/noFIT.cutoffs.csv", header = T, stringsAsFactors = F)

tax <- read.delim('data/process/followUps.final.an.unique_list.0.03.cons.taxonomy', 
                  sep='\t', header=T, row.names=1)
shared <- read.delim('data/process/followUps.final.an.unique_list.0.03.subsample.0.03.filter.shared', 
                     header=T, sep='\t')

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)
# This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components

# Generate COunt data for with and wo Fit
WF_modelCountList <- list(
  SRNModel = createCountTable(withFit_data, withFit_cutoffs$cutoff[1], 
                              "SRNlesion_rf_opt", c("initial", "followup"), c("Cancer", "adv Adenoma", "Adenoma")), 
  LesionModel = createCountTable(withFit_data, withFit_cutoffs$cutoff[3], 
                                 "lesion_rf_opt", c("initial", "followup"), c("Cancer", "adv Adenoma", "Adenoma"))
)


WoF_modelCountList <- list(
  SRNModel = createCountTable(woFit_data, woFit_cutoffs$cutoff[1], 
                              "SRNlesion_rf_opt", c("initial", "followup"), c("Cancer", "adv Adenoma", "Adenoma")), 
  LesionModel = createCountTable(woFit_data, woFit_cutoffs$cutoff[3], 
                                 "lesion_rf_opt", c("initial", "followup"), c("Cancer", "adv Adenoma", "Adenoma"))
)


groupList <- c("Cancer", "adv Adenoma", "Adenoma")
LM_WF_pvalues <- c()
LM_WoF_pvalues <- c()

for(i in 1:(length(groupList))){
  
  LM_WF_pvalues <- c(LM_WF_pvalues, runFisherTest(WF_modelCountList[["LesionModel"]], 
                                                  Group_inv = paste(groupList[i]), side = "greater")$p.value)
  
  LM_WoF_pvalues <- c(LM_WoF_pvalues, runFisherTest(WoF_modelCountList[["LesionModel"]], 
                                                   Group_inv = paste(groupList[i]), side = "greater")$p.value)
  
  
  
}

pvalue_table <- as.data.frame(cbind(c("Cancer", "advAdenoma", "Adenoma"), 
                      LM_WF_pvalues, WF_padjust = p.adjust(LM_WF_pvalues, method = "bonferroni"), 
                      LM_WoF_pvalues, WoF_padjust = p.adjust(LM_WoF_pvalues, method = "bonferroni")), stringsAsFactors = F)

colnames(pvalue_table)[1] <- "Disease"

write.csv(pvalue_table, "results/tables/lesion_WF_IF_pvalue_summary.csv", row.names = F)
write.csv(WF_modelCountList[["LesionModel"]], "results/tables/lesion_PosCount_WF_IF_summary.csv", row.names = F)
write.csv(WoF_modelCountList[["LesionModel"]], "results/tables/lesion_PosCount_WoF_IF_summary.csv", row.names = F)









