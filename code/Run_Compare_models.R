### Try to identify the important OTUs part 2
### Specifically compare initial follow up model to lesion model
### find similar OTUs and compare outcomes with multiple comparison testing
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", 
  "gridExtra", "scales", "wesanderson", "caret", "randomForest"))

load("exploratory/RF_model_Imp_OTU.RData")

imp_vars_list <- list()

for(i in 1:length(test_tune_list)){
  
  imp_vars_list[[paste("run_", i, sep = "")]] <- 
    varImp(test_tune_list[[paste("data_split", i, sep = "")]], 
      scale = FALSE)$importance %>% 
    mutate(Variable = rownames(.)) %>% arrange(desc(Overall))
}


# Calculate number of times an OTU is in the top 10% of overall importance

OTU_appearance_table <- as.data.frame(data_frame(
  Variable = imp_vars_list[["run_1"]]$Variable) %>% 
mutate(total_appearance = 0))

rownames(OTU_appearance_table) <- OTU_appearance_table$Variable

for(j in 1:length(imp_vars_list)){
  
  tempVars <- imp_vars_list[[j]][c(
    1:round(length(rownames(imp_vars_list[[j]]))*0.10)), ][, "Variable"]
  
  for(i in tempVars){
    
    OTU_appearance_table[i, "total_appearance"] <- 
    OTU_appearance_table[i, "total_appearance"] + 1
  }
}

OTU_appearance_table <- arrange(OTU_appearance_table, 
  desc(total_appearance))

# Keep Those over 50% of the total 100 runs of 80/20 splits
OTU_appearance_table <- filter(OTU_appearance_table, total_appearance > 50)

# Write out the important variables to a table
write.csv(OTU_appearance_table, 
  "results/tables/if_rf_wCV_imp_vars_summary.csv", row.names = F)

rm(list = setdiff(ls(), c("OTU_appearance_table")))

# Read in data tables
rf_imp_vars_summary <- read.csv(
  "results/tables/rf_wCV_imp_vars_summary.csv", 
  stringsAsFactors = F, header = T) %>% 
  filter(Variable != "fit_result" & Variable != "BMI" & 
    Variable != "Gender.mYes" & Variable != "Age")

shared <- read.delim('data/process/final.0.03.subsample.shared', 
  header=T, sep='\t') %>% select(-label, -numOtus)
rownames(shared) <- shared$Group

good_metaf <- read.csv(
  "results/tables/mod_metadata/good_metaf_final.csv", 
  stringsAsFactors = F, header = T) %>% 
  mutate(lesion_follow = ifelse(Disease_Free == "n", 1, 0)) %>% 
  filter(!is.na(fit_followUp))

mod_metaf <- read.csv("results/tables/follow_up_probability_summary.csv", 
  stringsAsFactors = F, header = T) %>% 
  mutate(sampleID = c(good_metaf$initial, good_metaf$followUp))

metaI <- read.csv("results/tables/mod_metadata/metaI_final.csv", 
  stringsAsFactors = F, header = T)

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

# Make Relative Abundance
subsample <- rowSums(select(shared, -Group))
shared <- as.data.frame(apply(select(shared, -Group), 2, 
  function(x) x/subsample))

# Find common OTUs
common_vars <- filter(OTU_appearance_table, 
  Variable %in% rf_imp_vars_summary$Variable)

# Filter taxonomy for only common OTUs
tax_df <- tax_df[common_vars[, "Variable"], ]

# Filter Shared for only follow up samples and commmon OTUs
train_test_shared <- shared[as.character(metaI$sample[metaI$dx == "normal"]), common_vars$Variable]

shared <- shared[as.character(mod_metaf$sampleID), 
common_vars$Variable] %>% 
mutate(Dx_Bin = mod_metaf$Dx_Bin, sampleType = mod_metaf$sampleType, 
  EDRN = rep(good_metaf$EDRN, 2))

write.csv(train_test_shared, "results/tables/train_initial_imp_vars.csv")
write.csv(shared, "results/tables/if_imp_vars_table.csv")

# Use a paired wilcoxson 

otu_pvalue <- as.data.frame(matrix(nrow = length(rownames(common_vars)), 
  ncol = 1, dimnames = list(
    nrow = common_vars[, "Variable"], ncol = "Pvalue")))

for(i in 1:length(rownames(common_vars))){
  
  temp_init <- filter(shared, sampleType == "initial") %>% 
                            select(contains(common_vars[i, "Variable"]))
  
  temp_follow <- filter(shared, sampleType == "followup") %>% 
                              select(contains(common_vars[i, "Variable"]))
  
  otu_pvalue[i, "Pvalue"] <- wilcox.test(
      temp_init[, common_vars[i, "Variable"]], 
      temp_follow[, common_vars[i, "Variable"]], 
      paired = TRUE)$p.value
}


# Create P-value table for later use
pvalue_table <- cbind(otu = rownames(tax_df), tax_df, 
  Pvalue = otu_pvalue, BH_corr = p.adjust(otu_pvalue$Pvalue, 
    method = "BH"))

write.csv(pvalue_table, "results/tables/pvalue_common_imp_vars.csv")

# Create Figure 5 graph of the one significant OTU

shared <- select(shared, contains(rownames(
  pvalue_table[pvalue_table$BH_corr < 0.05, ])), Dx_Bin, sampleType, EDRN)

init_follow_change <- ggplot(shared, aes(factor(sampleType, 
  levels = c("initial", "followup")), Otu000012, 
group = factor(EDRN))) + 
  geom_line(aes(
    color = factor(Dx_Bin, 
      levels = c("adenoma", "adv_adenoma", "cancer")))) + 
  geom_point(aes(
    color = factor(Dx_Bin, 
      levels = c("adenoma", "adv_adenoma", "cancer")))) + 
  theme_bw() + ylab("Relative Abundance") + xlab("") + 
  ggtitle(paste(
    rownames(pvalue_table[pvalue_table$BH_corr < 0.05, ]), " (", 
    tax_df[rownames(pvalue_table[pvalue_table$BH_corr < 0.05, ]), 
    "Genus"], ")", sep ="")) + 
  scale_colour_manual(
    name = "Lesion Type", 
    labels = c("Adenoma", "SRN", "Cancer"), 
    values = wes_palette("GrandBudapest")) + 
  scale_x_discrete(
    breaks = c("initial", "followup"), 
    labels = c("Initial", "Follow Up")) + 
  geom_hline(
    data = train_test_shared, aes(
      yintercept = mean(train_test_shared[, 
        rownames(pvalue_table[pvalue_table$BH_corr < 0.05, ])])), 
    linetype = 2) + 
  theme(
    legend.position=c(0.85,0.85), 
    plot.margin = unit(c(1, 1, 1, 1), "lines"), 
    plot.title = element_text(size=20, face="bold", hjust = 0.5), 
    legend.title = element_text(face="bold"))


ggsave(file = "results/figures/Figure5.pdf", init_follow_change, 
       width=6.5, height = 8, dpi = 300)
