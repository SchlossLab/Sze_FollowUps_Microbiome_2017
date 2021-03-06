### Check up on What Six samples that increase are doing
### Create data files to graph these differences
## Marc Sze

## Load in necessary libraries
source('code/functions.R')
loadLibs(c("dplyr", "tidyr"))


## Load needed data files
good_metaf <- read.csv('data/process/mod_metadata/good_metaf_final.csv', header = T, stringsAsFactors = F)
metaI <- read.csv('data/process/mod_metadata/metaI_final.csv', header = T, stringsAsFactors = F) %>% 
  filter(dx == "normal" | dx == "cancer")
crc_imp_vars <- read.csv('data/process/tables/crc_rf_otu_tax.csv', header = T, row.names = 1)
crc_probs <- read.csv('data/process/tables/crc_reduced_follow_up_probability_summary.csv', header = T)
shared <- read.delim('data/process/final.shared', header = T, stringsAsFactors = F) %>% select(-label, -numOtus)

## Get Rel. Abund
counts <- rowSums(select(shared, contains("Otu")))
otu.rel.abund <- shared
otu.rel.abund[, colnames(select(otu.rel.abund, contains("Otu")))] <- apply(
  select(otu.rel.abund, 
         contains("Otu")), 2, function(x) (x/counts)*100) 

## Get list of normal and cancer samples in training set
train_samples <- metaI$sample

## Get only the samples that increase on follow up
unique_edrn <- unique((crc_probs %>% 
  mutate(diffs = rep(filter(crc_probs, sampleType == "followup")[, 'Yes'] - filter(crc_probs, sampleType == "initial")[, 'Yes'], 2)) %>% 
  filter(diffs > 0) %>% select(EDRN))[, "EDRN"])

## shrink shared file down to needed samples and otus and combine into one data table
sample_ids <- unique(select(crc_probs, EDRN)[, "EDRN"])
crc_samples <- good_metaf %>% slice(match(sample_ids, EDRN))
crc_probs <- mutate(crc_probs, Group = c(crc_samples$initial, crc_samples$followUp))

combined_data <- otu.rel.abund %>% slice(match(c(good_metaf$initial, good_metaf$followUp), Group)) %>% 
  select(Group, one_of(rownames(crc_imp_vars))) %>% inner_join(crc_probs, by = "Group") %>% 
  mutate(probs_increase = ifelse(EDRN %in% unique_edrn, invisible("Yes"), invisible("No")))


test_data <- (filter(combined_data, sampleType == "followup") %>% select(contains("Otu")) - 
  filter(combined_data, sampleType == "initial") %>% select(contains("Otu"))) %>% 
  mutate(probs_increase = (select(combined_data, probs_increase) %>% slice(1:length(sample_ids)))[, "probs_increase"]) %>% 
  gather(key = otu, value = rel.abund, -probs_increase)


## shrink shared files down to only those used in the training
train_data <- shared %>% 
  slice(match(train_samples, Group)) %>% 
  select(Group, one_of(rownames(crc_imp_vars))) %>% 
  inner_join(metaI, by = c("Group" = "sample"))

summary_train <- train_data %>% group_by(dx) %>% select(dx, contains("Otu")) %>% 
  summarise_each(funs(median)) %>% gather(key = otu, value = median_value, -dx)

# log OTU medians higher in cancer and those higher in normal

test <- as.data.frame(matrix(nrow = length(rownames(crc_imp_vars)), 
                             ncol = 2, 
                             dimnames = list(rown = c(), coln = c("otu", "group")))) %>% 
  mutate(otu = rownames(crc_imp_vars))

for(i in 1:length(test$otu)){
  
  test[i, "group"] <- ifelse(
    filter(summary_train, otu == rownames(crc_imp_vars)[i] & dx == "normal")[, "median_value"] > 
      filter(summary_train, otu == rownames(crc_imp_vars)[i] & dx == "cancer")[, "median_value"], 
    invisible("normal"), invisible("cancer"))
  
  test[i, "group"] <- ifelse(
    filter(summary_train, otu == rownames(crc_imp_vars)[i] & dx == "normal")[, "median_value"] == 
      filter(summary_train, otu == rownames(crc_imp_vars)[i] & dx == "cancer")[, "median_value"], 
    invisible("neither"), invisible(test[i, "group"]))
}
  
# Store the otus for each specific group
crc_only_group <- filter(test, group == "cancer")[, "otu"]
normal_only_group <- filter(test, group == "normal")[, "otu"]

# Test if any increase or decrease change is significant
imp_otus <- rownames(crc_imp_vars)

change_pvalue_summary <- as.data.frame(
  matrix(nrow = length(imp_otus), 
         ncol = 7, 
         dimnames = list(rown = imp_otus, 
                         coln = c("otu", "pvalue", "bh", "median_Y", "median_N", "IQR_Y", "IQR_N")))) %>% 
  mutate(otu = imp_otus)

summary_test <- group_by(test_data, probs_increase, otu) %>% summarise_each(funs(median, IQR)) %>% 
  slice(match(change_pvalue_summary$otu, otu))

for(i in 1:length(imp_otus)){
  
  change_pvalue_summary[i, "pvalue"] <- wilcox.test(filter(test_data, otu == imp_otus[i] & probs_increase == "Yes")[, "rel.abund"], 
                                                    filter(test_data, otu == imp_otus[i] & probs_increase == "No")[, "rel.abund"])$p.value
  
  change_pvalue_summary[i, "median_Y"] <- filter(summary_test, probs_increase == "Yes" & otu == imp_otus[i])[, "median"]
  change_pvalue_summary[i, "median_N"] <- filter(summary_test, probs_increase == "No" & otu == imp_otus[i])[, "median"]
  change_pvalue_summary[i, "IQR_Y"] <- filter(summary_test, probs_increase == "Yes" & otu == imp_otus[i])[, "IQR"]
  change_pvalue_summary[i, "IQR_N"] <- filter(summary_test, probs_increase == "No" & otu == imp_otus[i])[, "IQR"]
  
}

change_pvalue_summary <- mutate(change_pvalue_summary, bh = p.adjust(pvalue, method = "BH"))

write.csv(change_pvalue_summary, "data/process/tables/inc_probs_crc_imp_otus_summary.csv", row.names = F)


## Check what the positivity is for the usual suspect OTUs for those with increased probs
usual_crc <- c("Otu000202", "Otu001273", "Otu000442")


## Add up Oral Pathogens and add up residents and compare
combined_residents <- tbl_df(combined_data$probs_increase) %>% rename(probs_increase = value) %>% 
  mutate(
    EDRN = combined_data$EDRN, 
    sampleType = combined_data$sampleType, 
    residents = (select(combined_data, one_of(imp_otus)) %>% select(-one_of(usual_crc)) %>% rowSums()), 
    oral_path = (select(combined_data, one_of(usual_crc)) %>% rowSums()), 
    crc_only = (select(combined_data, one_of(crc_only_group)) %>% rowSums()), 
    normal_only = (select(combined_data, one_of(normal_only_group)) %>% rowSums()))


write.csv(combined_residents, "data/process/tables/inc_probs_crc_oral_residents_data.csv", row.names = F)


## Test differences


test <- combined_residents %>% group_by(probs_increase, sampleType) %>% 
  select(-EDRN) %>% summarise_each(funs(median, IQR))

write.csv(test, "data/process/tables/inc_probs_crc_oral_residents_summary.csv", row.names = F)

pvalue_bacteria_groups <- tbl_df(c(
  wilcox.test(
  as.data.frame(filter(combined_residents, sampleType == "initial" & probs_increase == "No"))[, "oral_path"], 
  as.data.frame(filter(combined_residents, sampleType == "followup" & probs_increase == "No"))[, "oral_path"], 
  paired = TRUE)$p.value, 
  
  wilcox.test(
  as.data.frame(filter(combined_residents, sampleType == "initial" & probs_increase == "Yes"))[, "oral_path"], 
  as.data.frame(filter(combined_residents, sampleType == "followup" & probs_increase == "Yes"))[, "oral_path"], 
  paired = TRUE)$p.value, 
  
  wilcox.test(
  as.data.frame(filter(combined_residents, sampleType == "initial" & probs_increase == "No"))[, "crc_only"], 
  as.data.frame(filter(combined_residents, sampleType == "followup" & probs_increase == "No"))[, "crc_only"], 
  paired = TRUE)$p.value, 
  
  wilcox.test(
  as.data.frame(filter(combined_residents, sampleType == "initial" & probs_increase == "Yes"))[, "crc_only"], 
  as.data.frame(filter(combined_residents, sampleType == "followup" & probs_increase == "Yes"))[, "crc_only"], 
  paired = TRUE)$p.value, 
  
  wilcox.test(
    as.data.frame(filter(combined_residents, sampleType == "initial" & probs_increase == "No"))[, "normal_only"], 
    as.data.frame(filter(combined_residents, sampleType == "followup" & probs_increase == "No"))[, "normal_only"], 
    paired = TRUE)$p.value, 
  
  wilcox.test(
    as.data.frame(filter(combined_residents, sampleType == "initial" & probs_increase == "Yes"))[, "normal_only"], 
    as.data.frame(filter(combined_residents, sampleType == "followup" & probs_increase == "Yes"))[, "normal_only"], 
    paired = TRUE)$p.value)) %>% rename(pvalue = value) %>% 
  mutate(bh = p.adjust(pvalue, method = "BH"), 
         test = c("ivf_noinc_oralpath", "ivf_yesinc_oralpath", 
                  "ivf_noinc_crcassoc", "ivf_noinc_crcassoc", 
                  "ivf_noinc_normassoc", "ivf_noinc_normassoc"))

write.csv(pvalue_bacteria_groups, "data/process/tables/inc_probs_crc_oral_residents_pvalue_summary.csv", row.names = F)
