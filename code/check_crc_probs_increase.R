### Check up on What Six samples that increase are doing
### Create data files to graph these differences
## Marc Sze

## Load in necessary libraries
source('code/functions.R')
loadLibs(c("dplyr", "tidyr"))


## Load needed data files
good_metaf <- read.csv('data/process/mod_metadata/good_metaf_final.csv', header = T, stringsAsFactors = F)
crc_imp_vars <- read.csv('data/process/tables/crc_rf_otu_tax.csv', header = T, row.names = 1)
crc_probs <- read.csv('data/process/tables/crc_reduced_follow_up_probability_summary.csv', header = T)
shared <- read.delim('data/process/final.shared', header = T, stringsAsFactors = F) %>% select(-label, -numOtus)

## Get Rel. Abund
counts <- rowSums(select(shared, contains("Otu")))
otu.rel.abund <- shared
otu.rel.abund[, colnames(select(otu.rel.abund, contains("Otu")))] <- apply(
  select(otu.rel.abund, 
         contains("Otu")), 2, function(x) (x/counts)*100) 

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

# Test if any increase or decrease change is significant
change_pvalue_summary <- as.data.frame(
  matrix(nrow = length(imp_otus), 
         ncol = 7, 
         dimnames = list(rown = imp_otus, 
                         coln = c("otu", "pvalue", "bh", "median_Y", "median_N", "IQR_Y", "IQR_N")))) %>% 
  mutate(otu = imp_otus)

summary_test <- group_by(test, probs_increase, otu) %>% summarise_each(funs(median, IQR)) %>% slice(match(change_pvalue_summary$otu, otu))

for(i in 1:length(imp_otus)){
  
  change_pvalue_summary[i, "pvalue"] <- wilcox.test(filter(test_data, otu == imp_otus[i] & probs_increase == "Yes")[, "rel.abund"], 
                                                    filter(test_data, otu == imp_otus[i] & probs_increase == "No")[, "rel.abund"])$p.value
  
  change_pvalue_summary[i, "median_Y"] <- filter(summary_test, probs_increase == "Yes" & otu == imp_otus[i])[, "median"]
  change_pvalue_summary[i, "median_N"] <- filter(summary_test, probs_increase == "No" & otu == imp_otus[i])[, "median"]
  change_pvalue_summary[i, "IQR_Y"] <- filter(summary_test, probs_increase == "Yes" & otu == imp_otus[i])[, "IQR"]
  change_pvalue_summary[i, "IQR_N"] <- filter(summary_test, probs_increase == "No" & otu == imp_otus[i])[, "IQR"]
  
}

change_pvalue_summary <- mutate(change_pvalue_summary, bh = p.adjust(pvalue, method = "BH"))


## Check what the positivity is for the usual suspect OTUs for those with increased probs
usual_crc <- c("Otu000202", "Otu001273", "Otu000442")


## Add up Oral Pathogens and add up residents and compare
combined_residents <- tbl_df(combined_data$probs_increase) %>% rename(probs_increase = value) %>% 
  mutate(
    EDRN = combined_data$EDRN, 
    sampleType = combined_data$sampleType, 
    residents = (select(combined_data, one_of(imp_otus)) %>% select(-one_of(usual_crc)) %>% rowSums()), 
    oral_path = (select(combined_data, one_of(usual_crc)) %>% rowSums()))







