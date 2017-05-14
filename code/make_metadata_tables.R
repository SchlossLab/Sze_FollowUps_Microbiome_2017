### Make data analysis files 
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr"))

# Load needed data tables
dataI <- read.delim("data/process/metaI.txt", header = T, stringsAsFactors = F)
dataF <- read.delim("data/process/metaF.txt", header = T, stringsAsFactors = F)
extra_follow_samples <- read.csv("data/process/followup_samples.csv", header = T, stringsAsFactors = F)


# subset dataI and rename specific columns
metaI <- dataI %>% filter(grepl("mock", Sample_Name_s) == FALSE) %>% 
  mutate(dx = ifelse(diagnosis_s == "High Risk Normal" | diagnosis_s == "Normal", invisible("normal"), NA)) %>% 
  mutate(dx = ifelse(diagnosis_s == "Adenoma" | diagnosis_s == "adv Adenoma", invisible("adenoma"), dx)) %>% 
  mutate(dx = ifelse(is.na(dx), invisible("cancer"), dx)) %>% 
  select(Sample_Name_s, fit_result_s, diagnosis_s, dx, Hx_Prev_s, Hx_of_Polyps_s, Age_s, 
         Gender_s, Smoke_s, Diabetic_s, Hx_Fam_CRC_s, Height_s, Weight_s, BMI_s, White_s, 
         Native_s, Black_s, Pacific_s, Asian_s, Other_s, Ethnic_s, NSAID_s, Abx_s, 
         Diabetes_Med_s, cancer_stage_s) %>% 
  mutate(cancer_stage_s = ifelse(is.na(cancer_stage_s), 0, cancer_stage_s)) %>% 
  arrange(Sample_Name_s)

colnames(metaI) <- gsub("_s", "", colnames(metaI))

metaI <- metaI %>% 
  rename(sample = Sample_Name, Dx_Bin = diagnosis, stage = cancertage) %>% 
  mutate(sample = as.character(sample))

# create preliminary metaF data
metaF <- extra_follow_samples %>% 
  mutate(initial = as.character(initial), followUp = as.character(followUp)) %>% 
  inner_join(metaI, by = c("initial" = "sample")) %>% 
  select(-stage.y) %>% rename(stage = stage.x) %>% 
  distinct(initial, .keep_all = TRUE) %>% 
  mutate(Dx_Bin = ifelse(Dx_Bin == "Cancer", invisible("cancer"), Dx_Bin)) %>% 
  mutate(Dx_Bin = ifelse(Dx_Bin == "Adenoma", invisible("adenoma"), Dx_Bin)) %>% 
  mutate(Dx_Bin = ifelse(Dx_Bin == "adv Adenoma", invisible("adv_adenoma"), Dx_Bin))

#Select specific data from dataF file
mod_dataF <- dataF %>% 
  select(Sample_Name_s, fit_followUp_s, time_s, Disease_Free_s) %>% 
  distinct(Sample_Name_s, .keep_all = TRUE) %>% 
  mutate(Disease_Free_s = ifelse(Disease_Free_s == "y" | Disease_Free_s == "n", Disease_Free_s, invisible("unknown"))) %>% 
  mutate(Sample_Name_s = as.character(Sample_Name_s))
  
colnames(mod_dataF) <- gsub("_s", "", colnames(mod_dataF))

# Merge the two files together
metaF <- metaF %>% 
  inner_join(mod_dataF, by = c("followUp" = "Sample_Name"))

# Keep only distinct samples
metaI <- metaI %>% distinct(sample, .keep_all = TRUE)

### Organize tables for train and test sets (first focus is to look at non-cancer versus cancer)
metaF$cancer[metaF$dx =='normal' | metaF$dx =='adenoma'] <- 0
metaF$cancer[metaF$dx =='cancer'] <- 1
metaF$cancer <- factor(metaF$cancer)

metaI$cancer[metaI$dx == 'normal' | metaI$dx == 'adenoma'] <- 0
metaI$cancer[metaI$dx == 'cancer'] <- 1
metaI$cancer <- factor(metaI$cancer)

### Run Second organization to sort lesion vs. non-lesion
metaF$lesion[metaF$dx =='normal'] <- 0
metaF$lesion[metaF$dx =='adenoma' | metaF$dx == 'cancer'] <- 1
metaF$lesion <- factor(metaF$lesion)

metaI$lesion[metaI$dx == 'normal'] <- 0
metaI$lesion[metaI$dx == 'adenoma' | metaI$dx == 'cancer'] <- 1
metaI$lesion <- factor(metaI$lesion)

### Run third organization to sort SRN (Adv Adenoma) to cancer group and Adenoma to non-cancer group
metaF$SRNlesion[metaF$Dx_Bin =='normal' | metaF$Dx_Bin == 'Adenoma'] <- 0
metaF$SRNlesion[metaF$Dx_Bin =='adv Adenoma' | metaF$Dx_Bin == 'Cancer'] <- 1
metaF$SRNlesion <- factor(metaF$SRNlesion)

metaI$SRNlesion[metaI$Dx_Bin =='High Risk Normal' | metaI$Dx_Bin == 'Normal' | metaI$Dx_Bin == 'Adenoma'] <- 0
metaI$SRNlesion[metaI$Dx_Bin =='adv Adenoma' | metaI$Dx_Bin == 'Cancer'] <- 1
metaI$SRNlesion <- factor(metaI$SRNlesion)

### Run fourth organization to sort into three separate groups: normal, adenoma, cancer
metaF$threeway[metaF$dx == 'normal'] <- 0
metaF$threeway[metaF$dx == 'adenoma'] <- 1
metaF$threeway[metaF$dx == 'cancer'] <- 2
metaF$threeway <- factor(metaF$threeway)

metaI$threeway[metaI$dx == 'normal'] <- 0
metaI$threeway[metaI$dx == 'adenoma'] <- 1
metaI$threeway[metaI$dx == 'cancer'] <- 2
metaI$threeway <- factor(metaI$threeway)


### Create a fit positive or negative group (threshold set at 100)

metaI$fit_positive[metaI$fit_result > 100] <- 1
metaI$fit_positive[metaI$fit_result < 100] <- 0

metaF$fit_init_positive[metaF$fit_result > 100] <- 1
metaF$fit_init_positive[metaF$fit_result < 100] <- 0

metaF$fit_follow_positive[metaF$fit_followUp > 100] <- 1
metaF$fit_follow_positive[metaF$fit_followUp < 100] <- 0


### Add finalized column of finalized lesion call the good_metaf
metaF$lesionf[metaF$Disease_Free == 'n'] <- 1
metaF$lesionf[metaF$Disease_Free == 'y' | metaF$Disease_Free == 'unknown'] <- 0


write.csv(metaI, "data/process/mod_metadata/metaI_final.csv", row.names = F)
write.csv(metaF, "data/process/mod_metadata/metaF_final.csv", row.names = F)
write.csv(metaF, "data/process/mod_metadata/good_metaf_final.csv", row.names = F)



