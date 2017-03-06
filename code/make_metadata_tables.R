### Make data analysis files 
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr"))


### Organize tables for train and test sets (first focus is to look at non-cancer versus cancer)
metaF <- read.delim('data/raw/metadata/followUps_metadata.txt', 
	header=T, sep='\t') %>% mutate(lesion = factor(NA, levels=c(0,1)))
metaF$cancer[metaF$dx =='normal' | metaF$dx =='adenoma'] <- 0
metaF$cancer[metaF$dx =='cancer'] <- 1
metaF$cancer <- factor(metaF$cancer)

metaI <- read.delim('data/raw/metadata/initials_metadata.tsv', 
	header=T, sep='\t') %>% mutate(lesion = factor(NA, levels=c(0,1)))
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


### Need to amend and separate Adenoma and CRC
good_metaf <- read.csv('data/raw/metadata/followUp_outcome_data.csv', 
	header = T, stringsAsFactors = F) %>% inner_join(metaF, by="EDRN")

metaFConly <- filter(good_metaf, Diagnosis == "adenocarcinoma" | 
	Diagnosis == "N/D")

metaFAonly <- filter(good_metaf, Diagnosis == "adenoma")


### Add finalized column of finalized lesion call the good_metaf
good_metaf$lesionf[good_metaf$Disease_Free == 'n'] <- 1
good_metaf$lesionf[good_metaf$Disease_Free == 'y' | good_metaf$Disease_Free == 'unknown'] <- 0


write.csv(metaI, "data/process/mod_metadata/metaI_final.csv", row.names = F)
write.csv(metaF, "data/process/mod_metadata/metaF_final.csv", row.names = F)
write.csv(good_metaf, "data/process/mod_metadata/good_metaf_final.csv", row.names = F)



