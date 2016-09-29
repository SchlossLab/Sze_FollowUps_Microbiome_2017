## Run Previous Analysis align follow up data
## Code used and modified from Niel .Rmd file
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "foreach"))

### Read in necessary data 

tax <- read.delim(
	'data/process/old/followUps.final.an.unique_list.0.03.cons.taxonomy', 
                  sep='\t', header=T, row.names=1)
shared <- read.delim(
	'data/process/old/followUps.final.an.unique_list.0.03.subsample.0.03.filter.shared', 
                     header=T, sep='\t')

### Organize tables for train and test sets (first focus is to look at non-cancer versus cancer)
metaF <- read.delim('data/process/old/followUps_metadata.txt', header=T, sep='\t') %>% mutate(lesion = factor(NA, levels=c(0,1)))
metaF$cancer[metaF$dx =='normal' | metaF$dx =='adenoma'] <- 0
metaF$cancer[metaF$dx =='cancer'] <- 1
metaF$cancer <- factor(metaF$cancer)

metaI <- read.delim('data/process/old/initials_metadata.tsv', header=T, sep='\t') %>% mutate(lesion = factor(NA, levels=c(0,1)))
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

write.csv(metaI, "data/process/old/metaI_modified.csv")
write.csv(metaF, "data/process/old/metaF_modified.csv")

### Need to amend and separate Adenoma and CRC
good_metaf <- read.csv('data/process/old/followUp_outcome_data.csv', header = T, 
                       stringsAsFactors = F) %>% inner_join(metaF, by="EDRN")

metaFConly <- filter(good_metaf, Diagnosis == "adenocarcinoma" | Diagnosis == "N/D")
  
metaFAonly <- filter(good_metaf, Diagnosis == "adenoma")


###Pull out the thetayc distance between the initial and followup sample within the
###same person for all and split by adenoma or cancer

thetaCompTotal <- dissplit(
  'data/process/old/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist', metaF)

intra_ade <- as.vector(unlist(thetaCompTotal['intra_ade']))
intra_canc <- as.vector(unlist(thetaCompTotal['intra_canc']))
inter <- as.vector(unlist(thetaCompTotal['inter']))
rm(thetaCompTotal)

###Pull out the share OTUs between the initial and followup sample within the
###same person for all and split by adenoma or cancer

sobsCompTotal <- dissplit(
  'data/process/old/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist', metaF)

sobs_intra_ade <- as.vector(unlist(sobsCompTotal['intra_ade']))
sobs_intra_canc <- as.vector(unlist(sobsCompTotal['intra_canc']))
sobs_inter <- as.vector(unlist(sobsCompTotal['inter']))
rm(sobsCompTotal)

###Training and creating the models for prediction

# Create a list with all the data

testList <- list(
 cancer = train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
 filter(!sample %in% good_metaf$initial) %>% 
 select(cancer, fit_result, contains("Otu0")) %>% na.omit(), 
 lesion = train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
 filter(!sample %in% good_metaf$initial) %>% 
 select(lesion, fit_result, contains("Otu0")) %>% na.omit(),
 SRNlesion = train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
 filter(!sample %in% good_metaf$initial) %>% 
 select(SRNlesion, fit_result, contains("Otu0")) %>% na.omit(),
 threeway = train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
 filter(!sample %in% good_metaf$initial) %>% 
 select(threeway, fit_result, contains("Otu0")) %>% na.omit())


# Using the seperation as Cancer (ALL)
modelsToTest <- c("cancer", "lesion", "SRNlesion", "threeway")


# Set up parallelization
#doParallel::registerDoParallel(cores=5)

#finalDataResults <- foreach(i=1:length(modelsToTest)) %dopar% {

 # testdata(dataList, modelsToTest[i])
#}

### This doesnt work for the AUCRF function based on how it was written
### Can't seem to override the set functionality to look at global
### enivronment only.

# Sequentially 
cancer_AUC_data <- testdata(testList, modelsToTest[1])
lesion_AUC_data <- testdata(testList, modelsToTest[2])
SRNlesion_AUC_data <- testdata(testList, modelsToTest[3])
threeway_AUC_data <- testdata(testList, modelsToTest[4])

# Save the output for later
save.image("data/process/old/test2.RData")















