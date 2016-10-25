## Run Previous Analysis align follow up data
## Code used and modified from Niel .Rmd file
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))

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

# Using the seperation as Cancer (ALL)
train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  select(cancer, fit_result, contains("Otu0")) %>% na.omit()

set.seed(050416)
cancer_AUCRF_model<- AUCRF(cancer~., 
	data=train, pdel=0.05, ntree=500, ranking='MDA')
cancer_rf_opt <- cancer_AUCRF_model$RFopt
cancer_train_probs <- predict(cancer_rf_opt, type='prob')[,2]

cancer_train_roc <- roc(train$cancer ~ cancer_train_probs)

# Use cross validation to find the common OTUs
cancer_model1AUCRF_wCV <- AUCRFcv(cancer_AUCRF_model, 5, 5)
cancer_model1AUCRF_wCV_selected <- as.data.frame(
	cancer_model1AUCRF_wCV$Psel[cancer_model1AUCRF_wCV$Psel > 0.5])

write.csv(cancer_model1AUCRF_wCV_selected, 
	"data/process/old/cancer_AUCRF_wCV_selected.csv")

# Using the seperation as Lesion
train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  select(lesion, fit_result, contains("Otu0")) %>% na.omit()

set.seed(050416)
lesion_AUCRF_model <- AUCRF(lesion~., 
	data=train, pdel=0.05, ntree=500, ranking='MDA')
lesion_rf_opt <- lesion_AUCRF_model$RFopt
lesion_train_probs <- predict(lesion_rf_opt, type='prob')[,2]

lesion_train_roc <- roc(train$lesion ~ lesion_train_probs)


# Use cross validation to find the common OTUs
lesion_model1AUCRF_wCV <- AUCRFcv(lesion_AUCRF_model, 5, 5)
lesion_model1AUCRF_wCV_selected <- as.data.frame(
	lesion_model1AUCRF_wCV$Psel[lesion_model1AUCRF_wCV$Psel > 0.5])

write.csv(lesion_model1AUCRF_wCV_selected, 
	"data/process/old/lesion_AUCRF_wCV_selected.csv")


# Using the seperation as SRNlesion
train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  select(SRNlesion, fit_result, contains("Otu0")) %>% na.omit()

set.seed(050416)
SRNlesion_AUCRF_model <- AUCRF(SRNlesion~., 
	data=train, pdel=0.05, ntree=500, ranking='MDA')
SRNlesion_rf_opt <- SRNlesion_AUCRF_model$RFopt
SRNlesion_train_probs <- predict(SRNlesion_rf_opt, type='prob')[,2]

SRNlesion_train_roc <- roc(train$SRNlesion ~ SRNlesion_train_probs)


# Use cross validation to find the common OTUs


SRNlesion_model1AUCRF_wCV <- AUCRFcv(SRNlesion_AUCRF_model, 5, 5)
SRNlesion_model1AUCRF_wCV_selected <- as.data.frame(
	SRNlesion_model1AUCRF_wCV$Psel[SRNlesion_model1AUCRF_wCV$Psel > 0.5])

write.csv(SRNlesion_model1AUCRF_wCV_selected, 
	"data/process/old/SRNlesion_AUCRF_wCV_selected.csv")



# Using the seperation as threeway
train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  select(threeway, fit_result, contains("Otu0")) %>% na.omit()

set.seed(050416)
threeway_AUCRF_model <- AUCRF(threeway~., 
	data=train, pdel=0.05, ntree=500, ranking='MDA')
threeway_rf_opt <- threeway_AUCRF_model$RFopt

threeway_train_probs <- predict(threeway_rf_opt, type='prob')[,2]

threeway_train_roc <- roc(train$threeway ~ threeway_train_probs)


# Use cross validation to find the common OTUs


threeway_model1AUCRF_wCV <- AUCRFcv(threeway_AUCRF_model, 5, 5)
threeway_model1AUCRF_wCV_selected <- as.data.frame(
	threeway_model1AUCRF_wCV$Psel[threeway_model1AUCRF_wCV$Psel > 0.5])

write.csv(threeway_model1AUCRF_wCV_selected, 
	"data/process/old/threeway_AUCRF_wCV_selected.csv")

# Save the data image so I can upload it to R Studio later
save.image("data/process/old/test.RData")















