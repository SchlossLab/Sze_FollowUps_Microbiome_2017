## Run Previous Analysis align follow up data
## Code used and modified from Niel .Rmd file
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))

### Read in necessary data 

tax <- read.delim('data/process/followUps.final.an.unique_list.0.03.cons.taxonomy', 
                  sep='\t', header=T, row.names=1)
shared <- read.delim('data/process/followUps.final.an.unique_list.0.03.subsample.0.03.filter.shared', 
                     header=T, sep='\t')

### Organize tables for train and test sets (first focus is to look at non-cancer versus cancer)
metaF <- read.delim('data/process/followUps_metadata.txt', header=T, sep='\t') %>% mutate(lesion = factor(NA, levels=c(0,1)))
metaF$cancer[metaF$dx =='normal' | metaF$dx =='adenoma'] <- 0
metaF$cancer[metaF$dx =='cancer'] <- 1
metaF$cancer <- factor(metaF$cancer)

metaI <- read.delim('data/process/initials_metadata.tsv', header=T, sep='\t') %>% mutate(lesion = factor(NA, levels=c(0,1)))
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
good_metaf <- read.csv('data/process/followUp_outcome_data.csv', header = T, 
                       stringsAsFactors = F) %>% inner_join(metaF, by="EDRN")

metaFConly <- filter(good_metaf, Diagnosis == "adenocarcinoma" | Diagnosis == "N/D")
  
metaFAonly <- filter(good_metaf, Diagnosis == "adenoma")


write.csv(metaI, "results/tables/metaI_final.csv", row.names = F)
write.csv(metaF, "results/tables/metaF_final.csv", row.names = F)
write.csv(good_metaf, "results/tables/good_metaf_final.csv", row.names = F)

###Pull out the thetayc distance between the initial and followup sample within the
###same person for all and split by adenoma or cancer

thetaCompTotal <- dissplit(
  'data/process/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist', metaF)

intra_ade <- as.vector(unlist(thetaCompTotal['intra_ade']))
intra_canc <- as.vector(unlist(thetaCompTotal['intra_canc']))
inter <- as.vector(unlist(thetaCompTotal['inter']))
rm(thetaCompTotal)

###Pull out the share OTUs between the initial and followup sample within the
###same person for all and split by adenoma or cancer

sobsCompTotal <- dissplit(
  'data/process/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist', metaF)

sobs_intra_ade <- as.vector(unlist(sobsCompTotal['intra_ade']))
sobs_intra_canc <- as.vector(unlist(sobsCompTotal['intra_canc']))
sobs_inter <- as.vector(unlist(sobsCompTotal['inter']))
rm(sobsCompTotal)


###Training and creating the models for prediction

# createList with all data tables stored as a list

train <- list(cancer = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(cancer, fit_result, contains("Otu0")) %>% na.omit(), 
              lesion = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(lesion, fit_result, contains("Otu0")) %>% na.omit(), 
              SRNlesion = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(SRNlesion, fit_result, contains("Otu0")) %>% na.omit(), 
              threeway = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(threeway, fit_result, contains("Otu0")) %>% na.omit())


# Using the seperation as Cancer (ALL)
set.seed(050416)
cancer_AUCRF_model<- AUCRF(cancer~., data=train[["cancer"]], pdel=0.05, ntree=500, ranking='MDA')
cancer_rf_opt <- cancer_AUCRF_model$RFopt
cancer_train_probs <- predict(cancer_rf_opt, type='prob')[,2]

cancer_train_roc <- roc(train$cancer ~ cancer_train_probs)

# Use cross validation to find the common OTUs
cancer_model1AUCRF_wCV <- AUCRFcv(cancer_AUCRF_model, 5, 5)
cancer_model1AUCRF_wCV_selected <- as.data.frame(cancer_model1AUCRF_wCV$Psel[cancer_model1AUCRF_wCV$Psel > 0.5])


# ID important factors for Cancer using Boruta 

cancer_rf_model_used <- cancer_rf_opt$importance

selected_train <- select(train, cancer, one_of(rownames(cancer_rf_model_used)))

set.seed(050416)
cancer_impFactorData <- Boruta(cancer~., data=selected_train, mcAdj=TRUE, maxRuns=1000)
# Does not change after increasing runs to 2000

# Get the confirmed important variables
cancer_confirmed_vars <- as.data.frame(cancer_impFactorData['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

# write table to results to be used by other analysis components
write.csv(cancer_confirmed_vars, "results/tables/cancer_confirmed_vars.csv", row.names = F)

# Use the selected data set in AUCRF now
cancer_selected_train <- select(train, cancer, one_of(cancer_confirmed_vars[, 'otus']))
set.seed(050416)
cancer_selected_rf_opt <- AUCRF(cancer~., data=cancer_selected_train, pdel=0.99, ntree=500, ranking='MDA')$RFopt
cancer_selected_train_probs <- predict(cancer_selected_rf_opt, type='prob')[,2]

cancer_selected_train_roc <- roc(train$cancer ~ cancer_selected_train_probs)


# Using the seperation as Lesion
train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  select(lesion, fit_result, contains("Otu0")) %>% na.omit()

set.seed(050416)
lesion_rf_opt <- AUCRF(lesion~., data=train, pdel=0.05, ntree=500, ranking='MDA')$RFopt
lesion_train_probs <- predict(lesion_rf_opt, type='prob')[,2]

lesion_train_roc <- roc(train$lesion ~ lesion_train_probs)

# ID important factors for lesion using Boruta 
lesion_rf_model_used <- lesion_rf_opt$importance

selected_train <- select(train, lesion, one_of(rownames(lesion_rf_model_used)))

set.seed(050416)
lesion_impFactorData <- Boruta(lesion~., data=selected_train, mcAdj=TRUE, maxRuns=1000)
# Does not change after increasing runs to 2000

# Get the confirmed important variables
lesion_confirmed_vars <- as.data.frame(lesion_impFactorData['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

# write table to results to be used by other analysis components
write.csv(lesion_confirmed_vars, "results/tables/lesion_confirmed_vars.csv", row.names = F)

# Use the selected data set in AUCRF now
lesion_selected_train <- select(train, lesion, one_of(lesion_confirmed_vars[, 'otus']))
set.seed(050416)
lesion_selected_rf_opt <- AUCRF(lesion~., data=lesion_selected_train, pdel=0.99, ntree=500, ranking='MDA')$RFopt
lesion_selected_train_probs <- predict(lesion_selected_rf_opt, type='prob')[,2]

lesion_selected_train_roc <- roc(train$lesion ~ lesion_selected_train_probs)

# Using the seperation as SRNlesion
train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  select(SRNlesion, fit_result, contains("Otu0")) %>% na.omit()

set.seed(050416)
SRNlesion_rf_opt <- AUCRF(SRNlesion~., data=train, pdel=0.05, ntree=500, ranking='MDA')$RFopt
SRNlesion_train_probs <- predict(SRNlesion_rf_opt, type='prob')[,2]

SRNlesion_train_roc <- roc(train$SRNlesion ~ SRNlesion_train_probs)

# ID important factors for lesion using Boruta 
SRNlesion_rf_model_used <- SRNlesion_rf_opt$importance

selected_train <- select(train, SRNlesion, one_of(rownames(SRNlesion_rf_model_used)))

set.seed(050416)
SRNlesion_impFactorData <- Boruta(SRNlesion~., data=selected_train, mcAdj=TRUE, maxRuns=1000)
# Does not change after increasing runs to 2000

# Get the confirmed important variables
SRNlesion_confirmed_vars <- as.data.frame(SRNlesion_impFactorData['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

# write table to results to be used by other analysis components
write.csv(SRNlesion_confirmed_vars, "results/tables/SRNlesion_confirmed_vars.csv", row.names = F)

# Use the selected data set in AUCRF now
SRNlesion_selected_train <- select(train, SRNlesion, one_of(SRNlesion_confirmed_vars[, 'otus']))
set.seed(050416)
SRNlesion_selected_rf_opt <- AUCRF(SRNlesion~., data=SRNlesion_selected_train, pdel=0.99, ntree=500, ranking='MDA')$RFopt
SRNlesion_selected_train_probs <- predict(SRNlesion_selected_rf_opt, type='prob')[,2]

SRNlesion_selected_train_roc <- roc(train$SRNlesion ~ SRNlesion_selected_train_probs)


# Using the seperation as threeway
train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  select(threeway, fit_result, contains("Otu0")) %>% na.omit()

set.seed(050416)
threeway_rf_opt <- AUCRF(threeway~., data=train, pdel=0.05, ntree=500, ranking='MDA')$RFopt
threeway_train_probs <- predict(threeway_rf_opt, type='prob')[,2]

threeway_train_roc <- roc(train$threeway ~ threeway_train_probs)

# ID important factors for lesion using Boruta 
threeway_rf_model_used <- threeway_rf_opt$importance

selected_train <- select(train, threeway, one_of(rownames(threeway_rf_model_used)))

set.seed(050416)
threeway_impFactorData <- Boruta(threeway~., data=selected_train, mcAdj=TRUE, maxRuns=1000)
# Does not change after increasing runs to 2000

# Get the confirmed important variables
threeway_confirmed_vars <- as.data.frame(threeway_impFactorData['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

# write table to results to be used by other analysis components
write.csv(threeway_confirmed_vars, "results/tables/threeGroups_confirmed_vars.csv", row.names = F)

# Use the selected data set in AUCRF now
threeway_selected_train <- select(train, threeway, one_of(threeway_confirmed_vars[, 'otus']))
set.seed(050416)
threeway_selected_rf_opt <- AUCRF(threeway~., data=threeway_selected_train, pdel=0.99, ntree=500, ranking='MDA')$RFopt
threeway_selected_train_probs <- predict(threeway_selected_rf_opt, type='prob')[,2]

threeway_selected_train_roc <- roc(train$threeway ~ threeway_selected_train_probs)


### Graph the ROC curves for each of the different models and test for difference

# Created needed vectors and lists
rocNameList <- list(threeway_train_roc = threeway_train_roc, threeway_selected_train_roc = threeway_selected_train_roc, 
                    cancer_train_roc = cancer_train_roc, cancer_selected_train_roc = cancer_selected_train_roc, 
                    SRNlesion_train_roc = SRNlesion_train_roc, SRNlesion_selected_train_roc = SRNlesion_selected_train_roc, 
                    lesion_train_roc = lesion_train_roc, lesion_selected_train_roc = lesion_selected_train_roc)

variableList <- c("sensitivities", "specificities")
modelList <- c("ThreeGroupALL", "ThreeGroupSELECT", "cancerALL", "cancerSELECT", 
               "SRNlesionALL", "SRNlesionSELECT", "lesionALL", "lesionSELECT")

sens_specif_table <- makeSensSpecTable(rocNameList, variableList, modelList)
write.csv(sens_specif_table, "results/tables/ROCCurve_sens_spec.csv")

# Obtain the pvalue statistics as well as the bonferroni corrected values

corr_pvalue_ROC_table <- getROCPvalue(rocNameList, modelList, 8, multi = T)
write.csv(corr_pvalue_ROC_table, "results/tables/ROCCurve_corrected_pvalues.csv")

# Create the graph
ggplot(sens_specif_table, aes(sensitivities, specificities)) + 
  geom_line(aes(group = model, color = model), size = 1.5) + scale_x_continuous(trans = "reverse") + 
  theme_bw() + xlab("Sensitivity") + ylab("Specificity") + 
  theme(axis.title = element_text(face = "bold"))


### Investigate how these models do with initial and follow up samples

# Get threshold used by each of the different models to make the call

cutoffs <- as.data.frame(unlist(lapply(rocNameList, function(y) coords(y, x='best', ret='thr')))) %>% 
  cbind(rownames(.), .) %>%  mutate(model = c(rep("threeGroups", 2), rep("cancer", 2), rep("SRNlesion", 2), rep("lesion", 2))) %>% 
  mutate(dataset = rep(c("All", "Select"), length(rownames(.))/2))
            #take the threshold point used for the best sensitivity and specificty

colnames(cutoffs) <- c("L1", "cutoff", "model", "dataset")
cutoffs$model <- factor(cutoffs$model, levels = c("threeGroups", "cancer", "SRNlesion", "lesion"))
    
write.csv(cutoffs, "results/tables/withFIT.cutoffs.csv", row.names = F)

# Create data frames to be used for initial and follow up samples

initial <- inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
  select(threeway, cancer, lesion, SRNlesion, fit_result, contains("Otu0")) %>% rename(threeGroup = threeway)


good_metaf$lesionf[good_metaf$Disease_Free == 'n'] <- 1
good_metaf$lesionf[good_metaf$Disease_Free == 'y' | good_metaf$Disease_Free == 'unknown'] <- 0


followups <- inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
  select(lesionf, fit_result, contains("Otu0")) %>% rename(SRNlesion = lesionf)

followups <- cbind(good_metaf$lesionf, good_metaf$lesionf, good_metaf$lesionf, followups)
colnames(followups)[1:3] <- c("threeGroup", "cancer", "lesion")

# Get the prediction tables for each group using initial and follow up data

rf_opt_NameList <- list(threeway_rf_opt = threeway_rf_opt, threeway_selected_rf_opt = threeway_selected_rf_opt, 
                        cancer_rf_opt = cancer_rf_opt, cancer_selected_rf_opt = cancer_selected_rf_opt, 
                    SRNlesion_rf_opt = SRNlesion_rf_opt, SRNlesion_selected_rf_opt = SRNlesion_selected_rf_opt, 
                    lesion_rf_opt= lesion_rf_opt, lesion_selected_rf_opt = lesion_selected_rf_opt)

initial_predictions <- lapply(rf_opt_NameList, function(x) predict(x, initial, type='prob'))

followup_predictions <- lapply(rf_opt_NameList, function(x) predict(x, followups, type='prob'))

# Join good_metaf table with metaF table
combined_metaData <- inner_join(good_metaf, metaF, by = "EDRN")


######## can use full data set since predict function automatically selects variables from rf object to be used
######## Can use melt with defaults to create a ggplot useable table (0 is neg, 1 is pos)

df_initial_neg <- melt(initial_predictions) %>% filter(Var2 == 0)
df_initial_pos <- melt(initial_predictions) %>%  filter(Var2 == 1)

df_initial_preds <- cbind(rename(select(df_initial_neg, value), negative = value), 
                          rename(select(df_initial_pos, value, L1), positive = value))

rm(df_initial_neg, df_initial_pos)

df_initial_preds$model[grep("cancer_", df_initial_preds$L1)] <- "cancer"
df_initial_preds$model[grep("lesion", df_initial_preds$L1)] <- "lesion"
df_initial_preds$model[grep("SRN", df_initial_preds$L1)] <- "SRNlesion"
df_initial_preds$model[grep("three", df_initial_preds$L1)] <- "threeGroups"

df_initial_preds$dataset[grep("rf", df_initial_preds$L1)] <- "All"
df_initial_preds$dataset[grep("selected", df_initial_preds$L1)] <- "Select"
df_initial_preds <- mutate(df_initial_preds, 
                           diseaseFree = c(rep("n", 67*8))) %>% mutate(diagnosis = rep(good_metaf$Diagnosis, 8))

df_initial_preds$time_point <- "initial"

df_followups_neg <- melt(followup_predictions) %>% filter(Var2 == 0)
df_followups_pos <- melt(followup_predictions) %>% filter(Var2 == 1)

df_followups_preds <- cbind(rename(select(df_followups_neg, value), negative = value), 
                          rename(select(df_followups_pos, value, L1), positive = value))

rm(df_followups_neg, df_followups_pos)

df_followups_preds$model[grep("cancer_", df_followups_preds$L1)] <- "cancer"
df_followups_preds$model[grep("lesion", df_followups_preds$L1)] <- "lesion"
df_followups_preds$model[grep("SRN", df_followups_preds$L1)] <- "SRNlesion"
df_followups_preds$model[grep("three", df_followups_preds$L1)] <- "threeGroups"

df_followups_preds$dataset[grep("rf", df_followups_preds$L1)] <- "All"
df_followups_preds$dataset[grep("selected", df_followups_preds$L1)] <- "Select"
df_followups_preds <- mutate(df_followups_preds, diseaseFree = rep(good_metaf$Disease_Free, 8)) %>% 
  mutate(diagnosis = rep(good_metaf$Diagnosis, 8))

df_followups_preds$time_point <- "followup"

df_InitFollow_ALL <- rbind(df_initial_preds, df_followups_preds)
df_InitFollow_ALL$model <- factor(df_InitFollow_ALL$model, levels = c("threeGroups", "cancer", "SRNlesion", "lesion"))

df_InitFollow_ALL <- mutate(df_InitFollow_ALL, 
                            detailed_diagnosis = 
                              rep(combined_metaData$Dx_Bin.y, 
                                  length(rownames(df_InitFollow_ALL))/length(rownames(combined_metaData)))) %>% 
  mutate(EDRN = rep(good_metaf$EDRN, 16))

write.csv(df_InitFollow_ALL, "results/tables/withFIT.models.datatable.csv", row.names = F)


# Create labels for subset of data on graph

Names_facet <- c('SRNlesion' = "Classify SRN + Cancer", 'lesion' = "Classify Lesion")

# Graph the Adenoma data only

grid.arrange(
  # Graph the adenoma ALL data only
  filter(df_InitFollow_ALL, diagnosis == "adenoma" & dataset == "All" & model != "threeGroups" & model != "cancer") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive, group = factor(EDRN))) + 
    geom_point(aes(color=detailed_diagnosis)) + geom_line(aes(color = detailed_diagnosis)) + 
    scale_color_manual(name = "Polyp Type", values = c("cyan", "blue"), 
                       breaks = c("Adenoma", "adv Adenoma"), labels = c("Adenoma", "SRN")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(cutoffs, dataset == "All", model != "cancer", model != "threeGroups"), 
               aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Adenomas (Train on All Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")), 
  
  # Graph adenoma select data only
  filter(df_InitFollow_ALL, diagnosis == "adenoma" & dataset == "Select" & model != "threeGroups" & model != "cancer") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive, group = factor(EDRN))) + 
    geom_point(aes(color=detailed_diagnosis)) + geom_line(aes(color = detailed_diagnosis)) + 
    scale_color_manual(name = "Polyp Type", values = c("cyan", "blue"), 
                       breaks = c("Adenoma", "adv Adenoma"), labels = c("Adenoma", "SRN")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(cutoffs, dataset == "All", model != "cancer", model != "threeGroups"), 
               aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Adenomas (Train on Select Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")) 
  
)


# Graph the Cancer data only

grid.arrange(
  
  # Graph the cancer ALL data only
  filter(df_InitFollow_ALL, diagnosis != "adenoma" & dataset == "All" & model != "threeGroups" & model != "cancer") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive, group = factor(EDRN))) + 
    geom_point(aes(color=factor(diseaseFree, levels = c("n", "y", "unknown")))) + 
    geom_line(linetype = 2, alpha = 0.3) + 
    scale_color_manual(name = "Cancer\nFree", values = wes_palette("GrandBudapest")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(cutoffs, dataset == "All", model != "cancer", model != "threeGroups"), 
               aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Cancer (Train on All Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")), 
  
    # Graph cancer select data only
  filter(df_InitFollow_ALL, diagnosis != "adenoma" & dataset == "Select" & model != "threeGroups" & model != "cancer") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive, group = factor(EDRN))) + 
    geom_point(aes(color=factor(diseaseFree, levels = c("n", "y", "unknown")))) + 
    geom_line(linetype = 2, alpha = 0.3) + 
    scale_color_manual(name = "Cancer\nFree", values = wes_palette("GrandBudapest")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(cutoffs, dataset == "All", model != "cancer", model != "threeGroups"), 
               aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Cancer (Train on Select Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold"))
)





###Pull out the important OTUs from previous manuscript and make a ggplot format table 
###for those with follow up###

# Data tables that hold the information on the important features
# cancer_confirmed_vars, lesion_confirmed_vars, SRNlesion_confirmed_vars, threeway_confirmed_vars

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
    # This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components

# create needed labels for Boruta picked important variables for each model

cancer_selected_taxa <- tax_df[filter(cancer_confirmed_vars, otus != "fit_result")[, 'otus'], ]
cancer_selected_labs <- createTaxaLabeller(cancer_selected_taxa)

SRNlesion_selected_taxa <- tax_df[filter(SRNlesion_confirmed_vars, otus != "fit_result")[, 'otus'], ]
SRNlesion_selected_labs <- createTaxaLabeller(SRNlesion_selected_taxa)

lesion_selected_taxa <- tax_df[filter(lesion_confirmed_vars, otus != "fit_result")[, 'otus'], ]
lesion_selected_labs <- createTaxaLabeller(lesion_selected_taxa)

threeway_selected_taxa <- tax_df[filter(threeway_confirmed_vars, otus != "fit_result")[, 'otus'], ]
threeway_selected_labs <- createTaxaLabeller(threeway_selected_taxa)

cancer_positive <- as.data.frame(cbind(rep(as.character(good_metaf$Disease_Free), 36), rep(good_metaf$EDRN, 36))) %>% 
  rename(Disease_Free = V1, EDRN = V2)


# Graph the cancer model

cancer_model_impf_graphs <- follow_Abundance_Fit_Graph(cancer_selected_taxa, cancer_selected_labs, 
                                                       cancer_positive, initial, followups, good_metaf, 0.08, 67, "cancer", "Disease_Free", 
                                                       "Cancer")

grid.arrange(cancer_model_impf_graphs[['adenoma_OTUs']], cancer_model_impf_graphs[['cancer_OTUs']], 
             cancer_model_impf_graphs[['adenoma_fit']], cancer_model_impf_graphs[['cancer_fit']])

# Graph of the SRN lesion model

SRN_model_impf_graphs <- follow_Abundance_Fit_Graph(SRNlesion_selected_taxa, SRNlesion_selected_labs, 
                                                       cancer_positive, initial, followups, good_metaf, 0.08, 67, "SRNlesion", "Disease_Free", 
                                                    "Cancer")

grid.arrange(SRN_model_impf_graphs[['adenoma_OTUs']], SRN_model_impf_graphs[['cancer_OTUs']], 
             SRN_model_impf_graphs[['adenoma_fit']], SRN_model_impf_graphs[['cancer_fit']])



# Graph of the lesion model

lesion_model_impf_graphs <- follow_Abundance_Fit_Graph(lesion_selected_taxa, lesion_selected_labs, 
                                                    cancer_positive, initial, followups, good_metaf, 0.08, 67, "lesion")

grid.arrange(lesion_model_impf_graphs[['adenoma_OTUs']], lesion_model_impf_graphs[['cancer_OTUs']], 
             lesion_model_impf_graphs[['adenoma_fit']], lesion_model_impf_graphs[['cancer_fit']])


# Graph of the three group model


three_model_impf_graphs <- follow_Abundance_Fit_Graph(threeway_selected_taxa, threeway_selected_labs, 
                                                      cancer_positive, initial, followups, good_metaf %>% rename(threeGroup = threeway), 
                                                      0.08, 67, "threeGroup")

grid.arrange(three_model_impf_graphs[['adenoma_OTUs']], three_model_impf_graphs[['cancer_OTUs']], 
             three_model_impf_graphs[['adenoma_fit']], three_model_impf_graphs[['cancer_fit']])


# createList with all data tables stored as a list 
# use fit as a binary group rather than a continuous variable

fit_group_train <- list(cancer = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(cancer, fit_positive, contains("Otu0")) %>% na.omit(), 
              lesion = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(lesion, fit_positive, contains("Otu0")) %>% na.omit(), 
              SRNlesion = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(SRNlesion, fit_positive, contains("Otu0")) %>% na.omit(), 
              threeway = inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
                filter(!sample %in% good_metaf$initial) %>% 
                select(threeway, fit_positive, contains("Otu0")) %>% na.omit())



# Using the seperation as Lesion 
set.seed(050416)
fitGroup_lesion_AUCRF_model<- AUCRF(lesion~., data=fit_group_train[["lesion"]], pdel=0.05, ntree=500, ranking='MDA')
fitGroup_lesion_rf_opt <- fitGroup_lesion_AUCRF_model$RFopt
fitGroup_lesion_train_probs <- predict(fitGroup_lesion_rf_opt, type='prob')[,2]

fitGroup_lesion_train_roc <- roc(fit_group_train[["lesion"]]$lesion ~ fitGroup_lesion_train_probs)


# ID important factors for Cancer using Boruta 

fitGroup_lesion_rf_model_used <- fitGroup_lesion_rf_opt$importance

selected_train <- select(fit_group_train[["lesion"]], lesion, one_of(rownames(fitGroup_lesion_rf_model_used)))

set.seed(050416)
fitGroup_lesion_impFactorData <- Boruta(lesion~., data=selected_train, mcAdj=TRUE, maxRuns=1000)
# Does not change after increasing runs to 2000

# Get the confirmed important variables
fitGroup_lesion_confirmed_vars <- as.data.frame(fitGroup_lesion_impFactorData['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

# write table to results to be used by other analysis components
write.csv(fitGroup_lesion_confirmed_vars, "results/tables/fitGroup_lesion_confirmed_vars.csv", row.names = F)

# Use the selected data set in AUCRF now
fitGroup_lesion_selected_train <- select(fit_group_train[["lesion"]], 
                                         lesion, one_of(fitGroup_lesion_confirmed_vars[, 'otus']))
set.seed(050416)
fitGroup_lesion_selected_rf_opt <- AUCRF(lesion~., data=fitGroup_lesion_selected_train, 
                                         pdel=0.99, ntree=500, ranking='MDA')$RFopt

fitGroup_lesion_selected_train_probs <- predict(fitGroup_lesion_selected_rf_opt, type='prob')[,2]

fitGroup_lesion_selected_train_roc <- roc(fitGroup_lesion_selected_train$lesion ~ fitGroup_lesion_selected_train_probs)



# Using the seperation as SRNLesion 
set.seed(050416)
fitGroup_SRNlesion_AUCRF_model<- AUCRF(SRNlesion~., data=fit_group_train[["SRNlesion"]], pdel=0.05, ntree=500, ranking='MDA')
fitGroup_SRNlesion_rf_opt <- fitGroup_SRNlesion_AUCRF_model$RFopt
fitGroup_SRNlesion_train_probs <- predict(fitGroup_SRNlesion_rf_opt, type='prob')[,2]

fitGroup_SRNlesion_train_roc <- roc(fit_group_train[["SRNlesion"]]$SRNlesion ~ fitGroup_SRNlesion_train_probs)


# ID important factors for Cancer using Boruta 

fitGroup_SRNlesion_rf_model_used <- fitGroup_SRNlesion_rf_opt$importance

selected_train <- select(fit_group_train[["SRNlesion"]], SRNlesion, one_of(rownames(fitGroup_SRNlesion_rf_model_used)))

set.seed(050416)
fitGroup_SRNlesion_impFactorData <- Boruta(SRNlesion~., data=selected_train, mcAdj=TRUE, maxRuns=1000)
# Does not change after increasing runs to 2000

# Get the confirmed important variables
fitGroup_SRNlesion_confirmed_vars <- as.data.frame(fitGroup_SRNlesion_impFactorData['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

# write table to results to be used by other analysis components
write.csv(fitGroup_SRNlesion_confirmed_vars, "results/tables/fitGroup_SRNlesion_confirmed_vars.csv", row.names = F)

# Use the selected data set in AUCRF now
fitGroup_SRNlesion_selected_train <- select(fit_group_train[["SRNlesion"]], 
                                         SRNlesion, one_of(fitGroup_SRNlesion_confirmed_vars[, 'otus']))
set.seed(050416)
fitGroup_SRNlesion_selected_rf_opt <- AUCRF(SRNlesion~., data=fitGroup_SRNlesion_selected_train, 
                                         pdel=0.99, ntree=500, ranking='MDA')$RFopt

fitGroup_SRNlesion_selected_train_probs <- predict(fitGroup_SRNlesion_selected_rf_opt, type='prob')[,2]

fitGroup_SRNlesion_selected_train_roc <- roc(
  fitGroup_SRNlesion_selected_train$SRNlesion ~ fitGroup_SRNlesion_selected_train_probs)


# Using the seperation as Three Way 
set.seed(050416)
fitGroup_threeway_AUCRF_model<- AUCRF(threeway~., data=fit_group_train[["threeway"]], pdel=0.05, ntree=500, ranking='MDA')
fitGroup_threeway_rf_opt <- fitGroup_threeway_AUCRF_model$RFopt
fitGroup_threeway_train_probs <- predict(fitGroup_threeway_rf_opt, type='prob')[,2]

fitGroup_threeway_train_roc <- roc(fit_group_train[["threeway"]]$threeway ~ fitGroup_threeway_train_probs)


# ID important factors for three way using Boruta 

fitGroup_threeway_rf_model_used <- fitGroup_threeway_rf_opt$importance

selected_train <- select(fit_group_train[["threeway"]], threeway, one_of(rownames(fitGroup_threeway_rf_model_used)))

set.seed(050416)
fitGroup_threeway_impFactorData <- Boruta(threeway~., data=selected_train, mcAdj=TRUE, maxRuns=1000)
# Does not change after increasing runs to 2000

# Get the confirmed important variables
fitGroup_threeway_confirmed_vars <- as.data.frame(fitGroup_threeway_impFactorData['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

# write table to results to be used by other analysis components
write.csv(fitGroup_threeway_confirmed_vars, "results/tables/fitGroup_threeway_confirmed_vars.csv", row.names = F)

# Use the selected data set in AUCRF now
fitGroup_threeway_selected_train <- select(fit_group_train[["threeway"]], 
                                            threeway, one_of(fitGroup_threeway_confirmed_vars[, 'otus']))
set.seed(050416)
fitGroup_threeway_selected_rf_opt <- AUCRF(threeway~., data=fitGroup_threeway_selected_train, 
                                            pdel=0.99, ntree=500, ranking='MDA')$RFopt

fitGroup_threeway_selected_train_probs <- predict(fitGroup_threeway_selected_rf_opt, type='prob')[,2]

fitGroup_threeway_selected_train_roc <- roc(
  fitGroup_threeway_selected_train$threeway ~ fitGroup_threeway_selected_train_probs)


### Graph the ROC curves for each of the different models and test for difference

# Created needed vectors and lists
fitGroup_rocNameList <- list(lesion = fitGroup_lesion_train_roc, lesion_selected = fitGroup_lesion_selected_train_roc, 
                    SRNlesion = fitGroup_SRNlesion_train_roc, SRNlesion_selected = fitGroup_SRNlesion_selected_train_roc, 
                    threeway = fitGroup_threeway_train_roc, threeway_selected = fitGroup_threeway_selected_train_roc)

variableList <- c("sensitivities", "specificities")
modelList <- c("lesionALL", "lesionSELECT", "SRNlesionALL", "SRNlesionSELECT", "ThreeGroupALL", "ThreeGroupSELECT")

sens_specif_table <- makeSensSpecTable(fitGroup_rocNameList, variableList, modelList)
write.csv(sens_specif_table, "results/tables/fitGroup_ROCCurve_sens_spec.csv")

# Obtain the pvalue statistics as well as the bonferroni corrected values

corr_pvalue_ROC_table <- getROCPvalue(fitGroup_rocNameList, modelList, 6, multi = T)
write.csv(corr_pvalue_ROC_table, "results/tables/fitGroup_ROCCurve_corrected_pvalues.csv")

# Create the graph
ggplot(sens_specif_table, aes(sensitivities, specificities)) + 
  geom_line(aes(group = model, color = model), size = 1.5) + scale_x_continuous(trans = "reverse") + 
  theme_bw() + xlab("Sensitivity") + ylab("Specificity") + 
  theme(axis.title = element_text(face = "bold"))




### Investigate how these models do with initial and follow up samples with FIT as a group variable

# Get threshold used by each of the different models to make the call

fitgroup_cutoffs <- as.data.frame(unlist(lapply(fitGroup_rocNameList, function(y) coords(y, x='best', ret='thr')))) %>% 
  cbind(rownames(.), .) %>%  mutate(model = c(rep("threeGroups", 2), rep("SRNlesion", 2), rep("lesion", 2))) %>% 
  mutate(dataset = rep(c("All", "Select"), length(rownames(.))/2)) %>% filter(model != "threeGroups")
#take the threshold point used for the best sensitivity and specificty

colnames(fitgroup_cutoffs) <- c("L1", "cutoff", "model", "dataset")
fitgroup_cutoffs$model <- factor(fitgroup_cutoffs$model, levels = c("SRNlesion", "lesion"))

write.csv(fitgroup_cutoffs, "results/tables/withFITGroups.cutoffs.csv", row.names = F)

# Create data frames to be used for initial and follow up samples

initial <- inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
  select(threeway, lesion, SRNlesion, fit_init_positive, contains("Otu0")) %>% 
  rename(threeGroup = threeway, fit_positive = fit_init_positive)


good_metaf$lesionf[good_metaf$Disease_Free == 'n'] <- 1
good_metaf$lesionf[good_metaf$Disease_Free == 'y' | good_metaf$Disease_Free == 'unknown'] <- 0


followups <- inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
  select(lesionf, fit_follow_positive, contains("Otu0")) %>% 
  rename(SRNlesion = lesionf, fit_positive = fit_follow_positive)

followups <- cbind(good_metaf$lesionf, good_metaf$lesionf, good_metaf$lesionf, followups)
colnames(followups)[1:3] <- c("threeGroup", "cancer", "lesion")

# Get the prediction tables for each group using initial and follow up data

fitGroup_threeway_selected_rf_opt


fitGroup_rf_opt_NameList <- list(threeway = fitGroup_threeway_rf_opt, threeway_selected = fitGroup_threeway_selected_rf_opt, 
                        SRNlesion = fitGroup_SRNlesion_rf_opt, SRNlesion_selected = fitGroup_SRNlesion_selected_rf_opt, 
                        lesion = fitGroup_lesion_rf_opt, lesion_selected = fitGroup_lesion_selected_rf_opt)

initial_predictions <- lapply(fitGroup_rf_opt_NameList, function(x) predict(x, initial, type='prob'))

followup_predictions <- lapply(fitGroup_rf_opt_NameList, function(x) predict(x, followups, type='prob'))

# Join good_metaf table with metaF table
combined_metaData <- inner_join(good_metaf, metaF, by = "EDRN")


######## can use full data set since predict function automatically selects variables from rf object to be used
######## Can use melt with defaults to create a ggplot useable table (0 is neg, 1 is pos)

df_initial_neg <- melt(initial_predictions) %>% filter(Var2 == 0)
df_initial_pos <- melt(initial_predictions) %>%  filter(Var2 == 1)

df_initial_preds <- cbind(rename(select(df_initial_neg, value), negative = value), 
                          rename(select(df_initial_pos, value, L1), positive = value))

rm(df_initial_neg, df_initial_pos)

df_initial_preds$model[grepl("lesion", df_initial_preds$L1)] <- "lesion"
df_initial_preds$model[grepl("SRN", df_initial_preds$L1)] <- "SRNlesion"
df_initial_preds$model[grepl("three", df_initial_preds$L1)] <- "threeGroups"

df_initial_preds$dataset[!grepl("selected", df_initial_preds$L1)] <- "All"
df_initial_preds$dataset[grepl("selected", df_initial_preds$L1)] <- "Select"
df_initial_preds <- mutate(df_initial_preds, 
                           diseaseFree = c(rep("n", 67*6))) %>% mutate(diagnosis = rep(good_metaf$Diagnosis, 6))

    # Using grepl returns a logical (which I think should be preferred)

df_initial_preds$time_point <- "initial"

df_followups_neg <- melt(followup_predictions) %>% filter(Var2 == 0)
df_followups_pos <- melt(followup_predictions) %>% filter(Var2 == 1)

df_followups_preds <- cbind(rename(select(df_followups_neg, value), negative = value), 
                            rename(select(df_followups_pos, value, L1), positive = value))

rm(df_followups_neg, df_followups_pos)

df_followups_preds$model[grep("lesion", df_followups_preds$L1)] <- "lesion"
df_followups_preds$model[grep("SRN", df_followups_preds$L1)] <- "SRNlesion"
df_followups_preds$model[grep("three", df_followups_preds$L1)] <- "threeGroups"

df_followups_preds$dataset[!grepl("selected", df_followups_preds$L1)] <- "All"
df_followups_preds$dataset[grepl("selected", df_followups_preds$L1)] <- "Select"
df_followups_preds <- mutate(df_followups_preds, diseaseFree = rep(good_metaf$Disease_Free, 6)) %>% 
  mutate(diagnosis = rep(good_metaf$Diagnosis, 6))

df_followups_preds$time_point <- "followup"

df_InitFollow_ALL <- rbind(df_initial_preds, df_followups_preds)
df_InitFollow_ALL$model <- factor(df_InitFollow_ALL$model, levels = c("threeGroups", "SRNlesion", "lesion"))

df_InitFollow_ALL <- mutate(df_InitFollow_ALL, 
                            detailed_diagnosis = 
                              rep(combined_metaData$Dx_Bin.y, 
                                  length(rownames(df_InitFollow_ALL))/length(rownames(combined_metaData)))) %>% 
  mutate(EDRN = rep(good_metaf$EDRN, 12))

write.csv(df_InitFollow_ALL, "results/tables/FitGroups.models.datatable.csv", row.names = F)


# Create labels for subset of data on graph

Names_facet <- c('SRNlesion' = "Classify SRN + Cancer", 'lesion' = "Classify Lesion")

# Graph the Adenoma data only

grid.arrange(
  # Graph the adenoma ALL data only
  filter(df_InitFollow_ALL, diagnosis == "adenoma" & dataset == "All" & model != "threeGroups") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive, group = factor(EDRN))) + 
    geom_point(aes(color=detailed_diagnosis)) + geom_line(aes(color = detailed_diagnosis)) + 
    scale_color_manual(name = "Polyp Type", values = c("cyan", "blue"), 
                       breaks = c("Adenoma", "adv Adenoma"), labels = c("Adenoma", "SRN")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(fitgroup_cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Adenomas (Train on All Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")), 
  
  # Graph adenoma select data only
  filter(df_InitFollow_ALL, diagnosis == "adenoma" & dataset == "Select" & model != "threeGroups") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive, group = factor(EDRN))) + 
    geom_point(aes(color=detailed_diagnosis)) + geom_line(aes(color = detailed_diagnosis)) + 
    scale_color_manual(name = "Polyp Type", values = c("cyan", "blue"), 
                       breaks = c("Adenoma", "adv Adenoma"), labels = c("Adenoma", "SRN")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(fitgroup_cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Adenomas (Train on Select Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold"))
)


# Graph the Cancer data only

grid.arrange(
  
  # Graph the cancer ALL data only
  filter(df_InitFollow_ALL, diagnosis != "adenoma" & dataset == "All" & model != "threeGroups") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive, group= factor(EDRN))) + 
    geom_point(aes(color=factor(diseaseFree, levels = c("n", "y", "unknown")))) + 
    geom_line(linetype = 2, alpha = 0.3) + 
    scale_color_manual(name = "Cancer\nFree", values = wes_palette("GrandBudapest")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(fitgroup_cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Cancer (Train on All Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")),
  
  # Graph cancer select data only
  filter(df_InitFollow_ALL, diagnosis != "adenoma" & dataset == "Select" & model != "threeGroups") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive, group = factor(EDRN))) + 
    geom_point(aes(color=factor(diseaseFree, levels = c("n", "y", "unknown")))) + 
    geom_line(linetype = 2, alpha = 0.3) + 
    scale_color_manual(name = "Cancer\nFree", values = wes_palette("GrandBudapest")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(fitgroup_cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Cancer (Train on Select Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold"))
)





###Pull out theimportant OTUs from previous manuscript and make a ggplot format table 
###for those with follow up###

# Data tables that hold the information on the important features
# cancer_confirmed_vars, lesion_confirmed_vars, SRNlesion_confirmed_vars, threeway_confirmed_vars

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
# This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components

# create needed labels for Boruta picked important variables for each model

cancer_selected_taxa <- tax_df[filter(cancer_confirmed_vars, otus != "fit_result")[, 'otus'], ]
cancer_selected_labs <- createTaxaLabeller(cancer_selected_taxa)

SRNlesion_selected_taxa <- tax_df[filter(SRNlesion_confirmed_vars, otus != "fit_result")[, 'otus'], ]
SRNlesion_selected_labs <- createTaxaLabeller(SRNlesion_selected_taxa)

lesion_selected_taxa <- tax_df[filter(lesion_confirmed_vars, otus != "fit_result")[, 'otus'], ]
lesion_selected_labs <- createTaxaLabeller(lesion_selected_taxa)

threeway_selected_taxa <- tax_df[filter(threeway_confirmed_vars, otus != "fit_result")[, 'otus'], ]
threeway_selected_labs <- createTaxaLabeller(threeway_selected_taxa)

cancer_positive <- as.data.frame(cbind(rep(as.character(good_metaf$Disease_Free), 36), rep(good_metaf$EDRN, 36))) %>% 
  rename(Disease_Free = V1, EDRN = V2)


# Graph the cancer model

cancer_model_impf_graphs <- follow_Abundance_Fit_Graph(cancer_selected_taxa, cancer_selected_labs, 
                                                       cancer_positive, initial, followups, good_metaf, 0.08, 67, "cancer", "Disease_Free", 
                                                       "Cancer")

grid.arrange(cancer_model_impf_graphs[['adenoma_OTUs']], cancer_model_impf_graphs[['cancer_OTUs']], 
             cancer_model_impf_graphs[['adenoma_fit']], cancer_model_impf_graphs[['cancer_fit']])

# Graph of the SRN lesion model

SRN_model_impf_graphs <- follow_Abundance_Fit_Graph(SRNlesion_selected_taxa, SRNlesion_selected_labs, 
                                                    cancer_positive, initial, followups, good_metaf, 0.08, 67, "SRNlesion", "Disease_Free", 
                                                    "Cancer")

grid.arrange(SRN_model_impf_graphs[['adenoma_OTUs']], SRN_model_impf_graphs[['cancer_OTUs']], 
             SRN_model_impf_graphs[['adenoma_fit']], SRN_model_impf_graphs[['cancer_fit']])



# Graph of the lesion model

lesion_model_impf_graphs <- follow_Abundance_Fit_Graph(lesion_selected_taxa, lesion_selected_labs, 
                                                       cancer_positive, initial, followups, good_metaf, 0.08, 67, "lesion")

grid.arrange(lesion_model_impf_graphs[['adenoma_OTUs']], lesion_model_impf_graphs[['cancer_OTUs']], 
             lesion_model_impf_graphs[['adenoma_fit']], lesion_model_impf_graphs[['cancer_fit']])


# Graph of the three group model


three_model_impf_graphs <- follow_Abundance_Fit_Graph(threeway_selected_taxa, threeway_selected_labs, 
                                                      cancer_positive, initial, followups, good_metaf %>% rename(threeGroup = threeway), 
                                                      0.08, 67, "threeGroup")

grid.arrange(three_model_impf_graphs[['adenoma_OTUs']], three_model_impf_graphs[['cancer_OTUs']], 
             three_model_impf_graphs[['adenoma_fit']], three_model_impf_graphs[['cancer_fit']])



### Explore possibility of creating group variables for specific OTUs and whether that will improve classifications

# Create a list of the stored selected OTUs tables

boruta_selected_list_data <- list(lesion = fitGroup_lesion_selected_train, SRNlesion = fitGroup_SRNlesion_selected_train, 
                                  threeway = fitGroup_threeway_selected_train)

# Get the average value of 1% abundance (might be to strict)

abund_cutoff <- sum(
  select(shared, contains("Otu0")) %>% 
    rowSums())/length(select(shared, contains("Otu0")) %>% 
                        rowSums()) * 0.01

# First, look at simply using positive/negative as a cutoff
catLesion_table <- makeCat(boruta_selected_list_data[["lesion"]], cutoff = 0)

set.seed(050416)
catLesion_table_AUCRF_model<- AUCRF(lesion~., data=catLesion_table, pdel=0.99, ntree=500, ranking='MDA')
  # Throws out 14, 147, and 166 from equation
  # Lets convert these back to numeric

catLesion_table$Otu000014 <- boruta_selected_list_data[["lesion"]]$Otu000014
catLesion_table$Otu000147 <- boruta_selected_list_data[["lesion"]]$Otu000147
catLesion_table$Otu000166 <- boruta_selected_list_data[["lesion"]]$Otu000166

#try again to see if that improved the classification
set.seed(050416)
part_catLesion_table_AUCRF_model<- AUCRF(lesion~., data=catLesion_table, pdel=0.99, ntree=500, ranking='MDA')

#try a 1% cutoff
onePercent_catLesion_table <- makeCat(boruta_selected_list_data[["lesion"]], cutoff = abund_cutoff)
set.seed(050416)
onePercent_catLesion_table_AUCRF_model<- AUCRF(lesion~., data=onePercent_catLesion_table, pdel=0.99, ntree=500, ranking='MDA')
    # This definitely does not work


#Look at SRN Lesion model
SRNcatLesion_table <- makeCat(boruta_selected_list_data[["SRNlesion"]], cutoff = 0)
set.seed(050416)
catSRNLesion_table_AUCRF_model<- AUCRF(SRNlesion~., data=SRNcatLesion_table, pdel=0.99, ntree=500, ranking='MDA')
    # 166 and 162 are not useful as categories    
    # Lets convert these back to numeric

SRNcatLesion_table$Otu000162 <- boruta_selected_list_data[["SRNlesion"]]$Otu000162
SRNcatLesion_table$Otu000166 <- boruta_selected_list_data[["SRNlesion"]]$Otu000166
SRNcatLesion_table$Otu000003 <- boruta_selected_list_data[["SRNlesion"]]$Otu000003
SRNcatLesion_table$Otu000005 <- boruta_selected_list_data[["SRNlesion"]]$Otu000005
SRNcatLesion_table$Otu000008 <- boruta_selected_list_data[["SRNlesion"]]$Otu000008

#try again to see if that improved the classification
set.seed(050416)
part_catSRNLesion_table_AUCRF_model<- AUCRF(SRNlesion~., data=SRNcatLesion_table, pdel=0.99, ntree=500, ranking='MDA')

# convert the higher abundance OTUs back to continuous and remove 162
SRNcatLesion_table <- subset(SRNcatLesion_table, select = -c(Otu000162, Otu000003))

set.seed(050416)
part_catSRNLesion_table_AUCRF_model<- AUCRF(SRNlesion~., data=SRNcatLesion_table, pdel=0.99, ntree=500, ranking='MDA')

#try a 1% cutoff
onePercent_SRNcatLesion_table <- makeCat(boruta_selected_list_data[["SRNlesion"]], cutoff = abund_cutoff)
set.seed(050416)
onePercent_SRNcatLesion_table_AUCRF_model<- AUCRF(SRNlesion~., 
                                               data=onePercent_SRNcatLesion_table, pdel=0.99, ntree=500, ranking='MDA')
# This makes it worse but perhaps adding some as 1% cutoff will improve
    # Replace 3, 5, and 8 with 1% cutoffs  since in previous tests it results in worse predictions

SRNcatLesion_table <- makeCat(boruta_selected_list_data[["SRNlesion"]], cutoff = 0)
SRNcatLesion_table$Otu000003 <- onePercent_SRNcatLesion_table$Otu000003
SRNcatLesion_table$Otu000005 <- onePercent_SRNcatLesion_table$Otu000005
SRNcatLesion_table$Otu000008 <- onePercent_SRNcatLesion_table$Otu000008

onePercent_Part_SRNcatLesion_table_AUCRF_model<- AUCRF(SRNlesion~., 
                                                       data=SRNcatLesion_table, pdel=0.99, ntree=500, ranking='MDA')

# remove 162 and 147 due to worse predictions
SRNcatLesion_table <- subset(SRNcatLesion_table, select = -c(Otu000162, Otu000147))
set.seed(050416)
onePercent_Part_SRNcatLesion_table_AUCRF_model<- AUCRF(SRNlesion~., 
                                                       data=SRNcatLesion_table, pdel=0.99, ntree=500, ranking='MDA')

# use fit as continuous not categorical
SRNcatLesion_table$fit_positive <- metaI$fit_result
set.seed(050416)
onePercent_Part_SRNcatLesion_table_AUCRF_model<- AUCRF(SRNlesion~., 
                                                       data=SRNcatLesion_table, pdel=0.99, ntree=500, ranking='MDA')

# Perhaps adding number of known associated OTUs with CRC will help prediction
SRNlesion_selected_taxa <- tax_df[filter(SRNlesion_confirmed_vars, otus != "fit_positive")[, 'otus'], ]
SRNlesion_selected_labs <- createTaxaLabeller(SRNlesion_selected_taxa)


knownCRC_vars <- c("Otu000566", "Otu000126", "Otu000205", "Otu000397", "Otu000563")

SRNcatLesion_table$crc_bugs <- apply(SRNcatLesion_table[, knownCRC_vars], 1, function(x) getCount(x))
set.seed(050416)
onePercent_Part_SRNcatLesion_table_AUCRF_model<- AUCRF(SRNlesion~., 
                                                       data=SRNcatLesion_table, pdel=0.99, ntree=500, ranking='MDA')

  ## Still have 70 that are misclassified in the SRN + Cancer lesion group can this be improved with meta data on
  ## other risk factors?

# First convert crc_bugs back to sequence counts

SRNcatLesion_table$Otu000566 <- boruta_selected_list_data[["SRNlesion"]]$Otu000566
SRNcatLesion_table$Otu000126 <- boruta_selected_list_data[["SRNlesion"]]$Otu000126
SRNcatLesion_table$Otu000205 <- boruta_selected_list_data[["SRNlesion"]]$Otu000205
SRNcatLesion_table$Otu000397 <- boruta_selected_list_data[["SRNlesion"]]$Otu000397
SRNcatLesion_table$Otu000563 <- boruta_selected_list_data[["SRNlesion"]]$Otu000563


set.seed(050416)
onePercent_Part_SRNcatLesion_table_AUCRF_model<- AUCRF(SRNlesion~., 
                                                       data=SRNcatLesion_table, pdel=0.99, ntree=500, ranking='MDA')

# add metadata risk factors for CRC into model
# additions are 
  # overweight versus obese (normal (18.5-24.9) = 0, overweight (25.0 - 29.9) = 1, obese (> 30) = 2)
  metaI$BMIcat[metaI$BMI < 25] <- 0  
  metaI$BMIcat[metaI$BMI > 25 & metaI$BMI < 30] <- 1
  metaI$BMIcat[metaI$BMI > 30] <- 2

  # smoking (smoke variable), missing a variable so need to use na.action (already exists)
  
  # Age greater than 50 variable
  metaI$great50[metaI$Age <= 50] <- 0
  metaI$great50[metaI$Age > 50] <- 1
  
  # history of polyps (Hx_of_Polyps = 0, 1)  (already exists)
  # history of cancer (Hx_Prev) (already exists)
  # family history of crc (Hx_Fam_CRC) (already exists)
  # african american (Black)
  # diabetes (Diabetic)

# Lets build the data table now with all the different arguments in here
  allData_train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
    filter(!sample %in% good_metaf$initial) %>% 
    select(SRNlesion, fit_positive, BMIcat, great50, Hx_of_Polyps, Hx_Prev, Hx_Fam_CRC, Black, Diabetic, 
           one_of(fitGroup_SRNlesion_confirmed_vars[, 'otus'])) %>% na.omit()
  
  
UPdated_SRNcatLesion_table <- makeCat(allData_train, cutoff = 0, ignorecols = 9)
Updated_onePercent_SRNcatLesion_table <- makeCat(allData_train, cutoff = abund_cutoff, ignorecols = 9)
UPdated_SRNcatLesion_table$crc_bugs <- apply(UPdated_SRNcatLesion_table[, knownCRC_vars], 1, function(x) getCount(x))

#Replace specific high counts with greater or less than 1%
UPdated_SRNcatLesion_table$Otu000003 <- Updated_onePercent_SRNcatLesion_table$Otu000003
UPdated_SRNcatLesion_table$Otu000005 <- Updated_onePercent_SRNcatLesion_table$Otu000005
UPdated_SRNcatLesion_table$Otu000008 <- Updated_onePercent_SRNcatLesion_table$Otu000008

#Replace CRC OTUs with actual sequence counts
UPdated_SRNcatLesion_table$Otu000566 <- allData_train$Otu000566
UPdated_SRNcatLesion_table$Otu000126 <- allData_train$Otu000126
UPdated_SRNcatLesion_table$Otu000205 <- allData_train$Otu000205
UPdated_SRNcatLesion_table$Otu000397 <- allData_train$Otu000397
UPdated_SRNcatLesion_table$Otu000563 <- allData_train$Otu000563

set.seed(050416)
test_SRNlesionModel_wmeta <- AUCRF(SRNlesion~., data=UPdated_SRNcatLesion_table, pdel=0.05, ntree=500, ranking='MDA')

## Can I get a little better if I select specific measures (No I cannot, OOB error is 20.24% versus 18.55% - difference might be due to noise)

set.seed(050416)
Boruta_SRNlesionModel_wmeta <- Boruta(SRNlesion~., data=UPdated_SRNcatLesion_table, mcAdj=TRUE, maxRuns=1000)


test_SRNlesion_confirmed_vars <- as.data.frame(Boruta_SRNlesionModel_wmeta['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

# Use the selected data set in AUCRF now
updated_SRNlesion_selected_train <- select(UPdated_SRNcatLesion_table, 
                                            SRNlesion, one_of(test_SRNlesion_confirmed_vars[, 'otus']))

set.seed(050416)
selected_test_SRNlesionModel_wmeta <- AUCRF(SRNlesion~., data=updated_SRNlesion_selected_train, pdel=0.99, ntree=500, ranking='MDA')




# What happens if instead of the harder problem of SRN lesion I apply this to lesion (normal vs. adenoma + cancer)


# First, look at simply using positive/negative as a cutoff
catLesion_table <- makeCat(boruta_selected_list_data[["lesion"]], cutoff = 0)

set.seed(050416)
catLesion_table_AUCRF_model<- AUCRF(lesion~., data=catLesion_table, pdel=0.99, ntree=500, ranking='MDA')
# Throws out 14, 147, and 166 from equation
# Lets convert these back to numeric

catLesion_table$Otu000014 <- boruta_selected_list_data[["lesion"]]$Otu000014
catLesion_table$Otu000147 <- boruta_selected_list_data[["lesion"]]$Otu000147
catLesion_table$Otu000166 <- boruta_selected_list_data[["lesion"]]$Otu000166

#try again to see if that improved the classification
set.seed(050416)
part_catLesion_table_AUCRF_model<- AUCRF(lesion~., data=catLesion_table, pdel=0.99, ntree=500, ranking='MDA')

#try a 1% cutoff
onePercent_catLesion_table <- makeCat(boruta_selected_list_data[["lesion"]], cutoff = abund_cutoff)

catLesion_table$Otu000014 <- onePercent_catLesion_table$Otu000014

set.seed(050416)
onePercent_catLesion_table_AUCRF_model<- AUCRF(lesion~., data=catLesion_table, pdel=0.99, ntree=500, ranking='MDA')
# Convert back since the latter has a better OOB-AUCopt


# Perhaps adding number of known associated OTUs with CRC will help prediction


knownCRC_vars <- c("Otu000566", "Otu000126", "Otu000205", "Otu000397", "Otu000563")

catLesion_table$crc_bugs <- apply(SRNcatLesion_table[, knownCRC_vars], 1, function(x) getCount(x))
set.seed(050416)
onePercent_Part_catLesion_table_AUCRF_model<- AUCRF(lesion~., data=catLesion_table, pdel=0.99, ntree=500, ranking='MDA')

# First convert crc_bugs back to sequence counts

catLesion_table$Otu000566 <- boruta_selected_list_data[["lesion"]]$Otu000566
catLesion_table$Otu000126 <- boruta_selected_list_data[["lesion"]]$Otu000126


set.seed(050416)
onePercent_Part_catLesion_table_AUCRF_model<- AUCRF(lesion~., data=catLesion_table, pdel=0.99, ntree=500, ranking='MDA')


# Lets build the data table now with all the different arguments in here
allData_lesion_train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  mutate(crc_bugs = apply(.[, knownCRC_vars], 1, function(x) getCount(x))) %>% 
  select(lesion, fit_positive, BMIcat, great50, Hx_of_Polyps, Hx_Prev, Hx_Fam_CRC, Black, Diabetic, crc_bugs, 
         one_of(fitGroup_lesion_confirmed_vars[, 'otus'])) %>% na.omit()


UPdated_Lesion_table <- makeCat(allData_lesion_train, cutoff = 0, ignorecols = 10)
UPdated_onepercent_Lesion_table <- makeCat(allData_lesion_train, cutoff = abund_cutoff, ignorecols = 10)

#Replace specific high counts with greater or less than 1%
UPdated_Lesion_table$Otu000014 <- UPdated_onepercent_Lesion_table$Otu000014


#Replace CRC OTUs with actual sequence counts
UPdated_Lesion_table$Otu000566 <- allData_lesion_train$Otu000566
UPdated_Lesion_table$Otu000126 <- allData_lesion_train$Otu000126


set.seed(050416)
test_lesionModel_wmeta <- AUCRF(lesion~., data=UPdated_Lesion_table, pdel=0.05, ntree=500, ranking='MDA')

set.seed(050416)
Boruta_lesionModel_wmeta <- Boruta(lesion~., data=UPdated_Lesion_table, mcAdj=TRUE, maxRuns=1000)


test_lesion_confirmed_vars <- as.data.frame(Boruta_lesionModel_wmeta['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

# Use the selected data set in AUCRF now
updated_lesion_selected_train <- select(UPdated_Lesion_table, 
                                           lesion, one_of(test_lesion_confirmed_vars[, 'otus']))

set.seed(050416)
selected_test_lesionModel_wmeta <- AUCRF(lesion~., data=updated_lesion_selected_train, pdel=0.99, ntree=500, ranking='MDA')

## Can I get a little better if I select specific measures (No I cannot, OOB error is 20.24% versus 20% - difference might be due to noise)
## But very low false negative rate (22/245 or 0.08979592 ~ 9%)


## Need to remove history of polyps, cancer, and family crc to see how it performs without this information
## Closer to real life pre screening

# Lets build the data table now with all the different arguments in here
noHX_allData_train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% mutate(crc_bugs = apply(.[, knownCRC_vars], 1, function(x) getCount(x))) %>% 
  select(SRNlesion, fit_positive, BMIcat, great50, Black, Diabetic, crc_bugs, 
         one_of(fitGroup_SRNlesion_confirmed_vars[, 'otus'])) %>% na.omit()


noHx_SRNcat_table <- makeCat(noHX_allData_train, cutoff = 0, ignorecols = 7)
noHx_onePercent_SRNcat_table <- makeCat(noHX_allData_train, cutoff = abund_cutoff, ignorecols = 7)

#Replace specific high counts with greater or less than 1%
noHx_SRNcat_table$Otu000003 <- noHx_onePercent_SRNcat_table$Otu000003
noHx_SRNcat_table$Otu000005 <- noHx_onePercent_SRNcat_table$Otu000005
noHx_SRNcat_table$Otu000008 <- noHx_onePercent_SRNcat_table$Otu000008

#Replace CRC OTUs with actual sequence counts
noHx_SRNcat_table$Otu000566 <- noHX_allData_train$Otu000566
noHx_SRNcat_table$Otu000126 <- noHX_allData_train$Otu000126
noHx_SRNcat_table$Otu000205 <- noHX_allData_train$Otu000205
noHx_SRNcat_table$Otu000397 <- noHX_allData_train$Otu000397
noHx_SRNcat_table$Otu000563 <- noHX_allData_train$Otu000563

set.seed(050416)
noHx_SRNlesionModel_wmeta <- AUCRF(SRNlesion~., data=noHx_SRNcat_table, pdel=0.05, ntree=500, ranking='MDA')
    # OOB error is 21.24% about a 1% increase and an 0.02 decrease in AUC 

# Let's try lesion now

noHx_allData_lesion_train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  mutate(crc_bugs = apply(.[, knownCRC_vars], 1, function(x) getCount(x))) %>% 
  select(lesion, fit_positive, BMIcat, great50, Black, Diabetic, crc_bugs, 
         one_of(fitGroup_lesion_confirmed_vars[, 'otus'])) %>% na.omit()


noHx_Lesion_table <- makeCat(noHx_allData_lesion_train, cutoff = 0, ignorecols = 7)
noHx_onepercent_Lesion_table <- makeCat(noHx_allData_lesion_train, cutoff = abund_cutoff, ignorecols = 7)

#Replace specific high counts with greater or less than 1%
noHx_Lesion_table$Otu000014 <- noHx_onepercent_Lesion_table$Otu000014


#Replace CRC OTUs with actual sequence counts
noHx_Lesion_table$Otu000566 <- noHx_onepercent_Lesion_table$Otu000566
noHx_Lesion_table$Otu000126 <- noHx_onepercent_Lesion_table$Otu000126

set.seed(050416)
noHx_lesionModel_wmeta <- AUCRF(lesion~., data=noHx_Lesion_table, pdel=0.05, ntree=500, ranking='MDA')
  # OOB estimate of error becomes 26.49%

# So something about SRN lesion is actually descriminative versus keeping all adenomas together



