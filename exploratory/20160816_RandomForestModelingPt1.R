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



### Need to amend and separate Adenoma and CRC
good_metaf <- read.csv('data/process/followUp_outcome_data.csv', header = T, 
                       stringsAsFactors = F) %>% inner_join(metaF, by="EDRN")

metaFConly <- filter(good_metaf, Diagnosis == "adenocarcinoma" | Diagnosis == "N/D")
  
metaFAonly <- filter(good_metaf, Diagnosis == "adenoma")


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

# Using the seperation as Cancer (ALL)
train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  select(cancer, fit_result, contains("Otu0")) %>% na.omit()

set.seed(050416)
cancer_rf_opt <- AUCRF(cancer~., data=train, pdel=0.05, ntree=500, ranking='MDA')$RFopt
cancer_train_probs <- predict(cancer_rf_opt, type='prob')[,2]

cancer_train_roc <- roc(train$cancer ~ cancer_train_probs)

# ID important factors for Cancer using Boruta 
set.seed(050416)
cancer_impFactorData <- Boruta(cancer~., data=train, mcAdj=TRUE, maxRuns=1000)
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
set.seed(050416)
lesion_impFactorData <- Boruta(lesion~., data=train, mcAdj=TRUE, maxRuns=1000)
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
set.seed(050416)
SRNlesion_impFactorData <- Boruta(SRNlesion~., data=train, mcAdj=TRUE, maxRuns=1000)
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
set.seed(050416)
threeway_impFactorData <- Boruta(threeway~., data=train, mcAdj=TRUE, maxRuns=1000)
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

# Obtain the pvalue statistics as well as the bonferroni corrected values

corr_pvalue_ROC_table <- getROCPvalue(rocNameList, modelList, 8, multi = T)

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
write.csv(df_InitFollow_ALL, "results/tables/withFIT.models.datatable.csv", row.names = F)


# Create labels for subset of data on graph

Names_facet <- c('threeGroups' = "Classify N, A, C", 'cancer' = "Classify Cancer", 
                 'SRNlesion' = "Classify SRN + Cancer", 'lesion' = "Classify Lesion")

# Graph the ALL data only

grid.arrange(
  # Graph the adenoma ALL data only
  filter(df_InitFollow_ALL, diagnosis == "adenoma" & dataset == "All") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive)) + 
    geom_jitter(aes(color=diseaseFree), width = 0.3) + 
    scale_color_manual(name = "Adenoma\nFree", values = wes_palette("GrandBudapest")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + ylim(0, 1) + 
    geom_hline(data = filter(cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Adenomas (Train on All Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")), 
  
  # Graph the cancer ALL data only
  filter(df_InitFollow_ALL, diagnosis != "adenoma" & dataset == "All") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive)) + 
    geom_jitter(aes(color=factor(diseaseFree, levels = c("n", "y", "unknown"))), width = 0.3) + 
    scale_color_manual(name = "Cancer\nFree", values = wes_palette("GrandBudapest")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + ylim(0, 1) + 
    geom_hline(data = filter(cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Cancer (Train on All Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")),
  
  # Graph adenoma select data only
  filter(df_InitFollow_ALL, diagnosis == "adenoma" & dataset == "Select") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive)) + 
    geom_jitter(aes(color=diseaseFree), width = 0.3) + 
    scale_color_manual(name = "Adenoma\nFree", values = wes_palette("GrandBudapest")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + ylim(0, 1) + 
    geom_hline(data = filter(cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Adenomas (Train on Select Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")), 
  
  # Graph cancer select data only
  filter(df_InitFollow_ALL, diagnosis != "adenoma" & dataset == "Select") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive)) + 
    geom_jitter(aes(color=factor(diseaseFree, levels = c("n", "y", "unknown"))), width = 0.3) + 
    scale_color_manual(name = "Cancer\nFree", values = wes_palette("GrandBudapest")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + ylim(0, 1) + 
    geom_hline(data = filter(cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
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




