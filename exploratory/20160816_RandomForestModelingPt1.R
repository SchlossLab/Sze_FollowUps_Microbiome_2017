## Run Previous Analysis align follow up data
## Code used and modified from Niel .Rmd file
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", "wesanderson"))

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
lesion_confirmed_vars <- as.data.frame(cancer_impFactorData['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

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

# Use the selected data set in AUCRF now
SRNlesion_selected_train <- select(train, SRNlesion, one_of(SRNlesion_confirmed_vars[, 'otus']))
set.seed(050416)
SRNlesion_selected_rf_opt <- AUCRF(SRNlesion~., data=SRNlesion_selected_train, pdel=0.99, ntree=500, ranking='MDA')$RFopt
SRNlesion_selected_train_probs <- predict(SRNlesion_selected_rf_opt, type='prob')[,2]

SRNlesion_selected_train_roc <- roc(train$SRNlesion ~ SRNlesion_selected_train_probs)


### Graph the ROC curves for each of the different models and test for difference

# Created needed vectors and lists
rocNameList <- list(cancer_train_roc = cancer_train_roc, cancer_selected_train_roc = cancer_selected_train_roc, 
                    SRNlesion_train_roc = SRNlesion_train_roc, SRNlesion_selected_train_roc = SRNlesion_selected_train_roc, 
                    lesion_train_roc = lesion_train_roc, lesion_selected_train_roc = lesion_selected_train_roc)

variableList <- c("sensitivities", "specificities")
modelList <- c("cancerALL", "cancerSELECT", "SRNlesionALL", "SRNlesionSELECT", "lesionALL", "lesionSELECT")

sens_specif_table <- makeSensSpecTable(rocNameList, variableList, modelList)

# Create the graph
ggplot(sens_specif_table, aes(sensitivities, specificities)) + 
  geom_line(aes(group = model, color = model), size = 1.5) + scale_x_continuous(trans = "reverse") + 
  theme_bw() + xlab("Sensitivity") + ylab("Specificity") + 
  theme(axis.title = element_text(face = "bold"))

# Run the statistic test

test <- unname(unlist(roc.test(cancer_train_roc, cancer_selected_train_roc)['p.value']))





# Get OTU abundances and calculate probabilities for initial samples that have follow ups from full data set
initial <- inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
  select(lesion, fit_result, contains("Otu0"))

initial_probs <- predict(rf_opt, initial, type='prob')[,2] 

# create first new data column
good_metaf$cancerPos[good_metaf$Disease_Free == 'n'] <- "y"
good_metaf$cancerPos[good_metaf$Disease_Free == 'y'] <- "n"
good_metaf$cancerPos[good_metaf$Disease_Free == 'unknown'] <- "u"
good_metaf$cancerPos <- factor(good_metaf$cancerPos, levels = c("n", "y", "u"))

# create second new data column
good_metaf$lesionf[good_metaf$Disease_Free == 'n'] <- 1
good_metaf$lesionf[good_metaf$Disease_Free == 'y' | good_metaf$Disease_Free == 'unknown'] <- 0

followups <- inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
  select(lesionf, fit_result, contains("Otu0")) %>% rename(lesion = lesionf)

follow_probs <- predict(rf_opt, followups, type='prob')[,2]

# Create a ROC curve for full data set
train_roc <- roc(train$lesion~train_probs)
cutoff <- coords(train_roc, x='best', ret='thr')
    #take the threshold point used for the best sensitivity and specificty

# Get a table and graph from the results
followupsUpdate <- inner_join(good_metaf, shared, by = c("followUp" = "Group"))

probs_table <- as.data.frame(cbind(initial_probs, follow_probs)) %>% mutate(cancer = initial$lesion) %>% 
  mutate(cancerF = factor(followupsUpdate$cancerPos, levels = c("n", "y", "u")))


### Generate classification model with Boruta feature selection

# ID important factors using Boruta Independent of previous publications
set.seed(050416)
impFactorData <- Boruta(lesion~., data=train, mcAdj=TRUE, maxRuns=1000)
  # Does not change after increasing runs to 2000

# Get the confirmed important variables
imp_confirmed_otus <- as.data.frame(impFactorData['finalDecision'])  %>% 
                                  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

selected_train <- select(train, lesion, one_of(imp_confirmed_otus[, 'otus']))
selected_initial <- select(initial, lesion, one_of(imp_confirmed_otus[, 'otus']))
selected_followUp <- select(followups, lesion, one_of(imp_confirmed_otus[, 'otus']))

selected_probs_table <- as.data.frame(cbind(initial_probs, follow_probs)) %>% mutate(cancer = initial$lesion) %>% 
  mutate(cancerF = factor(followupsUpdate$cancerPos, levels = c("n", "y", "u")))


#Generate AUCRF model and prediction results on this selected data
set.seed(050416)
selected_rf_opt <- AUCRF(lesion~., data=selected_train, pdel=0.05, ntree=500, ranking='MDA')$RFopt
selected_train_probs <- predict(rf_opt, type='prob')[,2]

selected_initial_probs <- predict(selected_rf_opt, selected_initial, type='prob')[,2]
selected_follow_probs <- predict(selected_rf_opt, selected_followUp, type='prob')[,2]

selected_probs_table <- as.data.frame(cbind(selected_initial_probs, selected_follow_probs)) %>% 
  mutate(cancer = selected_initial$lesion) %>% 
  mutate(cancerF = factor(followupsUpdate$cancerPos, levels = c("n", "y", "u")))


selected_train_roc <- roc(selected_train$lesion~selected_train_probs)
selected_cutoff <- coords(selected_train_roc, x='best', ret='thr')
    #take the threshold point used for the best sensitivity and specificty


### Graph output of model for all features and selected features 
### Trainset is based on excluding samples with follow up data

grid.arrange(
  # Graph the initial data with cancer diagnosis uses all features (A)
  ggplot(probs_table, aes(cancer, initial_probs)) + geom_jitter(aes(color = cancer), width = 0.25, size = 4) + theme_bw() + 
    scale_color_manual(values = wes_palette("GrandBudapest"), labels = c("No", "Yes"))  + ylim(0, 1) + 
    geom_hline(yintercept = cutoff, size = 2, linetype = 2) + scale_x_discrete(labels = c("No", "Yes")) + 
    xlab("Initial Cancer Diagnosis") + ylab("Probability") + labs(colour = "Cancer at\nStool\nSampling") + 
    ggtitle("A)    Initial Samples with All Features") + 
    theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"), 
          plot.title = element_text(hjust = 0)), 
  
  # Graph the follow up data with cancer diagnosis uses all features (B)
  ggplot(probs_table, aes(cancer, follow_probs)) + geom_jitter(aes(color = cancerF), width = 0.25, size = 4) + theme_bw() + 
    scale_color_manual(values = wes_palette("GrandBudapest"), labels = c("No", "Yes"))  + ylim(0, 1) + 
    geom_hline(yintercept = cutoff, size = 2, linetype = 2) + labs(colour = "Cancer at\nStool\nSampling") + 
    scale_x_discrete(labels = c("No", "Yes")) + ylab("Probability") + xlab("Initial Cancer Diagnosis") + 
    ggtitle("B)    Follow Up Samples with All Features") + 
    theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"), 
          plot.title = element_text(hjust = 0)), 
  
  # Graph the initial data with cancer diagnosis selected features only (c)
  ggplot(selected_probs_table, aes(cancer, selected_initial_probs)) + 
    geom_jitter(aes(color = cancer), width = 0.25, size = 4) + theme_bw() + 
    scale_color_manual(values = wes_palette("GrandBudapest"), labels = c("No", "Yes"))  + ylim(0, 1) + 
    labs(colour = "Cancer at\nStool\nSampling") + xlab("Initial Cancer Diagnosis") + 
    geom_hline(yintercept = selected_cutoff, size = 2, linetype = 2) + scale_x_discrete(labels = c("No", "Yes")) + 
    ylab("Probability") + ggtitle("C)    Initial Samples with Selected Features") + 
    theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"), 
          plot.title = element_text(hjust = 0)),  
  
  # Graph the follow up data with cancer diagnosis selected features only (D)
  ggplot(selected_probs_table, aes(cancer, selected_follow_probs)) + 
    geom_jitter(aes(color = cancerF), width = 0.25, size = 4) + theme_bw() + 
    scale_color_manual(values = wes_palette("GrandBudapest"), labels = c("No", "Yes", "Unknown"))  + ylim(0, 1) + 
    geom_hline(yintercept = selected_cutoff, size = 2, linetype = 2) + 
    labs(colour = "Cancer at\nStool\nSampling") + scale_x_discrete(labels = c("No", "Yes")) + 
    ylab("Probability") + xlab("Initial Cancer Diagnosis") + 
    ggtitle("D)    Follow Up Samples with Selected Features") + 
    theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"), 
          plot.title = element_text(hjust = 0))
  )


###Pull out the 4 important OTUs from previous manuscript and make a ggplot format table 
###for those with follow up###

otus_int <- filter(imp_confirmed_otus, otus != "fit_result")[, 'otus']
selected_taxa <- tax[otus_int, ] %>% select(Taxonomy)

taxa_int <- c('Otu000034' = "Collinsella", 
              'Otu000059' = "Clostridiales(1)", 
              'Otu000062' = "Odoribacter", 
              'Otu000084' = "Bilophila", 
              'Otu000103' = "Lachnospiraceae", 
              'Otu000126' = "P.asaccharolytica",
              'Otu000145' = "Enterorhabdus", 
              'Otu000147' = "Prevotella", 
              'Otu000205' = "F.nucleatum", 
              'Otu000217' = "Coriobacteriaceae", 
              'Otu000333' = "Fusobacterium", 
              'Otu000381' = "Firmicutes", 
              'Otu000397' = "P. micra", 
              'Otu000563' = "Porphyromonas", 
              'Otu000566' = "P.stomatis",
              'Otu000671' = "Dialister", 
              'Otu000902' = "Clostridiales(2)")



#delta_index <- c(rep("delta_Otu000126", 2), rep("delta_Otu000566", 2), rep("delta_Otu000205", 2), 
#                rep("delta_Otu000397", 2), rep("delta_fit", 2))

jitter_values <- jitter(rep(0, 67), amount = 0.08)

cancer_positive <- as.data.frame(cbind(rep(as.character(good_metaf$cancerPos), 36), rep(good_metaf$EDRN, 36))) %>% 
  rename(cancerPos = V1, EDRN = V2)


#creating the data tables to be used
select_otu_data <- melt(
  rbind(select(initial, lesion, one_of(otus_int)) %>% 
          mutate(sampleType = "initial", EDRN = good_metaf$EDRN, Diagnosis = good_metaf$Diagnosis), 
        select(followups, lesion, one_of(otus_int)) %>% 
          mutate(sampleType = "follow_up", EDRN = good_metaf$EDRN, Diagnosis = good_metaf$Diagnosis)), 
  id=c("sampleType", "Diagnosis", "EDRN")) %>% filter(variable != "lesion") %>% 
  mutate(jitteredValues = as.numeric(value) + rep(jitter_values, 36)) %>% cbind(select(cancer_positive, cancerPos))


fit_results <- select(good_metaf, EDRN, Diagnosis, fit_result, fit_followUp) %>% 
  rename(initial = fit_result, follow_up = fit_followUp) %>% melt(id=c("EDRN", "Diagnosis")) %>% 
  mutate(jitteredValues = as.numeric(value) + rep(jitter_values, 2)) %>% mutate(cancerPos = rep(good_metaf$cancerPos, 2))

average_change <- select(good_metaf, EDRN, Diagnosis) %>% 
  mutate(delta_Otu000126 = followups$Otu000126 - initial$Otu000126, 
         delta_Otu000566 = followups$Otu000566 - initial$Otu000566, 
         delta_Otu000205 = followups$Otu000205 - initial$Otu000205, 
         delta_Otu000397 = followups$Otu000397 - initial$Otu000397, 
         delta_fit = good_metaf$fit_followUp - good_metaf$fit_result) %>% 
  melt(id=c("EDRN", "Diagnosis")) %>% filter(value != "NA")
average_change$reduc_value[which(average_change$value >= 0)] <- "no"
average_change$reduc_value[which(average_change$value < 0)] <- "yes"

adenoma_count <- rbind(
  filter(average_change, Diagnosis == "adenoma", variable == "delta_Otu000126") %>% count(reduc_value), 
  filter(average_change, Diagnosis == "adenoma", variable == "delta_Otu000566") %>% count(reduc_value), 
  filter(average_change, Diagnosis == "adenoma", variable == "delta_Otu000205") %>% count(reduc_value), 
  filter(average_change, Diagnosis == "adenoma", variable == "delta_Otu000397") %>% count(reduc_value), 
  filter(average_change, Diagnosis == "adenoma", variable == "delta_fit") %>% count(reduc_value)) %>% 
  mutate(delta_index) %>% mutate(fraction = 
           n/c(rep(length(rownames(metaFAonly)), 8), rep(length(rownames(
             filter(average_change, Diagnosis == "adenoma", variable == "delta_fit"))), 2)))

cancer_count <- rbind(
  filter(average_change, Diagnosis == "adenocarcinoma", variable == "delta_Otu000126") %>% count(reduc_value), 
  filter(average_change, Diagnosis == "adenocarcinoma", variable == "delta_Otu000566") %>% count(reduc_value), 
  filter(average_change, Diagnosis == "adenocarcinoma", variable == "delta_Otu000205") %>% count(reduc_value), 
  filter(average_change, Diagnosis == "adenocarcinoma", variable == "delta_Otu000397") %>% count(reduc_value), 
  filter(average_change, Diagnosis == "adenocarcinoma", variable == "delta_fit") %>% count(reduc_value)) %>% 
  mutate(delta_index) %>% mutate(fraction = 
            n/(length(rownames(filter(average_change, Diagnosis == "adenocarcinoma", variable == "delta_Otu000397")))))



###Create graph of the data

grid.arrange(
  # Make Adenoma OTU part of graph
  filter(select_otu_data, Diagnosis == "adenoma") %>% 
    ggplot(aes(factor(sampleType, levels = c("initial", "follow_up")), 
               log10(jitteredValues + 1.1), group = factor(EDRN))) + 
    geom_line(aes(color = factor(cancerPos))) + geom_point(aes(color = factor(cancerPos))) +  
    facet_wrap(~variable, labeller = as_labeller(taxa_int)) + scale_y_continuous(limits = c(0, 3)) + 
    theme_bw() + ylab("Log Total Sequences") + xlab("") + ggtitle("Adenoma") + 
    scale_colour_manual(values = "blue") +  
    theme(legend.position="none", plot.title = element_text(size=20, face="bold"), 
          strip.text.x = element_text(size = 5), axis.text.x = element_text(size = 6)), 
  
  # Make Cancer OTU part of graph
  filter(select_otu_data, Diagnosis == "adenocarcinoma") %>% 
    ggplot(aes(factor(sampleType, levels = c("initial", "follow_up")), 
               log10(jitteredValues+1.1), group = factor(EDRN))) + 
    geom_line(aes(color = factor(cancerPos, levels = c("n", "y", "u")))) + 
    geom_point(aes(color = factor(cancerPos, levels = c("n", "y", "u")))) + 
    facet_wrap(~variable, labeller = as_labeller(taxa_int)) + scale_y_continuous(limits = c(0, 3)) + 
    theme_bw() + ylab("Log Total Sequences") + xlab("") + ggtitle("Cancer") + 
    scale_colour_manual(values = c("pink", "darkred", "red")) + 
    theme(legend.position="none", plot.title = element_text(size=20, face="bold"), 
          strip.text.x = element_text(size = 5), axis.text.x = element_text(size = 6)), 
  
  # Make Adenoma FIT part of graph
  filter(fit_results, Diagnosis == "adenoma") %>% 
    ggplot(aes(factor(variable, levels = c("initial", "follow_up")), 
               log10(jitteredValues+1.1), group = factor(EDRN))) + 
    geom_line(aes(color = factor(cancerPos))) + geom_point(aes(color = factor(cancerPos))) + 
    theme_bw() + ylab("Log FIT Result") + xlab("") + scale_y_continuous(limits = c(0, 3.5)) + 
    scale_colour_manual(values = "blue") + theme(legend.position="none"), 
  
  # Make Cancer OTU part of graph
  filter(fit_results, Diagnosis == "adenocarcinoma") %>% 
    ggplot(aes(factor(variable, levels = c("initial", "follow_up")), 
               log10(jitteredValues+1.1), group = factor(EDRN))) + 
    geom_line(aes(color = factor(cancerPos, levels = c("n", "y", "u")))) + 
    geom_point(aes(color = factor(cancerPos, levels = c("n", "y", "u")))) + 
    theme_bw() + ylab("Log FIT Result") + xlab("") + scale_y_continuous(limits = c(0, 3.5)) + 
    scale_colour_manual(values = c("pink", "darkred", "red")) + theme(legend.position="none")
  
  # Make bar plot of adenomas with reduction
#  ggplot(adenoma_count, aes(factor(delta_index), fraction, fill = reduc_value)) + geom_bar(stat="identity") + 
#    scale_y_continuous(labels=percent_format()) + theme_bw() + ylab("Percent of Total") + xlab("") + 
#    scale_fill_manual(values = c("darkblue", "lightblue"), name = "Reduction\nObserved") + labs(color = "Reduction Observed") + 
#    scale_x_discrete(labels = c("FIT", "P.\nasaccharolytica", "P.\nstomatis", "F.\nnucleatum", "P.\nmicra")) + 
#    theme(axis.text.x = element_text(size = 5)), 
  
  # Make a bar plot of cancers with reduction
#  ggplot(cancer_count, aes(factor(delta_index), fraction, fill = reduc_value)) + geom_bar(stat="identity") + 
#    scale_y_continuous(labels=percent_format()) + theme_bw() + ylab("Percent of Total") + xlab("") + 
#    scale_fill_manual(values = c("darkred", "pink"), name = "Reduction\nObserved") + labs(color = "Reduction Observed") + 
#    scale_x_discrete(labels = c("FIT", "P.\nasaccharolytica", "P.\nstomatis", "F.\nnucleatum", "P.\nmicra")) + 
#    theme(axis.text.x = element_text(size = 5))
  
)










