## Use Random Forest but run with only microbiome
## Focus strictly on Lesion, SRNLesion, and three groups (Normal, Adenoma, Cancer) classifications
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

metaF <- read.delim('data/process/followUps_metadata.txt', header=T, sep='\t') %>% mutate(lesion = factor(NA, levels=c(0,1)))

metaI <- read.delim('data/process/initials_metadata.tsv', header=T, sep='\t') %>% mutate(lesion = factor(NA, levels=c(0,1)))

good_metaf <- read.csv('data/process/followUp_outcome_data.csv', header = T, 
                       stringsAsFactors = F) %>% inner_join(metaF, by="EDRN")

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
    # This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components
rm(tax)


### Run First organization to sort lesion vs. non-lesion
metaF$lesion[metaF$dx =='normal'] <- 0
metaF$lesion[metaF$dx =='adenoma' | metaF$dx == 'cancer'] <- 1
metaF$lesion <- factor(metaF$lesion)

metaI$lesion[metaI$dx == 'normal'] <- 0
metaI$lesion[metaI$dx == 'adenoma' | metaI$dx == 'cancer'] <- 1
metaI$lesion <- factor(metaI$lesion)


### Run Second organization to sort SRN (Adv Adenoma) to cancer group and Adenoma to non-cancer group
metaF$SRNlesion[metaF$Dx_Bin =='normal' | metaF$Dx_Bin == 'Adenoma'] <- 0
metaF$SRNlesion[metaF$Dx_Bin =='adv Adenoma' | metaF$Dx_Bin == 'Cancer'] <- 1
metaF$SRNlesion <- factor(metaF$SRNlesion)

metaI$SRNlesion[metaI$Dx_Bin =='High Risk Normal' | metaI$Dx_Bin == 'Normal' | metaI$Dx_Bin == 'Adenoma'] <- 0
metaI$SRNlesion[metaI$Dx_Bin =='adv Adenoma' | metaI$Dx_Bin == 'Cancer'] <- 1
metaI$SRNlesion <- factor(metaI$SRNlesion)


### Run Third organization to sort into three separate groups: normal, adenoma, cancer
metaF$threeGroup[metaF$dx == 'normal'] <- 0
metaF$threeGroup[metaF$dx == 'adenoma'] <- 1
metaF$threeGroup[metaF$dx == 'cancer'] <- 2
metaF$threeGroup <- factor(metaF$threeGroup)

metaI$threeGroup[metaI$dx == 'normal'] <- 0
metaI$threeGroup[metaI$dx == 'adenoma'] <- 1
metaI$threeGroup[metaI$dx == 'cancer'] <- 2
metaI$threeGroup <- factor(metaI$threeGroup)


# Using the seperation as Lesion
train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  select(lesion, contains("Otu0")) %>% na.omit()

set.seed(050416)
lesion_rf_opt_OTUs <- AUCRF(lesion~., data=train, pdel=0.05, ntree=500, ranking='MDA')$RFopt
lesion_train_probs_OTUs <- predict(lesion_rf_opt_OTUs, type='prob')[,2]

lesion_train_roc_OTUs <- roc(train$lesion ~ lesion_train_probs_OTUs)

# ID important factors for lesion using Boruta 
lesion_rf_model_OTUs_used <- lesion_rf_opt_OTUs$importance

selected_train <- select(train, lesion, one_of(rownames(lesion_rf_model_OTUs_used)))

set.seed(050416)
lesion_impFactorData_OTUs <- Boruta(lesion~., data=selected_train, mcAdj=TRUE, maxRuns=1000)
    # Does not change after increasing runs to 2000

# Get the confirmed important variables
lesion_confirmed_vars_OTUs <- as.data.frame(lesion_impFactorData_OTUs['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

# write table to results to be used by other analysis components
write.csv(lesion_confirmed_vars_OTUs, "results/tables/lesion_confirmed_vars_NOFIT.csv", row.names = F)

# Use the selected data set in AUCRF now
lesion_selected_train_OTUs <- select(train, lesion, one_of(lesion_confirmed_vars_OTUs[, 'otus']))
set.seed(050416)
lesion_selected_rf_opt_OTUs <- AUCRF(lesion~., data=lesion_selected_train_OTUs, pdel=0.99, ntree=500, ranking='MDA')$RFopt
lesion_selected_train_probs_OTUs <- predict(lesion_selected_rf_opt_OTUs, type='prob')[,2]

lesion_selected_train_roc_OTUs <- roc(train$lesion ~ lesion_selected_train_probs_OTUs)

# Using the seperation as SRNlesion
train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  select(SRNlesion, contains("Otu0")) %>% na.omit()

set.seed(050416)
SRNlesion_rf_opt_OTUs <- AUCRF(SRNlesion~., data=train, pdel=0.05, ntree=500, ranking='MDA')$RFopt
SRNlesion_train_probs_OTUs <- predict(SRNlesion_rf_opt_OTUs, type='prob')[,2]

SRNlesion_train_roc_OTUs <- roc(train$SRNlesion ~ SRNlesion_train_probs_OTUs)

# ID important factors for lesion using Boruta 
SRNlesion_rf_model_OTUs_used <- SRNlesion_rf_opt_OTUs$importance

selected_train <- select(train, SRNlesion, one_of(rownames(SRNlesion_rf_model_OTUs_used)))

set.seed(050416)
SRNlesion_impFactorData_OTUs <- Boruta(SRNlesion~., data=selected_train, mcAdj=TRUE, maxRuns=1000)
  # Does not change after increasing runs to 2000

# Get the confirmed important variables
SRNlesion_confirmed_vars_OTUs <- as.data.frame(SRNlesion_impFactorData_OTUs['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

# write table to results to be used by other analysis components
write.csv(SRNlesion_confirmed_vars_OTUs, "results/tables/SRNlesion_confirmed_vars_NOFIT.csv", row.names = F)

# Use the selected data set in AUCRF now
SRNlesion_selected_train_OTUs <- select(train, SRNlesion, one_of(SRNlesion_confirmed_vars_OTUs[, 'otus']))
set.seed(050416)
SRNlesion_selected_rf_opt_OTUs <- AUCRF(SRNlesion~., 
                                        data=SRNlesion_selected_train_OTUs, pdel=0.99, ntree=500, ranking='MDA')$RFopt
SRNlesion_selected_train_probs_OTUs <- predict(SRNlesion_selected_rf_opt_OTUs, type='prob')[,2]

SRNlesion_selected_train_roc_OTUs <- roc(train$SRNlesion ~ SRNlesion_selected_train_probs_OTUs)


# Using the seperation as threeGroup
train <- inner_join(metaI, shared, by = c("sample" = "Group")) %>% 
  filter(!sample %in% good_metaf$initial) %>% 
  select(threeGroup, contains("Otu0")) %>% na.omit()

set.seed(050416)
threeGroup_rf_opt_OTUs <- AUCRF(threeGroup~., data=train, pdel=0.05, ntree=500, ranking='MDA')$RFopt
threeGroup_train_probs_OTUs <- predict(threeGroup_rf_opt_OTUs, type='prob')[,2]

threeGroup_train_roc_OTUs <- roc(train$threeGroup ~ threeGroup_train_probs_OTUs)

# ID important factors for lesion using Boruta 
threeGroup_rf_model_OTUs_used <- threeGroup_rf_opt_OTUs$importance

selected_train <- select(train, threeGroup, one_of(rownames(threeGroup_rf_model_OTUs_used)))

set.seed(050416)
threeGroup_impFactorData_OTUs <- Boruta(threeGroup~., data=selected_train, mcAdj=TRUE, maxRuns=1000)
  # Does not change after increasing runs to 2000

# Get the confirmed important variables
threeGroup_confirmed_vars_OTUs <- as.data.frame(threeGroup_impFactorData_OTUs['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)

# write table to results to be used by other analysis components
write.csv(threeGroup_confirmed_vars_OTUs, "results/tables/threeGroups_confirmed_vars_NOFIT.csv", row.names = F)

# Use the selected data set in AUCRF now
threeGroup_selected_train_OTUs <- select(train, threeGroup, one_of(threeGroup_confirmed_vars_OTUs[, 'otus']))
set.seed(050416)
threeGroup_selected_rf_opt_OTUs <- AUCRF(
  threeGroup~., data=threeGroup_selected_train_OTUs, pdel=0.99, ntree=500, ranking='MDA')$RFopt
threeGroup_selected_train_probs_OTUs <- predict(threeGroup_selected_rf_opt_OTUs, type='prob')[,2]

threeGroup_selected_train_roc_OTUs <- roc(train$threeGroup ~ threeGroup_selected_train_probs_OTUs)


### Graph the ROC curves for each of the different models and test for difference

# Created needed vectors and lists
rocNameList_OTUs <- list(threeGroup_train_roc = threeGroup_train_roc_OTUs, 
                         threeGroup_selected_train_roc = threeGroup_selected_train_roc_OTUs, 
                    SRNlesion_train_roc = SRNlesion_train_roc_OTUs, 
                    SRNlesion_selected_train_roc = SRNlesion_selected_train_roc_OTUs, 
                    lesion_train_roc = lesion_train_roc_OTUs, lesion_selected_train_roc = lesion_selected_train_roc_OTUs)

variableList <- c("sensitivities", "specificities")
modelList <- c("ThreeGroupALL", "ThreeGroupSELECT", "SRNlesionALL", "SRNlesionSELECT", "lesionALL", "lesionSELECT")

sens_specif_table <- makeSensSpecTable(rocNameList_OTUs, variableList, modelList)
write.csv(sens_specif_table, "results/tables/ROCCurve_sens_spec_NOFIT.csv")

# Obtain the pvalue statistics as well as the bonferroni corrected values

corr_pvalue_ROC_table <- getROCPvalue(rocNameList_OTUs, modelList, 6, multi = T)
write.csv(corr_pvalue_ROC_table, "results/tables/ROCCurve_NOFIT_corrected_pvalues.csv")


# Create the graph
ggplot(sens_specif_table, aes(sensitivities, specificities)) + 
  geom_line(aes(group = model, color = model), size = 1.5) + scale_x_continuous(trans = "reverse") + 
  theme_bw() + xlab("Sensitivity") + ylab("Specificity") + 
  theme(axis.title = element_text(face = "bold"))


### Investigate how these models do with initial and follow up samples

# Get threshold used by each of the different models to make the call

cutoffs <- as.data.frame(unlist(lapply(rocNameList_OTUs, function(y) coords(y, x='best', ret='thr')))) %>% 
  cbind(rownames(.), .) %>%  mutate(model = c(rep("threeGroups", 2), rep("SRNlesion", 2), rep("lesion", 2))) %>% 
  mutate(dataset = rep(c("All", "Select"), length(rownames(.))/2))
    #take the threshold point used for the best sensitivity and specificty

colnames(cutoffs) <- c("L1", "cutoff", "model", "dataset")
cutoffs$model <- factor(cutoffs$model, levels = c("threeGroups", "SRNlesion", "lesion"))

write.csv(cutoffs, "results/tables/noFIT.cutoffs.csv", row.names = F)

# Join good_metaf table with metaF table
combined_metaData <- inner_join(good_metaf, metaF, by = "EDRN")

# Create data frames to be used for initial and follow up samples

initial <- inner_join(metaF, shared, by = c("initial" = "Group")) %>% 
  select(threeGroup, lesion, SRNlesion, contains("Otu0"))

good_metaf$lesionf[good_metaf$Disease_Free == 'n'] <- 1
good_metaf$lesionf[good_metaf$Disease_Free == 'y' | good_metaf$Disease_Free == 'unknown'] <- 0

followups <- inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
  select(lesionf, contains("Otu0")) %>% rename(SRNlesion = lesionf)

followups <- cbind(good_metaf$lesionf, good_metaf$lesionf, followups)
colnames(followups)[1:2] <- c("threeGroup", "lesion")


# Get the prediction tables for each group using initial and follow up data

rf_opt_NameList_OTUs <- list(threeGroup_rf_opt = threeGroup_rf_opt_OTUs, 
                             threeGroup_selected_rf_opt = threeGroup_selected_rf_opt_OTUs, 
                        SRNlesion_rf_opt = SRNlesion_rf_opt_OTUs, SRNlesion_selected_rf_opt = SRNlesion_selected_rf_opt_OTUs, 
                        lesion_rf_opt= lesion_rf_opt_OTUs, lesion_selected_rf_opt = lesion_selected_rf_opt_OTUs)

initial_predictions <- lapply(rf_opt_NameList_OTUs, function(x) predict(x, initial, type='prob'))

followup_predictions <- lapply(rf_opt_NameList_OTUs, function(x) predict(x, followups, type='prob'))

df_initial_neg <- melt(initial_predictions) %>% filter(Var2 == 0)
df_initial_pos <- melt(initial_predictions) %>%  filter(Var2 == 1)

df_initial_preds <- cbind(rename(select(df_initial_neg, value), negative = value), 
                          rename(select(df_initial_pos, value, L1), positive = value))

rm(df_initial_neg, df_initial_pos)

df_initial_preds$model[grep("lesion", df_initial_preds$L1)] <- "lesion"
df_initial_preds$model[grep("SRN", df_initial_preds$L1)] <- "SRNlesion"
df_initial_preds$model[grep("three", df_initial_preds$L1)] <- "threeGroups"

df_initial_preds$dataset[grep("rf", df_initial_preds$L1)] <- "All"
df_initial_preds$dataset[grep("selected", df_initial_preds$L1)] <- "Select"
df_initial_preds <- mutate(df_initial_preds, 
                           diseaseFree = c(rep("n", 67*6))) %>% mutate(diagnosis = rep(good_metaf$Diagnosis, 6))

df_initial_preds$time_point <- "initial"

df_followups_neg <- melt(followup_predictions) %>% filter(Var2 == 0)
df_followups_pos <- melt(followup_predictions) %>% filter(Var2 == 1)

df_followups_preds <- cbind(rename(select(df_followups_neg, value), negative = value), 
                            rename(select(df_followups_pos, value, L1), positive = value))

rm(df_followups_neg, df_followups_pos)

df_followups_preds$model[grep("lesion", df_followups_preds$L1)] <- "lesion"
df_followups_preds$model[grep("SRN", df_followups_preds$L1)] <- "SRNlesion"
df_followups_preds$model[grep("three", df_followups_preds$L1)] <- "threeGroups"

df_followups_preds$dataset[grep("rf", df_followups_preds$L1)] <- "All"
df_followups_preds$dataset[grep("selected", df_followups_preds$L1)] <- "Select"
df_followups_preds <- mutate(df_followups_preds, diseaseFree = rep(good_metaf$Disease_Free, 6)) %>% 
  mutate(diagnosis = rep(good_metaf$Diagnosis, 6))

df_followups_preds$time_point <- "followup"

df_InitFollow_ALL <- rbind(df_initial_preds, df_followups_preds)
df_InitFollow_ALL$model <- factor(df_InitFollow_ALL$model, levels = c("threeGroups", "SRNlesion", "lesion"))

df_InitFollow_ALL <- mutate(df_InitFollow_ALL, 
                            detailed_diagnosis = rep(combined_metaData$Dx_Bin.y, 
                                                     length(rownames(df_InitFollow_ALL))/length(rownames(combined_metaData))))

write.csv(df_InitFollow_ALL, "results/tables/noFIT.models.datatable.csv", row.names=F)

# Create labels for subset of data on graph

Names_facet <- c('threeGroups' = "Classify N, A, C", 'SRNlesion' = "Classify SRN + Cancer", 'lesion' = "Classify Lesion")


# Graph the ALL data only

grid.arrange(
  # Graph the adenoma ALL data only
  filter(df_InitFollow_ALL, diagnosis == "adenoma" & dataset == "All") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive)) + 
    geom_jitter(aes(color=detailed_diagnosis), width = 0.3) + 
    scale_color_manual(name = "Polyp Type", values = c("cyan", "blue"), 
                       breaks = c("Adenoma", "adv Adenoma"), labels = c("Adenoma", "SRN")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Adenomas (Train on All Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")), 
  
  # Graph the cancer ALL data only
  filter(df_InitFollow_ALL, diagnosis != "adenoma" & dataset == "All") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive)) + 
    geom_jitter(aes(color=factor(diseaseFree, levels = c("n", "y", "unknown"))), width = 0.3) + 
    scale_color_manual(name = "Cancer\nFree", values = wes_palette("GrandBudapest")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Cancer (Train on All Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")),
  
  # Graph adenoma select data only
  filter(df_InitFollow_ALL, diagnosis == "adenoma" & dataset == "Select") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive)) + 
    geom_jitter(aes(color=detailed_diagnosis), width = 0.3) + 
    scale_color_manual(name = "Polyp Type", values = c("cyan", "blue"), 
                       breaks = c("Adenoma", "adv Adenoma"), labels = c("Adenoma", "SRN")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Adenomas (Train on Select Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")), 
  
  # Graph cancer select data only
  filter(df_InitFollow_ALL, diagnosis != "adenoma" & dataset == "Select") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive)) + 
    geom_jitter(aes(color=factor(diseaseFree, levels = c("n", "y", "unknown"))), width = 0.3) + 
    scale_color_manual(name = "Cancer\nFree", values = wes_palette("GrandBudapest")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Cancer (Train on Select Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold"))
)


### Create a new group with normals and follow ups without cancer: find OTUs that are still different
followUp_Select <- rbind(filter(metaI, dx == "normal") %>% select(sample, dx) %>% mutate(dx = as.character(dx)) %>% 
                           mutate(disease = rep(0, length(rownames(.)))), 
                         select(metaF, followUp) %>% mutate(dx = rep("follow_up", 
                                                                     length(rownames(metaF)))) %>% rename(sample = followUp) %>% 
                           mutate(disease = rep(1, length(rownames(metaF))))) %>% 
  inner_join(shared, by = c("sample" = "Group")) %>% select(disease, contains("Otu0")) %>% mutate(disease = factor(disease))


set.seed(050416)
follow_v_norm_diffs <- Boruta(disease~., data=followUp_Select, mcAdj=TRUE, maxRuns=1000)
follow_v_norm_diffs_rf_opt <- AUCRF(disease~., data=followUp_Select, pdel=0.05, ntree=500, ranking='MDA')$RFopt

# Get the confirmed important variables
follow_v_norm_diffs_confirmed_vars_OTUs <- as.data.frame(follow_v_norm_diffs['finalDecision'])  %>% 
  mutate(otus = rownames(.))  %>% filter(finalDecision == "Confirmed")  %>% select(otus)


# create table for p-value test
pValue_test_data <- as.data.frame(cbind(c(metaI$sample[metaI$dx == "normal"], metaF$initial, metaF$followUp), 
                                                          c(rep("normal", length(rownames(metaI[metaI$dx == "normal", ]))), 
                                                            rep("initial", length(rownames(metaF))), 
                                                            rep("follow_up", length(rownames(metaF)))))) %>% 
  rename(Group = V1, sampleType = V2) %>% mutate(Group = as.integer(as.character((Group)))) %>% 
  inner_join(shared, by = "Group") %>% select(Group, sampleType, one_of(follow_v_norm_diffs_confirmed_vars_OTUs$otus)) %>% 
  mutate(Otutotal8OTUs = rowSums(.[, follow_v_norm_diffs_confirmed_vars_OTUs$otus]))


test <- getCorrectedPvalue(pValue_test_data, c(follow_v_norm_diffs_confirmed_vars_OTUs$otus, "Otutotal8OTUs"), 
                           "Group", wilcox = F, kruskal = T)

graph_data <- melt(pValue_test_data, id = c("Group", "sampleType"))

# Graph

# Convert taxa table to a data frame with columns for each taxonomic division
tax <- read.delim('data/process/followUps.final.an.unique_list.0.03.cons.taxonomy', 
                  sep='\t', header=T, row.names=1)
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
# This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components
rm(tax)

diffvNorm_selected_taxa <- tax_df[follow_v_norm_diffs_confirmed_vars_OTUs[, 'otus'], ]
diffvNorm_selected_labs <- c(createTaxaLabeller(diffvNorm_selected_taxa), Otutotal8OTUs = "Total of OTUs")


ggplot(graph_data, aes(factor(sampleType, levels = c("normal", "initial", "follow_up")), log(value+1), group = 1)) + 
  geom_jitter(aes(color=factor(sampleType, levels = c("normal", "initial", "follow_up"))), width = 0.5) + 
  scale_color_discrete(name = "Sample Type", breaks = c("normal", "initial", "follow_up"), 
                       labels = c("normal", "initial", "follow_up")) + 
  stat_summary(fun.y = mean, colour = "black", geom = "line") + ggtitle("Lesion Model") + 
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +  
  facet_wrap(~variable, labeller = as_labeller(diffvNorm_selected_labs)) + ylab("Subsampled Sequence Reads") + 
  xlab("") + theme_bw() + theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), 
                                axis.ticks.x = element_blank(), plot.title = element_text(face = "bold")) 















