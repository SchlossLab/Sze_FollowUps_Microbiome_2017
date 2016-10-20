### Test the first Random Forest without Fit
### Base analysis without metadata or OTU groupings
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))


# Read in data and remove unneeded tables and lists
load("exploratory/RFwoFit.RData")
rm(selected_train, sens_specif_table, corr_pvalue_ROC_table, orig_probs, impfactor_Data_List, orig_RF_run, 
   selected_probs, selected_RF_run, variableList, modelList)

# Get threshold used by each of the different models to make the call

cutoffs <- as.data.frame(unlist(lapply(rocNameList[c("SRNlesion_train_roc", "SRNlesion_selected_train_roc", 
                                                     "lesion_train_roc", "lesion_selected_train_roc")], 
                                       function(y) coords(y, x='best', ret='thr')))) %>% 
  cbind(rownames(.), .) %>%  mutate(model = c(rep("SRNlesion", 2), rep("lesion", 2))) %>% 
  mutate(dataset = rep(c("All", "Select"), length(rownames(.))/2))
#take the threshold point used for the best sensitivity and specificty

colnames(cutoffs) <- c("L1", "cutoff", "model", "dataset")
cutoffs$model <- factor(cutoffs$model, levels = c("SRNlesion", "lesion"))

write.csv(cutoffs, "results/tables/noFIT.cutoffs.csv", row.names = F)

initial <- inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
  select(SRNlesion, lesion, fit_result, contains("Otu0"))


good_metaf$lesionf[good_metaf$Disease_Free == 'n'] <- 1
good_metaf$lesionf[good_metaf$Disease_Free == 'y' | good_metaf$Disease_Free == 'unknown'] <- 0

followups <- inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
  select(lesionf, fit_followUp, contains("Otu0")) %>% rename(SRNlesion = lesionf) %>% rename(fit_result = fit_followUp)

followups <- cbind(good_metaf$lesionf, followups)
colnames(followups)[1] <- c("lesion")

rf_opt_NameList <- list(SRNlesion_rf_opt = orig_rf_opt[["SRNlesion"]], 
                        SRNlesion_selected_rf_opt = selected_rf_opt[["SRNlesion"]], 
                        lesion_rf_opt= orig_rf_opt[["lesion"]], lesion_selected_rf_opt = selected_rf_opt[["lesion"]])

initial_predictions <- lapply(rf_opt_NameList, function(x) predict(x, initial, type='prob'))

followup_predictions <- lapply(rf_opt_NameList, function(x) predict(x, followups, type='prob'))


#Create ROC data for lesion model

AUCData <- read.csv("results/tables/testSetAUC.csv", header = T, stringsAsFactors = F)

lesion_combined_ROC <- roc(c(initial$lesion, followups$lesion) ~ 
                             c(initial_predictions$lesion_rf_opt[, 2], followup_predictions$lesion_rf_opt[, 2]))

write.csv(rbind(AUCData, c("woFit", lesion_combined_ROC$auc)), "results/tables/testSetAUC.csv", row.names = F)

variableList <- c("sensitivities", "specificities")
sens_specif_table <- makeSensSpecTable(list(lesion = lesion_combined_ROC), variableList, modelList = "woFitTest")

write.csv(sens_specif_table, "results/tables/woFit_test_ROCCurve_sens_spec.csv")



# Join good_metaf table with metaF table
combined_metaData <- inner_join(good_metaf, metaF, by = "EDRN")


df_initial_neg <- melt(initial_predictions) %>% filter(Var2 == 0)
df_initial_pos <- melt(initial_predictions) %>%  filter(Var2 == 1)

df_initial_preds <- cbind(rename(select(df_initial_neg, value), negative = value), 
                          rename(select(df_initial_pos, value, L1), positive = value))

rm(df_initial_neg, df_initial_pos)

df_initial_preds$model[grep("lesion", df_initial_preds$L1)] <- "lesion"
df_initial_preds$model[grep("SRN", df_initial_preds$L1)] <- "SRNlesion"

df_initial_preds$dataset[grep("rf", df_initial_preds$L1)] <- "All"
df_initial_preds$dataset[grep("selected", df_initial_preds$L1)] <- "Select"
df_initial_preds <- mutate(df_initial_preds, 
                           diseaseFree = c(rep("n", 67*4))) %>% mutate(diagnosis = rep(good_metaf$Diagnosis, 4))

df_initial_preds$time_point <- "initial"

df_followups_neg <- melt(followup_predictions) %>% filter(Var2 == 0)
df_followups_pos <- melt(followup_predictions) %>% filter(Var2 == 1)

df_followups_preds <- cbind(rename(select(df_followups_neg, value), negative = value), 
                            rename(select(df_followups_pos, value, L1), positive = value))

rm(df_followups_neg, df_followups_pos)

df_followups_preds$model[grep("lesion", df_followups_preds$L1)] <- "lesion"
df_followups_preds$model[grep("SRN", df_followups_preds$L1)] <- "SRNlesion"


df_followups_preds$dataset[grep("rf", df_followups_preds$L1)] <- "All"
df_followups_preds$dataset[grep("selected", df_followups_preds$L1)] <- "Select"
df_followups_preds <- mutate(df_followups_preds, diseaseFree = rep(good_metaf$Disease_Free, 4)) %>% 
  mutate(diagnosis = rep(good_metaf$Diagnosis, 4))

df_followups_preds$time_point <- "followup"

df_InitFollow_ALL <- rbind(df_initial_preds, df_followups_preds)
df_InitFollow_ALL$model <- factor(df_InitFollow_ALL$model, levels = c("SRNlesion", "lesion"))

df_InitFollow_ALL <- mutate(df_InitFollow_ALL, 
                            detailed_diagnosis = 
                              rep(combined_metaData$Dx_Bin.y, 
                                  length(rownames(df_InitFollow_ALL))/length(rownames(combined_metaData)))) %>% 
  mutate(EDRN = rep(good_metaf$EDRN, 8))

write.csv(df_InitFollow_ALL, "results/tables/noFIT.models.datatable.csv", row.names = F)

Names_facet <- c('SRNlesion' = "Classify SRN + Cancer", 'lesion' = "Classify Lesion")


# Graph the Adenoma data only

adn_graph <- grid.arrange(
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

ggsave(file = "results/figures/adn_predictions_noFit.tiff", adn_graph, width=8, height = 8, dpi = 300)


# Graph the Cancer data only

crc_graph <- grid.arrange(
  
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


ggsave(file = "results/figures/crc_predictions_noFit.tiff", crc_graph, width=8, height = 8, dpi = 300)


rm(list=setdiff(ls(), "lesion_combined_ROC"))
save.image("exploratory/woFit_Test_ROC.RData")









