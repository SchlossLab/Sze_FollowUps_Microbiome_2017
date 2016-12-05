### Prediction of Follow Ups Analysis
### How do the with and without FIT models do in correctly calling initial and follow up samples
### P-value table of comparison and Figure 4 graph
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "caret"))


follow_up_probability <- read.csv(
  "results/tables/follow_up_probability_summary.csv", 
  header = T, stringsAsFactors = F)

# Read in data tables
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", 
                       stringsAsFactors = F, header = T) %>% 
filter(!is.na(fit_followUp))


# create table to hold wilcoxson paired tests pvalues

wilcox_pvalue_summary <- as.data.frame(matrix(
  nrow = 4, ncol = 1, dimnames = list(
    rows = c("lesion", "all_adenoma", "carcinoma_only", "SRN_only"), 
    cols = "Pvalue")))

# Set up variable vector
lesion_type <- c("all", "cancer", "adenoma", "cancer")
filter_diagnosis <- c("none", "none", "adv_adenoma", "adenoma")

for(i in 1:length(lesion_type)){
  
  wilcox_pvalue_summary[i, "Pvalue"] <- wilcox.test(
    filter(follow_up_probability, sampleType == "initial" & 
             Dx_Bin != paste(lesion_type[i]) & 
             Dx_Bin != paste(filter_diagnosis[i]))[, "Yes"], 
    filter(follow_up_probability, sampleType == "followup" & 
             Dx_Bin != paste(lesion_type[i]) & 
             Dx_Bin != paste(filter_diagnosis[i]))[, "Yes"], 
    paired = TRUE)$p.value
  
}


# Add Benjamini-Hochberg correction
wilcox_pvalue_summary <- cbind(
  wilcox_pvalue_summary, 
  BH_correction = p.adjust(wilcox_pvalue_summary$Pvalue, 
    method = "BH")) 

# Create a confusion matrix
follow_up_probability <- mutate(
  follow_up_probability, 
  predict_call = factor(ifelse(Yes > 0.5, "Yes", "No"))) %>% 
  mutate(
    initial_call = factor(c(rep("Yes", length(rownames(good_metaf))), 
      ifelse(good_metaf$Disease_Free == "n", "Yes", "No"))))

confusion_initial <- confusionMatrix(
  follow_up_probability$predict_call[1:length(rownames(good_metaf))], 
  follow_up_probability$initial_call[1:length(rownames(good_metaf))], 
  positive = "Yes")

confusion_follow <- confusionMatrix(
  follow_up_probability$predict_call[(length(rownames(good_metaf))+1):
  length(rownames(follow_up_probability))], 
  follow_up_probability$initial_call[(length(rownames(good_metaf))+1):
  length(rownames(follow_up_probability))], positive = "Yes")

confusion_summary <- cbind(
  initial = c(confusion_initial$overall, confusion_initial$byClass), 
  followup = c(confusion_follow$overall, confusion_follow$byClass))
    #McNemar's P-value gives information on whether the 
    #prediction is significantly different than the actual



c_initial_table <- matrix(confusion_initial$table, nrow = 2, ncol = 2, 
               dimnames = list(nrow = c("pred_no", "pred_yes"), 
                ncol = c("ref_no", "ref_yes")))

c_follow_table <- matrix(confusion_follow$table, nrow = 2, ncol = 2, 
                          dimnames = list(nrow = c("pred_no", "pred_yes"), 
                            ncol = c("ref_no", "ref_yes")))


# Create Figure 4 
# Visual summery of the pvalues obtained

accuracy_plot <- grid.arrange(
  # Graph the carcinoma wFIT data only
  filter(follow_up_probability, diagnosis != "adenoma") %>% 
    ggplot(aes(factor(sampleType, 
      levels = c("initial", "followup")), Yes, group = factor(EDRN))) + 
    geom_point(aes(color=factor(disease_free, 
      levels = c("n", "y", "unknown"))), size = 2) + 
    geom_line(linetype = 2, alpha = 0.3) + 
    scale_color_manual(name = "Cancer Free", 
                       label = c("No", "Yes", "Unknown"),  
                       values = wes_palette("GrandBudapest")) + 
    scale_x_discrete(breaks = c("initial", "followup"), 
                     labels = c("Initial", "Follow Up")) + 
    coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(aes(yintercept = 0.5), linetype = 2) + 
    ggtitle("A") + ylab("Postive Probability") + 
    xlab("") + theme_bw() + theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = c(0.20, 0.15), 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the adenoma wFIT data only
  filter(follow_up_probability, diagnosis == "adenoma") %>%
    ggplot(aes(factor(sampleType, levels = c("initial", "followup")), 
               Yes, group = factor(EDRN))) + 
    geom_point(aes(color=factor(Dx_Bin))) + 
    geom_line(aes(color = factor(Dx_Bin))) + 
    scale_color_manual(name = "Polyp Type", values = c("cyan", "blue"), 
                       breaks = c("adenoma", "adv_adenoma"), 
                       labels = c("Adenoma", "SRN")) + 
    scale_x_discrete(
      breaks = c("initial", "followup"), 
      labels = c("Initial", "Follow Up")) + 
    coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(aes(yintercept = 0.5), linetype = 2) + 
    ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = c(0.20, 0.15), 
      plot.title = element_text(face="bold", hjust = 0))
  
)


# Save figures and write necessary tables

ggsave(file = "results/figures/Figure4.pdf", accuracy_plot, 
       width=8.5, height = 11, dpi = 300)

write.csv(wilcox_pvalue_summary, 
          "results/tables/IF_model_wicox_paired_pvalue_summary.csv")

write.csv(confusion_summary, "results/tables/confusion_summary.csv")

write.csv(c_initial_table, "results/tables/initial_confusion_matrix.csv")
write.csv(c_follow_table, "results/tables/followup_confusion_matrix.csv")


