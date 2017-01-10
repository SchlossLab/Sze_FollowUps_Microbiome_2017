### Prediction of Follow Ups Analysis
### How do the with and without FIT models do in correctly calling initial and follow up samples
### P-value table of comparison and Figure 4 graph
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "caret"))

# Load lesion model data
follow_up_probability <- read.csv(
  "results/tables/follow_up_probability_summary.csv", 
  header = T, stringsAsFactors = F)

red_follow_up_probability <- read.csv("results/tables/reduced_follow_up_probability_summary.csv", 
                                      stringsAsFactors = F, header = T)

# Load IF model data
IF_follow_up_probability <- read.csv(
  "results/tables/IF_follow_up_probability_summary.csv", 
  header = T, stringsAsFactors = F)

IF_red_follow_up_probability <- read.csv("results/tables/reduced_IF_follow_up_probability_summary.csv", 
                                      stringsAsFactors = F, header = T)


# Read in meta data tables
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", 
                       stringsAsFactors = F, header = T) %>% 
filter(!is.na(fit_followUp))


# create tables to hold wilcoxson paired tests pvalues with BH correction
lesion_wilcox_pvalue_summary <- getProb_PairedWilcox(follow_up_probability)
lesion_red_wilcox_pvalue_summary <- getProb_PairedWilcox(red_follow_up_probability)
IF_wilcox_pvalue_summary <- getProb_PairedWilcox(IF_follow_up_probability)
IF_red_wilcox_pvalue_summary <- getProb_PairedWilcox(IF_red_follow_up_probability)

# Create a confusion matrix
general_summary <- get_confusion_data(follow_up_probability, good_metaf)

adenoma_summary <- get_confusion_data(filter(follow_up_probability, Dx_Bin != "cancer"), 
                                      filter(good_metaf, Dx_Bin != "cancer"))

crc_summary <- get_confusion_data(filter(follow_up_probability, Dx_Bin == "cancer"), 
                                      filter(good_metaf, Dx_Bin == "cancer"))


# Create confusion tables
general_I <- make_confusionTable(follow_up_probability, good_metaf)
general_F <- make_confusionTable(follow_up_probability, good_metaf, n = 67, m = 132)

# Create Figure 4 
# Visual summery of the pvalues obtained

accuracy_plot <- grid.arrange(
  # Graph the lesion carcinoma data only
  filter(red_follow_up_probability, diagnosis != "adenoma") %>% 
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
      legend.position = c(0.25, 0.15), 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the lesion adenoma data only
  filter(red_follow_up_probability, diagnosis == "adenoma") %>%
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
    ggtitle("C") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = c(0.20, 0.15), 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the IF model CRC data only
  filter(IF_red_follow_up_probability, diagnosis != "adenoma") %>% 
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
    ggtitle("B") + ylab("Postive Probability") + 
    xlab("") + theme_bw() + theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = c(0.25, 0.15), 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the IF model adenoma data only
  filter(IF_red_follow_up_probability, diagnosis == "adenoma") %>%
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
    ggtitle("D") + ylab("Postive Probability") + xlab("") + theme_bw() + 
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

write.csv(general_I, "results/tables/initial_confusion_matrix.csv")
write.csv(general_F, "results/tables/followup_confusion_matrix.csv")


