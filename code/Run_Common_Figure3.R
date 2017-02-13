### Create Figure 3
### Graph to show AUC of full data with chosen model
### Graph to show how min, middle, and max performed on 20% test data
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", 
           "wesanderson", "caret", "pROC"))

### Load in needed data tables
common_lesion_model_roc <- read.csv("results/tables/common_test_data_roc.csv", header = T)
common_lesion_auc_info <- read.csv("results/tables/common_auc_summary.csv", 
  header = T, row.names = 1)

# Load common follow up data
common_follow_up_probability <- read.csv("results/tables/common_follow_up_probability_summary.csv", 
                                      stringsAsFactors = F, header = T)
# Read in meta data tables
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", 
                       stringsAsFactors = F, header = T) %>% 
  filter(!is.na(fit_followUp))


# create tables to hold wilcoxson paired tests pvalues with BH correction
common_wilcox_pvalue_summary <- getProb_PairedWilcox(common_follow_up_probability)

write.csv(common_wilcox_pvalue_summary, 
          "results/tables/common_model_wilcox_paired_pvalue_summary.csv", row.names = F)



common_graph <- grid.arrange(
  # Graph showing all variable and reduced variable lesion model
  filter(common_lesion_model_roc, run != "full_roc") %>% 
    ggplot(aes(x = sensitivities, y = specificities)) + 
    geom_polygon(alpha = 0.25) + 
    geom_line(aes(group = run), size = 1.25, color = "black", alpha = 0.5) + 
    geom_line(data = filter(common_lesion_model_roc, run == "full_roc"), 
              size = 1.5, color = "cornflowerblue") + 
    scale_x_continuous(trans = "reverse") + theme_bw() + geom_abline(intercept = 1, linetype = 2, size = 1) + 
    ggtitle("A") + xlab("Sensitivity") + ylab("Specificity") + theme(plot.title = element_text(face= "bold")), 
  
  # Graph the common carcinoma data only
  filter(common_follow_up_probability, diagnosis != "adenoma") %>% 
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
    ggtitle("B") + ylab("Lesion Postive Probability") + 
    xlab("") + theme_bw() + theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = c(0.25, 0.15), 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the common adenoma data only
  filter(common_follow_up_probability, diagnosis == "adenoma") %>%
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
    ggtitle("C") + ylab("Lesion Postive Probability") + xlab("") + theme_bw() + 
    theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = c(0.20, 0.12), 
      plot.title = element_text(face="bold", hjust = 0)))
  
  

# Save graph image as a pdf
ggsave(file = "results/figures/Figure3.tiff", roc_curve_graph, device = "tiff", 
       width=8, height = 8, dpi = 300)



