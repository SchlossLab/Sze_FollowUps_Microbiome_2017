### Figure S4
### How do the full models of lesion and IF look
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

# Load IF model data
IF_follow_up_probability <- read.csv(
  "results/tables/IF_follow_up_probability_summary.csv", 
  header = T, stringsAsFactors = F)

# Read in meta data tables
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", 
                       stringsAsFactors = F, header = T) %>% 
  filter(!is.na(fit_followUp))


accuracy_plot <- grid.arrange(
  # Graph the lesion carcinoma data only
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
      legend.position = c(0.25, 0.15), 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the lesion adenoma data only
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
    ggtitle("C") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = c(0.20, 0.15), 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the IF model CRC data only
  filter(IF_follow_up_probability, diagnosis != "adenoma") %>% 
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
  filter(IF_follow_up_probability, diagnosis == "adenoma") %>%
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

ggsave(file = "results/figures/FigureS4.pdf", accuracy_plot, 
       width=8.5, height = 11, dpi = 300)










