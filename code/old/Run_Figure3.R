### Create Figure 2 graph
### Show results for initial and follow up samples based on testing
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))

### Load in needed data tables


adn_red_follow_up_probability <- read.csv("data/process/tables/adn_reduced_follow_up_probability_summary.csv", 
                                      stringsAsFactors = F, header = T)
srn_red_follow_up_probability <- read.csv("data/process/tables/srn_reduced_follow_up_probability_summary.csv", 
                                         stringsAsFactors = F, header = T)
crc_red_follow_up_probability <- read.csv("data/process/tables/crc_reduced_follow_up_probability_summary.csv", 
                                          stringsAsFactors = F, header = T) %>% 
  mutate(followup_crc = ifelse(disease_free == "n", "Yes", "No"))

# Create Figure
Lesion_plot <- grid.arrange(
  
  # Graph the adenoma data only
  ggplot(adn_red_follow_up_probability, 
         aes(factor(sampleType, 
                    levels = c("initial", "followup"), labels = c("Pre", "Post")), 
             Yes, group = factor(EDRN))) + 
    geom_point(color = '#006400', size = 2) + 
    geom_line(color = '#66CD00') + 
    coord_cartesian(ylim = c(0, 0.75)) + 
    geom_hline(aes(yintercept = 0.5), linetype = 2) + 
    ggtitle("A") + ylab("Adenoma Postive Probability") + xlab("") + theme_bw() + 
    theme(
      axis.title = element_text(face="bold"), 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the SRN data only
  ggplot(srn_red_follow_up_probability, 
         aes(factor(sampleType, 
                    levels = c("initial", "followup"), labels = c("Pre", "Post")), 
             Yes, group = factor(EDRN))) + 
    geom_point(color = '#8B7500', size = 2) + 
    geom_line(color = '#FFD700') + 
    coord_cartesian(ylim = c(0, 0.75)) + 
    geom_hline(aes(yintercept = 0.5), linetype = 2) + 
    ggtitle("B") + ylab("SRN Postive Probability") + xlab("") + theme_bw() + 
    theme(
      axis.title = element_text(face="bold"), 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the CRC data only
  ggplot(crc_red_follow_up_probability, 
         aes(factor(sampleType, 
                    levels = c("initial", "followup"), labels = c("Pre", "Post")), 
             Yes, group = factor(EDRN))) + 
    geom_line(color = '#CD1076') + 
    geom_point(aes(color=factor(followup_crc, 
                                levels = c("No", "Yes"))), size = 2) + 
    scale_color_manual(name = "Cancer on\nFollow Up", 
                       label = c("No", "Yes"),  
                       values = c('#CD919E', '#DC143C')) + 
    coord_cartesian(ylim = c(0, 0.75)) + 
    geom_hline(aes(yintercept = 0.5), linetype = 2) + 
    ggtitle("C") + ylab("Carcinoma Postive Probability") + xlab("") + theme_bw() + 
    theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = "none", 
      plot.title = element_text(face="bold", hjust = 0)), nrow = 1)
  
# Save figures and write necessary tables
ggsave(file = "results/figures/Figure2.tiff", Lesion_plot, 
       width=13, height = 10, dpi = 300)







