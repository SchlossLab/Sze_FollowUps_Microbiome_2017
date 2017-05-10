### Build Figure 1 
### Differences in measurements and in beta-diversity as visualized by NMDS
## Marc Sze

# Load needed functions and libraries
source('code/functions.R')

loadLibs(c("dplyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))

# Load needed data
difference_table_treatment <- read.csv("data/process/tables/difference_table.csv", header = T)
diff_summary_table <- read.csv("data/process/tables/change_theta_fit_summary.csv", header = T, row.names = 1)
thetayc_adn_IF <- read.csv("data/process/tables/thetayc_adn_IF.csv", header = T)
thetayc_srn_IF <- read.csv("data/process/tables/thetayc_srn_IF.csv", header = T)
thetayc_crc_IF <- read.csv("data/process/tables/thetayc_crc_IF.csv", header = T)
thetayc_summary_table <- read.csv("data/process/tables/beta_diver_summary.csv", header = T, row.names = 1)

# mean_cl_boot:
# very fast implementation of the basic nonparametric bootstrap for obtaining confidence limits 
# for the population mean without assuming normality.

#Difference between initial and follow up for thetayc and fit broken down by adenoma and cancer
differences_graph <- grid.arrange(
  # Difference from bacterial community structure between adenoma and cancer
  ggplot(difference_table_treatment, 
         aes(factor(Dx_Bin, levels = c("adenoma", "adv_adenoma", "cancer"), 
                    labels = c("Adenoma", "Advanced\nAdenoma", "Carcinoma")), distance, group = 1)) + 
    geom_jitter(aes(color=Dx_Bin), width = 0.3, size = 4, alpha = 0.5) + 
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
                 colour = "black", geom = "crossbar", size = 0.5, width = 0.5) + 
    scale_color_manual(name = "Lesion Type", 
                       values = c('#228B22', '#FFD700', '#DC143C'), 
                       breaks = c("adenoma", "adv_adenoma", "cancer"), 
                       labels = c("Adenoma", "Advanced\nAdenoma", "Carcinoma")) + 
    coord_cartesian(ylim = c(0, 1)) + ylab(expression(paste(theta[YC], " Distance"))) + 
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), 
                       labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),limits = c(-0.06, 1.06)) + 
    xlab("") + theme_bw() + ggtitle("A") + 
    theme(axis.title = element_text(face="bold", hjust = 0.5), 
          axis.text.x = element_text(size = 12), 
          legend.title = element_text(face="bold"), 
          legend.position = "none", 
          plot.title = element_text(face="bold", hjust = -0.065, size = 20), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()), 
  
  # Adenoma Only initial and follow up
  mutate(thetayc_adn_IF, samples = factor(samples, 
                                          levels = c("initial", "follow_up"), 
                                          labels = c("Pre", "Post"))) %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
    geom_point(aes(color=samples, shape = samples), size = 3) + 
    theme_bw() + coord_equal() + ggtitle("B") + ylim(-0.7, 0.7) + xlim(-0.7, 0.7) + 
    scale_color_manual(values = c('#00EE76', '#228B22')) + scale_fill_manual(values = c('#00EE76', '#228B22')) + 
    scale_shape_manual(values = c(16, 1)) + 
    theme(plot.title = element_text(face="bold", hjust = -0.28, size = 20), 
          legend.position = c(0.15, 0.11), 
          legend.title = element_blank(), 
          legend.key = element_blank(), 
          legend.background = element_rect(color = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    coord_fixed(ratio = 1.40), 
  
  # SRN Only initial and follow up
  mutate(thetayc_srn_IF, samples = factor(samples, 
                                          levels = c( "initial", "follow_up"), 
                                          labels = c("Pre", "Post"))) %>% 
    ggplot(aes(x=NMDS1, y=NMDS2)) + 
    geom_point(aes(color=samples, shape = samples), size = 3) + 
    theme_bw() + coord_equal() + ggtitle("C") + ylim(-0.7, 0.7) + xlim(-0.7, 0.7) + 
    scale_color_manual(values = c('#CDAD00', '#8B7500')) + scale_fill_manual(values = c('#CDAD00', '#8B7500')) + 
    scale_shape_manual(values = c(16, 1)) + 
    theme(plot.title = element_text(face="bold", hjust = -0.32, size = 20), 
          legend.position = c(0.15, 0.11), 
          legend.title = element_blank(), 
          legend.key = element_blank(), 
          legend.background = element_rect(color = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    coord_fixed(ratio = 1.40), 
  
  # CRC Only initial and follow up
  mutate(thetayc_crc_IF, 
         samples = factor(samples, 
                          levels = c("initial", "follow_up"), 
                          labels = c("Pre", "Post"))) %>%
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
    geom_point(aes(color=samples, shape = samples), size = 3) + theme_bw() + 
    coord_equal() + ggtitle("D") + ylim(-0.7, 0.7) + xlim(-0.7, 0.7) +
    scale_color_manual(values = c('#EEA2AD', '#DC143C')) + scale_fill_manual(values = c('#EEA2AD', '#DC143C')) + 
    scale_shape_manual(values = c(16, 1)) + 
    theme(plot.title = element_text(face="bold", hjust = -0.32, size = 20), 
          legend.position = c(0.15, 0.11), 
          legend.title = element_blank(), 
          legend.background = element_rect(color = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    coord_fixed(ratio = 1.40), layout_matrix = cbind(c(1, 2), c(1, 3), c(1, 4)))



ggsave(file = "results/figures/Figure1.pdf", differences_graph, 
       width=10, height = 8.5, dpi = 300)


