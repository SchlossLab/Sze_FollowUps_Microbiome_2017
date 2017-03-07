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
thetayc_crc_IF <- read.csv("data/process/tables/thetayc_crc_IF.csv", header = T)
thetayc_summary_table <- read.csv("data/process/tables/beta_diver_summary.csv", header = T, row.names = 1)

# mean_cl_boot:
# very fast implementation of the basic nonparametric bootstrap for obtaining confidence limits 
# for the population mean without assuming normality.

#Difference between initial and follow up for thetayc and fit broken down by adenoma and cancer
differences_graph <- grid.arrange(
  # Difference from bacterial community structure between adenoma and cancer
  ggplot(difference_table_treatment, 
         aes(factor(dx, levels = c("adenoma", "cancer"), labels = c("Adenoma", "Carcinoma")), 
             distance, group = 1)) + 
    geom_jitter(aes(color=Dx_Bin), width = 0.3, size = 4, alpha = 0.5) + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + 
    scale_color_manual(name = "Lesion Type", 
                       values = c('#00FFFF', '#0000FF', '#F1BB7B'), 
                       breaks = c("adenoma", "adv_adenoma", "cancer"), 
                       labels = c("Adenoma", "SRN", "Carcinoma")) + 
    coord_cartesian(ylim = c(0, 1)) + ylab(expression(paste(theta[YC], " Distance"))) + 
    xlab("") + theme_bw() + ggtitle("A") + 
    theme(axis.title = element_text(face="bold", hjust = 0.5), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), 
          legend.title = element_text(face="bold"), 
          legend.position = c(0.75, 0.2), 
          title = element_text(face="bold", hjust = 0), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()), 
  
  # Difference from fit between adenoma and cancer
  ggplot(difference_table_treatment, 
         aes(factor(dx, levels = c("adenoma", "cancer"), labels = c("Adenoma", "Carcinoma")), 
             fit_difference*-1, group = 1)) + 
    geom_jitter(aes(color=Dx_Bin), width = 0.3, size = 4, alpha = 0.5) + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + 
    scale_color_manual(name = "Lesion Type", 
                       values = c('#00FFFF', '#0000FF', '#F1BB7B'), 
                       breaks = c("adenoma", "adv_adenoma", "cancer"), 
                       labels = c("Adenoma", "SRN", "Carcinoma")) + 
    ylab("Change in FIT from Follow Up to Initial") + 
    xlab("") + theme_bw() + ggtitle("B") + 
    theme(axis.title = element_text(face="bold", hjust = 0.5), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), 
          legend.title = element_text(face="bold"), 
          legend.position = "none", 
          title = element_text(face="bold", hjust = 0), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()), 
  
  # Adenoma Only initial and follow up
  mutate(thetayc_adn_IF, samples = factor(samples, levels = c("follow_up", "initial"), labels = c("Follow Up", "Initial"))) %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) + geom_point(aes(color=samples)) + 
    theme_bw() + coord_equal() + ggtitle("C") + ylim(-0.7, 0.7) + xlim(-0.7, 0.7) + 
    stat_ellipse(aes(group = samples, color = samples, fill = samples), 
                 alpha = 0.25, geom = "polygon", show.legend = FALSE) + 
    scale_color_manual(values = c('#483D8B', '#00BFFF')) + scale_fill_manual(values = c('#483D8B', '#00BFFF')) +
    theme(plot.title = element_text(face = "bold", hjust = 0), 
          legend.position = c(0.20, 0.12), 
          legend.title = element_blank(), 
          legend.key = element_blank(), 
          legend.background = element_rect(color = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    coord_fixed(ratio = 1.40), 
  
  mutate(thetayc_crc_IF, 
         samples = factor(samples, 
                          levels = c("follow_up", "initial"), 
                          labels = c("Follow Up", "Initial"))) %>%
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
    geom_point(aes(color=samples)) + theme_bw() + 
    coord_equal() + ggtitle("D") + ylim(-0.7, 0.7) + xlim(-0.7, 0.7) +
    stat_ellipse(aes(group = samples, color = samples, fill = samples), 
                 alpha = 0.25, geom = "polygon", show.legend = FALSE) + 
    scale_color_manual(values = c('#FD6467', '#F1BB7B')) + scale_fill_manual(values = c('#FD6467', '#F1BB7B')) + 
    theme(plot.title = element_text(face = "bold", hjust = 0), 
          legend.position = c(0.20, 0.12), 
          legend.title = element_blank(), 
          legend.background = element_rect(color = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    coord_fixed(ratio = 1.4), ncol = 2, nrow = 2)


ggsave(file = "results/figures/Figure1.pdf", differences_graph, 
       width=8.5, height = 9, dpi = 300)
















