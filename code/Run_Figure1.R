### Build Figure 1 
### Differences in measurements and in beta-diversity as visualized by NMDS
## Marc Sze

# Load needed functions and libraries
source('code/functions.R')

loadLibs(c("dplyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))

# Load needed data
difference_table_treatment <- read.csv("results/tables/difference_table.csv", header = T)
diff_summary_table <- read.csv("results/tables/change_theta_fit_summary.csv", header = T, row.names = 1)
thetayc_adn_IF <- read.csv("results/tables/thetayc_adn_IF.csv", header = T)
thetayc_crc_IF <- read.csv("results/tables/thetayc_crc_IF.csv", header = T)
thetayc_summary_table <- read.csv("results/tables/beta_diver_summary.csv", header = T, row.names = 1)

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
                       values = wes_palette("GrandBudapest"), 
                       breaks = c("adenoma", "adv_adenoma", "cancer"), 
                       labels = c("Adenoma", "SRN", "Carcinoma")) + 
    coord_cartesian(ylim = c(0, 1)) + ylab("Thetayc Distance") + 
    xlab("") + theme_bw() + ggtitle("A") + 
    theme(axis.title = element_text(face="bold", hjust = 0.5), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), 
          legend.title = element_text(face="bold"), 
          legend.position = c(0.75, 0.2), 
          title = element_text(face="bold", hjust = 0)) + 
    annotate("text", label = paste("P-value = ", 
              round(diff_summary_table["thetayc", "Pvalue"], digits = 4)), 
                x = 0.85, y = 0, size = 2.5), 
  
  # Difference from fit between adenoma and cancer
  ggplot(difference_table_treatment, 
         aes(factor(dx, levels = c("adenoma", "cancer"), labels = c("Adenoma", "Carcinoma")), 
             fit_difference*-1, group = 1)) + 
    geom_jitter(aes(color=Dx_Bin), width = 0.3, size = 4, alpha = 0.5) + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + 
    scale_color_manual(name = "Lesion Type", 
                       values = wes_palette("GrandBudapest"), 
                       breaks = c("adenoma", "adv_adenoma", "cancer"), 
                       labels = c("Adenoma", "SRN", "Carcinoma")) + 
    ylab("Change in Fit from Follow up to Initial") + 
    xlab("") + theme_bw() + ggtitle("B") + 
    theme(axis.title = element_text(face="bold", hjust = 0.5), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), 
          legend.title = element_text(face="bold"), 
          legend.position = "none", 
          title = element_text(face="bold", hjust = 0)) + 
    annotate("text", label = paste("P-value = ", 
            format(round(diff_summary_table["fit", "Pvalue"], digits = 7))), 
              x = 0.90, y = -2250, size = 2.5), 
  
  # Adenoma Only initial and follow up
  ggplot(data = thetayc_adn_IF, 
         aes(x=NMDS1, y=NMDS2)) + geom_point(aes(color=samples)) + 
    theme_bw() + coord_equal() + ggtitle("C") + 
    stat_ellipse(aes(group = samples, color = samples, fill = samples), 
                 alpha = 0.25, geom = "polygon") + 
    theme(plot.title = element_text(face = "bold", hjust = 0), 
          legend.position = "none") + 
    annotate("text", label = paste("PERMANOVA = ", 
                                   round(thetayc_summary_table["Adenoma", "PERMANOVA"], digits = 3)), 
             x = 0.20, y = -0.50, size = 2.5) + coord_fixed(ratio = 1), 
  
  ggplot(data = thetayc_crc_IF, 
         aes(x=NMDS1, y=NMDS2)) + 
    geom_point(aes(color=samples)) + theme_bw() + 
    coord_equal() + ggtitle("D") + 
    stat_ellipse(aes(group = samples, color = samples, fill = samples), 
                 alpha = 0.25, geom = "polygon") + 
    theme(plot.title = element_text(face = "bold", hjust = 0), 
          legend.position = "none") + 
    annotate("text", label = paste("PERMANOVA = ", 
                                   round(thetayc_summary_table["Carcinoma", "PERMANOVA"], digits = 3)), 
             x = 0.40, y = -0.50, size = 2.5) + coord_fixed(ratio = 1.20), ncol = 2, nrow = 2)


ggsave(file = "results/figures/Figure1.pdf", differences_graph, 
       width=8, height = 8, dpi = 300)
















