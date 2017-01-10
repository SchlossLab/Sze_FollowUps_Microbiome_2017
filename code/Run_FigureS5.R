### Figure S4
### Distribution of P-values from Paired Wilcoxson Analysis
## Marc Sze

#Load needed code and packages
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", "wesanderson", "knitr", "rmarkdown"))

#Read data needed
diff_table_time <- read.csv("results/tables/time_datatable.csv", header = T)

# Plot the change in distance versus time
time_graph <- ggplot(diff_table_time, aes(x = time, y = distance)) + 
  geom_point(aes(color = dx), size = 4) + theme_bw() + 
  scale_color_manual(name = "Lesion Type", values = wes_palette("GrandBudapest"), 
                     breaks = c("adenoma", "cancer"), labels = c("Adenoma", "Cancer")) + 
  coord_cartesian(ylim = c(0, 1)) + ylab("Thetayc Distance") + xlab("Time from Initial (Days)") + 
  theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"))

ggsave(file = "results/figures/FigureS5.pdf", time_graph, 
       width=6, height = 8, dpi = 300)