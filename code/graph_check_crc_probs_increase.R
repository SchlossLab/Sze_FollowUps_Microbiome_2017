### Graph increased versus decreased probs
### create graphs to visualise differences (if there are any)
## Marc Sze

## Load in necessary libraries
source('code/functions.R')
loadLibs(c("dplyr", "tidyr", "ggplot2", "gridExtra"))


## Load needed data files
residents <- read.csv("data/process/tables/inc_probs_crc_oral_residents_data.csv", 
                      header = T, stringsAsFactors = F)


## Graph the data

oral_p <- ggplot(residents, aes(factor(sampleType, 
                             levels = c("initial", "followup"), 
                             labels = c("Pre", "Post")), oral_path)) + 
  geom_jitter(aes(color = probs_increase), width = 0.1) + 
  scale_color_manual(name = "Carcinoma\nProbabiliy\nIncrease", 
                     values = c('#8968CD', '#FF00FF')) + 
  ggtitle("A") + theme_bw() +  xlab("") + ylab("% Relative Abundance") + 
  stat_summary(aes(group = probs_increase), fun.y = median, 
               colour = c('#8968CD', '#8968CD', '#FF00FF', '#FF00FF'), geom = "point", size = 4.5) + 
  coord_cartesian(ylim = c(0, 15)) + 
  theme(plot.title = element_text(face = "bold", hjust = -0.08, size = 20), 
        legend.position = c(0.9, 0.4)) + 
  annotate("text", 
          label = paste("Oral Pathogens"),x = 1.5, y = 15, size = 4)



non_oral_p <- ggplot(residents, aes(factor(sampleType, 
                             levels = c("initial", "followup"), 
                             labels = c("Pre", "Post")), residents)) + 
  geom_jitter(aes(color = probs_increase), width = 0.1) + 
  scale_color_manual(name = "Carcinoma\nProbabiliy\nIncrease", 
                     values = c('#8968CD', '#FF00FF')) + 
  ggtitle("B") + theme_bw() +  xlab("") + ylab("% Relative Abundance") + 
  stat_summary(aes(group = probs_increase), fun.y = median, 
               colour = c('#8968CD', '#8968CD', '#FF00FF', '#FF00FF'), geom = "point", size = 4.5) + 
  coord_cartesian(ylim = c(0, 100)) + 
  theme(plot.title = element_text(face = "bold", hjust = -0.09, size = 20), 
        legend.position = "none") + 
  annotate("text", 
           label = paste("Non-Oral Pathogens"),x = 1.5, y = 100, size = 4)


crc_only_group <- ggplot(residents, aes(factor(sampleType, 
                                           levels = c("initial", "followup"), 
                                           labels = c("Pre", "Post")), crc_only)) + 
  geom_jitter(aes(color = probs_increase), width = 0.1) + 
  scale_color_manual(name = "Carcinoma\nProbabiliy\nIncrease", 
                     values = c('#8968CD', '#FF00FF')) + 
  ggtitle("B") + theme_bw() +  xlab("") + ylab("% Relative Abundance") + 
  stat_summary(aes(group = probs_increase), fun.y = median, 
               colour = c('#8968CD', '#8968CD', '#FF00FF', '#FF00FF'), geom = "point", size = 4.5) + 
  coord_cartesian(ylim = c(0, 75)) + 
  theme(plot.title = element_text(face = "bold", hjust = -0.09, size = 20), 
        legend.position = "none") + 
  annotate("text", 
           label = paste("Carcinoma Higher"),x = 1.5, y = 75, size = 4)


normal_only_group <- ggplot(residents, aes(factor(sampleType, 
                                               levels = c("initial", "followup"), 
                                               labels = c("Pre", "Post")), normal_only)) + 
  geom_jitter(aes(color = probs_increase), width = 0.1) + 
  scale_color_manual(name = "Carcinoma\nProbabiliy\nIncrease", 
                     values = c('#8968CD', '#FF00FF')) + 
  ggtitle("C") + theme_bw() +  xlab("") + ylab("% Relative Abundance") + 
  stat_summary(aes(group = probs_increase), fun.y = median, 
               colour = c('#8968CD', '#8968CD', '#FF00FF', '#FF00FF'), geom = "point", size = 4.5) + 
  coord_cartesian(ylim = c(0, 75)) + 
  theme(plot.title = element_text(face = "bold", hjust = -0.09, size = 20), 
        legend.position = "none") + 
  annotate("text", 
           label = paste("Normal Higher"),x = 1.5, y = 75, size = 4)
  
  
combined_graph <- grid.arrange(oral_p, crc_only_group, normal_only_group, ncol = 1)

ggsave(file = "results/figures/FigureS4.pdf", combined_graph, 
       width = 6, height = 8, dpi = 300)


