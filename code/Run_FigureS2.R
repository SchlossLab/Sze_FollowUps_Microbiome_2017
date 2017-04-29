### Figure S2
### ROC Curves summary for three models used
## Marc Sze

#Load needed code and packages
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", "wesanderson", "knitr", "rmarkdown"))

#Read data needed
adn_model_roc <- read.csv("data/process/tables/adn_reduced_test_data_roc.csv", header = T)
srn_model_roc <- read.csv("data/process/tables/srn_reduced_test_data_roc.csv", header = T)
crc_model_roc <- read.csv("data/process/tables/crc_reduced_test_data_roc.csv", header = T)


# Create the graph

rocs_graph <- grid.arrange(
  # Adenoma ROC curve information
  filter(adn_model_roc, run != "middle_roc" & run != "full_roc") %>% 
    ggplot(aes(x = sensitivities, y = specificities)) + 
    geom_polygon(data = filter(adn_model_roc, run != "full_roc"), alpha = 0.5, fill = '#76EE00') + 
    geom_line(data = filter(adn_model_roc, run != "full_roc"), 
              aes(group = run), size = 1.25, color = '#76EE00') + 
    geom_line(data = filter(adn_model_roc, run == "full_roc"), 
              size = 1.5, color = '#006400') + 
    scale_x_continuous(trans = "reverse") + theme_bw() + 
    ggtitle("A") + xlab("Sensitivity") + ylab("Specificity") + 
    theme(plot.title = element_text(face= "bold", hjust = -0.19, size = 20), 
          axis.title = element_text(face = "bold")), 
  
  # SRN ROC curve information
  filter(srn_model_roc, run != "middle_roc" & run != "full_roc") %>% 
    ggplot(aes(x = sensitivities, y = specificities)) + 
    geom_polygon(data = filter(srn_model_roc, run != "full_roc"), alpha = 0.5, fill = '#F0E68C') + 
    geom_line(data = filter(srn_model_roc, run != "full_roc"), 
              aes(group = run), size = 1.25, color = '#F0E68C') + 
    geom_line(data = filter(srn_model_roc, run == "full_roc"), 
              size = 1.5, color = '#EEC900') + 
    scale_x_continuous(trans = "reverse") + theme_bw() +  
    ggtitle("B") + xlab("Sensitivity") + ylab("Specificity") + 
    theme(plot.title = element_text(face= "bold", hjust = -0.19, size = 20), 
          axis.title = element_text(face = "bold")), 
  
  # CRC ROC curve information
  filter(crc_model_roc, run != "middle_roc" & run != "full_roc") %>% 
    ggplot(aes(x = sensitivities, y = specificities)) + 
    geom_polygon(data = filter(crc_model_roc, run != "full_roc"), alpha = 0.5, fill = '#FFB6C1') + 
    geom_line(data = filter(crc_model_roc, run != "full_roc"), 
              aes(group = run), size = 1.25, color = '#FFB6C1') + 
    geom_line(data = filter(crc_model_roc, run == "full_roc"), 
              size = 1.5, color = '#DC143C') + 
    scale_x_continuous(trans = "reverse") + theme_bw() +  
    ggtitle("C") + xlab("Sensitivity") + ylab("Specificity") + 
    theme(plot.title = element_text(face= "bold", hjust = -0.19, size = 20), 
          axis.title = element_text(face = "bold")))


ggsave(file = "results/figures/FigureS2.tiff", rocs_graph, 
       width=4, height = 8, dpi = 300)
