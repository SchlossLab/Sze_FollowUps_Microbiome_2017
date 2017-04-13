# Rough comparison of common OTUs
### Most concerned with commonalities between models (treatment and normal)
### Marc Sze


source('code/functions.R')

loadLibs(c("tidyr", "dplyr", "scales", "VennDiagram", "ggplot2"))


# Load necessary data
adn_data <- read.csv("data/process/tables/adn_reduced_model_top_vars_MDA_Summary.csv", 
                     header = T, stringsAsFactors = F)
srn_data <- read.csv("data/process/tables/srn_reduced_model_top_vars_MDA_Summary.csv", 
                      header = T, stringsAsFactors = F)
crc_data <- read.csv("data/process/tables/reduced_crc_model_top_vars_MDA_Summary.csv", 
                      header = T, stringsAsFactors = F)

#Generate vectors of OTUs for each diagnosis
adn_OTUs <- adn_data$otu
srn_OTUs <- srn_data$otu
crc_data <- crc_data$otu


test <- calculate.overlap(list(adn_OTUs, srn_OTUs, crc_data))

test2 <- as.data.frame.list(lapply(test, function(x) length(x)))


venn_graph <- as.data.frame(cbind(g1 = c(4.5, 5.4, 6.2), g2 = c(7, 5.5, 7))) %>% 
  mutate(group = c("adn", "srn", "crc")) %>% 
  ggplot(aes(x = g1, y = g2)) + 
  geom_point(aes(color = group), size = 105, alpha = 0.5) + 
  scale_color_manual(values = c('#228B22', '#FFD700', '#DC143C')) + 
  coord_cartesian(xlim = c(1:10), ylim = c(1:10)) + xlab("") + ylab("") + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.ticks = element_blank()) + 
  annotate("text", 
           label = 'atop(bold("Adenoma"))', x = 2.5, y = 8.7, size = 4, parse = TRUE) + 
  annotate("text", 
           label = 'atop(bold("Advanced Adenoma"))', x = 8, y = 8.7, size = 4, parse = TRUE) + 
  annotate("text", 
           label = 'atop(bold("Carcinoma"))', x = 5.5, y = 3.4, size = 4, parse = TRUE) + 
  annotate("text", 
           label = test2$a1, x = 3, y = 7.2, size = 14, parse = TRUE) + 
  annotate("text", 
           label = test2$a2, x = 5.4, y = 7.7, size = 14, parse = TRUE) + 
  annotate("text", 
           label = test2$a3, x = 7.6, y = 7.2, size = 14, parse = TRUE) + 
  annotate("text", 
           label = test2$a4, x = 4, y = 5.8, size = 14, parse = TRUE) + 
  annotate("text", 
           label = test2$a5, x = 5.4, y = 6.3, size = 14, parse = TRUE) + 
  annotate("text", 
           label = test2$a6, x = 7, y = 5.8, size = 14, parse = TRUE) + 
  annotate("text", 
           label = test2$a7, x = 5.4, y = 4.8, size = 14, parse = TRUE)


ggsave(file = "results/figures/venn.tiff", venn_graph, 
       width=7, height = 9, dpi = 300)


## Extra venn stuff
## Uses default venn diagram to plot for base
#test <- venn.diagram(x = list(adn_OTUs, srn_OTUs, crc_data), 
#                     category.names = c("Adenoma", "Advanced\nAdenoma", "Carcinoma"), 
#                     height = 6000, width = 9000, 
#                     resolution = 500, 
#                     filename = "results/figures/test.tiff", 
#                     compression = "lzw",
#                     imagetype = "tiff", 
#                     fill = c('#006400', '#EEC900', '#DC143C'), 
#                     lwd = 2, 
#                     lty = 'blank', 
#                     cex = 3, 
#                     fontface = "bold", 
#                     cat.cex = 4, 
#                     cat.fontface = "bold")


### treatment model random info on model similarities to normal vs diagnosis models
#adn_treat <- read.csv("data/process/tables/adn_treatment_imp_vars_summary.csv")
#srn_treat <- read.csv("data/process/tables/srn_treatment_imp_vars_summary.csv")
#crc_treat <- read.csv("data/process/tables/crc_treatment_imp_vars_summary.csv")
#common_norm <- read.csv('data/process/tables/pvalue_adn_srn_crc_common_imp_vars.csv', 
#                                                       header = T, row.names = 1, stringsAsFactors = F)


#adn_otus <- as.character(adn_treat$Variable[1:round(length(adn_treat$Variable)*0.10)])
#srn_otus <- as.character(srn_treat$Variable[1:round(length(srn_treat$Variable)*0.10)])
#crc_otus <- as.character(crc_treat$Variable[1:round(length(crc_treat$Variable)*0.10)])

#length(rownames(filter(common_norm, otu %in% adn_otus) %>% select(otu)))
#length(rownames(filter(common_norm, otu %in% srn_otus) %>% select(otu)))
#length(rownames(filter(common_norm, otu %in% crc_otus) %>% select(otu)))


