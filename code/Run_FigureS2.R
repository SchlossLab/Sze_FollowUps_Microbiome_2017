### Figure S2
### Summary of Important OTUs within the model and their MDAs
## Marc Sze

#Load needed code and packages
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", "wesanderson", "knitr", "rmarkdown"))

#Read data needed
adn_MDA_data_summary <- read.csv("data/process/tables/adn_MDA_Summary.csv", 
                                    header = T, stringsAsFactors = F) %>% slice(1:50)

adn_MDA_full <- read.csv("data/process/tables/adn_raw_mda_values.csv", 
                            header = T, stringsAsFactors = F) %>% filter(Variable %in% adn_MDA_data_summary$Variable)

srn_MDA_data_summary <- read.csv("data/process/tables/srn_MDA_Summary.csv", 
                                 header = T, stringsAsFactors = F) %>% slice(1:50)

srn_MDA_full <- read.csv("data/process/tables/srn_raw_mda_values.csv", 
                         header = T, stringsAsFactors = F) %>% filter(Variable %in% srn_MDA_data_summary$Variable)

crc_MDA_data_summary <- read.csv("data/process/tables/crc_MDA_Summary.csv", 
                                 header = T, stringsAsFactors = F) %>% slice(1:50)

crc_MDA_full <- read.csv("data/process/tables/crc_raw_mda_values.csv", 
                         header = T, stringsAsFactors = F) %>% filter(Variable %in% crc_MDA_data_summary$Variable)

# Generate the files needed for the MDA graph (ADN)
otu_num <- as.numeric(gsub("Otu", "", adn_MDA_data_summary$Variable))
low_tax_ID <- adn_MDA_data_summary$tax_ID

adn_graph_labels <- paste(gsub("2", "", adn_MDA_data_summary$tax_ID), " (OTU", otu_num, ")", sep = "")

test <- c()
for(i in 1:length(low_tax_ID)){
  
  test <- c(test, bquote(paste(italic(.(low_tax_ID[i])) ~ "(OTU", .(otu_num[i]), ")", sep = "")))
  
}

adn_labels <- do.call(expression, test)


# Generate the files needed for the MDA graph (SRN)
otu_num <- as.numeric(gsub("Otu", "", srn_MDA_data_summary$Variable))
low_tax_ID <- srn_MDA_data_summary$tax_ID

srn_graph_labels <- paste(gsub("2", "", srn_MDA_data_summary$tax_ID), " (OTU", otu_num, ")", sep = "")

test <- c()
for(i in 1:length(low_tax_ID)){
  
  test <- c(test, bquote(paste(italic(.(low_tax_ID[i])) ~ "(OTU", .(otu_num[i]), ")", sep = "")))
}

srn_labels <- do.call(expression, test)

# Generate the files needed for the MDA graph (CRC)
otu_num <- as.numeric(gsub("Otu", "", crc_MDA_data_summary$Variable))
low_tax_ID <- crc_MDA_data_summary$tax_ID

crc_graph_labels <- paste(gsub("2", "", crc_MDA_data_summary$tax_ID), " (OTU", otu_num, ")", sep = "")

test <- c()
for(i in 1:length(low_tax_ID)){
  
  test <- c(test, bquote(paste(italic(.(low_tax_ID[i])) ~ "(OTU", .(otu_num[i]), ")", sep = "")))
}

crc_labels <- do.call(expression, test)

# Create decimal control
fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}


# MDA Graph
adn_MDA_graph <- adn_MDA_full %>% filter(Overall > 0) %>% ggplot(aes(factor(Variable, 
                                   levels = rev(unique(adn_MDA_data_summary$Variable)), 
                                   labels = rev(adn_graph_labels)), 
                            log10(Overall))) + 
  geom_point(color = '#76EE00') + stat_summary(fun.y = "median", colour = '#006400', geom = "point", size = 2.5) + 
  coord_flip(ylim = c(-1.5, 1)) + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("A") + 
  scale_x_discrete(labels = rev(adn_labels)) + scale_y_continuous(labels = fmt_dcimals(0)) + 
  theme(plot.title = element_text(face = "bold", hjust = -0.56, size = 20), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 6))


srn_MDA_graph <- srn_MDA_full %>% filter(Overall > 0) %>% ggplot(aes(factor(Variable, 
                                                 levels = rev(unique(srn_MDA_data_summary$Variable)), 
                                                 labels = rev(srn_graph_labels)), 
                                          log10(Overall))) + 
  geom_point(color = '#F0E68C') + stat_summary(fun.y = "median", colour = '#EEC900', geom = "point", size = 2.5) + 
  coord_flip(ylim = c(-1.5, 1)) + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("B") + 
  scale_x_discrete(labels = rev(srn_labels)) + 
  theme(plot.title = element_text(face = "bold", hjust = -0.56, size = 20), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 6))

crc_MDA_graph <- crc_MDA_full %>% filter(Overall > 0) %>% ggplot(aes(factor(Variable, 
                                                 levels = rev(unique(crc_MDA_data_summary$Variable)), 
                                                 labels = rev(crc_graph_labels)), 
                                          log10(Overall))) + 
  geom_point(color = '#FFB6C1') + stat_summary(fun.y = "median", colour = '#DC143C', geom = "point", size = 2.5) + 
  coord_flip(ylim = c(-1.5, 1)) + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("C") + 
  scale_x_discrete(labels = rev(crc_labels)) + 
  theme(plot.title = element_text(face = "bold", hjust = -0.56, size = 20), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 6))


full_figS2 <- grid.arrange(adn_MDA_graph, srn_MDA_graph, crc_MDA_graph, nrow = 1)



ggsave(file = "results/figures/FigureS2.pdf", full_figS2, 
       width=12, height = 8, dpi = 300)



