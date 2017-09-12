### Plot for FHCRC
### Common OTUs
### Marc Sze

source('code/functions.R')

loadLibs(c("tidyr", "dplyr", "scales", "ggplot2", "gridExtra"))

# Load needed data and take only top 10
adn_treat_varsMDA <- read.csv('data/process/tables/adn_treatment_MDA_Summary.csv', 
                        header = T, stringsAsFactors = F) %>% slice(1:10)

srn_treat_varsMDA <- read.csv('data/process/tables/srn_treatment_MDA_Summary.csv', 
                        header = T, stringsAsFactors = F) %>% slice(1:10)

crc_treat_varsMDA <- read.csv('data/process/tables/crc_treatment_MDA_Summary.csv', 
                        header = T, stringsAsFactors = F) %>% slice(1:10)


adn_data <- read.csv('data/process/tables/adn_treatment_raw_mda_values.csv', 
                     header = T, stringsAsFactors = F) %>% 
  filter(Variable %in% adn_treat_varsMDA$Variable)

srn_data <- read.csv('data/process/tables/srn_treatment_raw_mda_values.csv', 
                     header = T, stringsAsFactors = F) %>% 
  filter(Variable %in% srn_treat_varsMDA$Variable)

crc_data <- read.csv('data/process/tables/crc_treatment_raw_mda_values.csv', 
                     header = T, stringsAsFactors = F) %>% 
  filter(Variable %in% crc_treat_varsMDA$Variable)

# Create adn labels
otu_num <- as.numeric(gsub("Otu", "", adn_treat_varsMDA$Variable))

test <- c()
for(i in 1:length(otu_num)){
  
  test <- c(test, 
            bquote(paste(italic(.(gsub("_", " ", 
                                       adn_treat_varsMDA$tax_ID[i]))) ~ "(OTU", .(otu_num[i]), ")", 
                         sep = "")))
}

adn_labels <- do.call(expression, test)
adn_axis <- paste(adn_treat_varsMDA$tax_ID, " (OTU", otu_num, ")", sep="")


# Create srn labels
otu_num <- as.numeric(gsub("Otu", "", srn_treat_varsMDA$Variable))

test <- c()
for(i in 1:length(otu_num)){
  
  test <- c(test, 
            bquote(paste(italic(.(gsub("_", " ", 
                                       srn_treat_varsMDA$tax_ID[i]))) ~ "(OTU", .(otu_num[i]), ")", 
                         sep = "")))
}

srn_labels <- do.call(expression, test)
srn_axis <- paste(adn_treat_varsMDA$tax_ID, " (OTU", otu_num, ")", sep="")


# Create crc labels
otu_num <- as.numeric(gsub("Otu", "", crc_treat_varsMDA$Variable))

test <- c()
for(i in 1:length(otu_num)){
  
  test <- c(test, 
            bquote(paste(italic(.(gsub("_", " ", 
                                       crc_treat_varsMDA$tax_ID[i]))) ~ "(OTU", .(otu_num[i]), ")", 
                         sep = "")))
}

crc_labels <- do.call(expression, test)
crc_axis <- paste(adn_treat_varsMDA$tax_ID, " (OTU", otu_num, ")", sep="")


# Create decimal control
fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}


# Create Graph
adn_treat_MDA_graph <- adn_data %>% filter(Overall > 0) %>% 
  ggplot(aes(factor(Variable, 
                    levels = rev(unique(adn_data$Variable)), 
                    labels = rev(adn_axis)), log10(Overall))) + 
  geom_point(color = '#76EE00') + stat_summary(fun.y = "median", colour = '#006400', geom = "point", size = 2.5) + 
  coord_flip(ylim = c(-3, 0.8)) + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("A") + 
  scale_x_discrete(labels = rev(adn_labels)) + scale_y_continuous(labels = fmt_dcimals(0)) + 
  theme(plot.title = element_text(face = "bold", hjust = -0.60, size = 20), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 10))


srn_treat_MDA_graph <- srn_data %>% filter(Overall > 0) %>% 
  ggplot(aes(factor(Variable, 
                    levels = rev(unique(srn_data$Variable)), 
                    labels = rev(srn_axis)), log10(Overall))) + 
  geom_point(color = '#F0E68C') + stat_summary(fun.y = "median", colour = '#EEC900', geom = "point", size = 2.5) + 
  coord_flip(ylim = c(-0.3, 0.8)) + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("B") + 
  scale_x_discrete(labels = rev(srn_labels)) + 
  theme(plot.title = element_text(face = "bold", hjust = -0.60, size = 20), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 10))


crc_treat_MDA_graph <- ggplot(crc_data, aes(factor(otu, 
                                                 levels = rev(unique(crc_data$otu)), 
                                                 labels = rev(crc_axis)), log10(value))) + 
  geom_point(color = '#FFB6C1') + stat_summary(fun.y = "median", colour = '#DC143C', geom = "point", size = 2.5) + 
  coord_flip(ylim = c(-0.3, 0.8)) + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("C") + 
  scale_x_discrete(labels = rev(crc_labels)) + 
  theme(plot.title = element_text(face = "bold", hjust = -0.60, size = 20), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 10))


full_fig2 <- grid.arrange(adn_treat_MDA_graph, srn_treat_MDA_graph, crc_treat_MDA_graph, nrow = 1)

ggsave(file = "results/figures/Figure2.pdf", full_fig2, 
       width=17, height = 6, dpi = 300)
