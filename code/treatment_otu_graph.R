### Plot for FHCRC
### Common OTUs
### Marc Sze

source('code/functions.R')

loadLibs(c("tidyr", "dplyr", "scales", "ggplot2", "gridExtra"))

# Load needed data
adn_treat_varsMDA <- read.csv('data/process/tables/reduced_adn_treatment_top_vars_MDA_Summary.csv', 
                        header = T, stringsAsFactors = F) %>% 
  mutate(rank = seq(1:10), tax_id = gsub("_unclassified", "", tax_id))

srn_treat_varsMDA <- read.csv('data/process/tables/reduced_srn_treatment_top_vars_MDA_Summary.csv', 
                        header = T, stringsAsFactors = F) %>% 
  mutate(rank = seq(1:10), tax_id = gsub("_unclassified", "", tax_id))

crc_treat_varsMDA <- read.csv('data/process/tables/reduced_crc_treatment_top_vars_MDA_Summary.csv', 
                        header = T, stringsAsFactors = F) %>% 
  mutate(rank = seq(1:10), tax_id = gsub("_unclassified", "", tax_id))


adn_data <- read.csv('data/process/tables/reduced_adn_treatment_top_vars_MDA_full_data.csv', 
                     header = T, stringsAsFactors = F)

srn_data <- read.csv('data/process/tables/reduced_srn_treatment_top_vars_MDA_full_data.csv', 
                     header = T, stringsAsFactors = F)
crc_data <- read.csv('data/process/tables/reduced_crc_treatment_top_vars_MDA_full_data.csv', 
                     header = T, stringsAsFactors = F)

# Create adn labels
otu_num <- as.numeric(gsub("Otu", "", adn_treat_varsMDA$otu))

test <- c()
for(i in 1:length(otu_num)){
  
  test <- c(test, 
            bquote(paste(italic(.(gsub("_", " ", adn_treat_varsMDA$tax_id[i]))) ~ "(OTU", .(otu_num[i]), ")", sep = "")))
}

adn_labels <- do.call(expression, test)
adn_axis <- paste(adn_treat_varsMDA$tax_id, " (OTU", otu_num, ")", sep="")


# Create srn labels
otu_num <- as.numeric(gsub("Otu", "", srn_treat_varsMDA$otu))

test <- c()
for(i in 1:length(otu_num)){
  
  test <- c(test, 
            bquote(paste(italic(.(gsub("_", " ", srn_treat_varsMDA$tax_id[i]))) ~ "(OTU", .(otu_num[i]), ")", sep = "")))
}

srn_labels <- do.call(expression, test)
srn_axis <- paste(adn_treat_varsMDA$tax_id, " (OTU", otu_num, ")", sep="")


# Create crc labels
otu_num <- as.numeric(gsub("Otu", "", crc_treat_varsMDA$otu))

test <- c()
for(i in 1:length(otu_num)){
  
  test <- c(test, 
            bquote(paste(italic(.(gsub("_", " ", crc_treat_varsMDA$tax_id[i]))) ~ "(OTU", .(otu_num[i]), ")", sep = "")))
}

crc_labels <- do.call(expression, test)
crc_axis <- paste(adn_treat_varsMDA$tax_id, " (OTU", otu_num, ")", sep="")


# Create decimal control
fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}


# Create Graph
adn_treat_MDA_graph <- ggplot(adn_data, aes(factor(otu, 
                                levels = rev(unique(adn_data$otu)), 
                                labels = rev(adn_axis)), log10(value))) + 
  geom_point(color = '#76EE00') + stat_summary(fun.y = "median", colour = '#006400', geom = "point", size = 2.5) + 
  coord_flip() + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("A") + 
  scale_x_discrete(labels = rev(adn_labels)) + scale_y_continuous(labels = fmt_dcimals(0)) + 
  theme(plot.title = element_text(face = "bold"), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 6))


srn_treat_MDA_graph <- ggplot(srn_data, aes(factor(otu, 
                                                 levels = rev(unique(srn_data$otu)), 
                                                 labels = rev(srn_axis)), log10(value))) + 
  geom_point(color = '#F0E68C') + stat_summary(fun.y = "median", colour = '#EEC900', geom = "point", size = 2.5) + 
  coord_flip() + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("B") + 
  scale_x_discrete(labels = rev(srn_labels)) + 
  theme(plot.title = element_text(face = "bold"), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 6))


crc_treat_MDA_graph <- ggplot(crc_data, aes(factor(otu, 
                                                 levels = rev(unique(crc_data$otu)), 
                                                 labels = rev(crc_axis)), log10(value))) + 
  geom_point(color = '#FFB6C1') + stat_summary(fun.y = "median", colour = '#DC143C', geom = "point", size = 2.5) + 
  coord_flip() + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("C") + 
  scale_x_discrete(labels = rev(crc_labels)) + 
  theme(plot.title = element_text(face = "bold"), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 6))


full_fig2 <- grid.arrange(adn_treat_MDA_graph, srn_treat_MDA_graph, crc_treat_MDA_graph, nrow = 1)

ggsave(file = "results/figures/Figure2.pdf", full_fig2, 
       width=11, height = 6, dpi = 300)
