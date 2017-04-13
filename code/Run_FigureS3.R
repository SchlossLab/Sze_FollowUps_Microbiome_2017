### Figure S3
### Summary of Important OTUs within the model and their MDAs
## Marc Sze

#Load needed code and packages
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", "wesanderson", "knitr", "rmarkdown"))

#Read data needed
adn_MDA_data_summary <- read.csv("data/process/tables/adn_reduced_model_top_vars_MDA_Summary.csv", 
                                    header = T, stringsAsFactors = F)

adn_MDA_full <- read.csv("data/process/tables/adn_reduced_lesion_model_top_vars_MDA.csv", 
                            header = T, stringsAsFactors = F)

srn_MDA_data_summary <- read.csv("data/process/tables/srn_reduced_model_top_vars_MDA_Summary.csv", 
                                 header = T, stringsAsFactors = F)

srn_MDA_full <- read.csv("data/process/tables/srn_reduced_model_top_vars_MDA.csv", 
                         header = T, stringsAsFactors = F)

crc_MDA_data_summary <- read.csv("data/process/tables/reduced_crc_model_top_vars_MDA_Summary.csv", 
                                 header = T, stringsAsFactors = F)

crc_MDA_full <- read.csv("data/process/tables/crc_reduced_lesion_model_top_vars_MDA.csv", 
                         header = T, stringsAsFactors = F)


tax <- read.delim('data/process/final.taxonomy', sep='\t', header=T, row.names=1)

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)


# Generate the files needed for the MDA graph (ADN)

# Order the data file from most to least important based on mean MDA
adn_MDA_data_summary <- adn_MDA_data_summary[order(adn_MDA_data_summary$median_MDA, decreasing = TRUE), ]

# Pull OTUs that are only in the MDA data
OTU_IDs <- unique(filter(adn_MDA_data_summary, otu != "fit_result")[, "otu"])
select_tax_df <- tax_df[OTU_IDs, ]
low_tax_ID <- gsub("_", " ", gsub("2", "", gsub("_unclassified", "", createTaxaLabeller(select_tax_df))))
otu_num <- as.numeric(gsub("Otu", "", OTU_IDs))

# create labels for factor values with low taxonomy
adn_graph_labels <- paste(gsub("2", "", low_tax_ID), " (OTU", otu_num, ")", sep = "")

test <- c()
for(i in 1:length(low_tax_ID)){
  
  test <- c(test, bquote(paste(italic(.(low_tax_ID[i])) ~ "(OTU", .(otu_num[i]), ")", sep = "")))
  
}

adn_labels <- do.call(expression, test)


# Generate the files needed for the MDA graph (SRN)

# Order the data file from most to least important based on mean MDA
srn_MDA_data_summary <- srn_MDA_data_summary[order(srn_MDA_data_summary$median_MDA, decreasing = TRUE), ]

# Pull OTUs that are only in the MDA data
OTU_IDs <- unique(filter(srn_MDA_data_summary, otu != "fit_result")[, "otu"])
select_tax_df <- tax_df[OTU_IDs, ]
low_tax_ID <- gsub("_", " ", gsub("2", "", gsub("_unclassified", "", createTaxaLabeller(select_tax_df))))
otu_num <- as.numeric(gsub("Otu", "", OTU_IDs))

# create labels for factor values with low taxonomy
srn_graph_labels <- paste(gsub("2", "", low_tax_ID), " (OTU", otu_num, ")", sep = "")

test <- c()
for(i in 1:length(low_tax_ID)){
  
  test <- c(test, bquote(paste(italic(.(low_tax_ID[i])) ~ "(OTU", .(otu_num[i]), ")", sep = "")))
}

srn_labels <- do.call(expression, test)

# Generate the files needed for the MDA graph (CRC)

# Order the data file from most to least important based on mean MDA
crc_MDA_data_summary <- crc_MDA_data_summary[order(crc_MDA_data_summary$median_MDA, decreasing = TRUE), ]

# Pull OTUs that are only in the MDA data
OTU_IDs <- unique(filter(crc_MDA_data_summary, otu != "fit_result")[, "otu"])
select_tax_df <- tax_df[OTU_IDs, ]
low_tax_ID <- gsub("_", " ", gsub("2", "", gsub("_unclassified", "", createTaxaLabeller(select_tax_df))))
otu_num <- as.numeric(gsub("Otu", "", OTU_IDs))

# create labels for factor values with low taxonomy
crc_graph_labels <- paste(gsub("2", "", low_tax_ID), " (OTU", otu_num, ")", sep = "")

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
adn_MDA_graph <- ggplot(adn_MDA_full, aes(factor(otu, 
                                   levels = rev(unique(adn_MDA_full$otu)), 
                                   labels = rev(adn_graph_labels)), 
                            log10(value))) + 
  geom_point(color = '#76EE00') + stat_summary(fun.y = "median", colour = '#006400', geom = "point", size = 2.5) + 
  coord_flip() + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("A") + 
  scale_x_discrete(labels = rev(adn_labels)) + scale_y_continuous(labels = fmt_dcimals(0)) + 
  theme(plot.title = element_text(face = "bold"), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 6))


srn_MDA_graph <- ggplot(srn_MDA_full, aes(factor(otu, 
                                                 levels = rev(unique(srn_MDA_full$otu)), 
                                                 labels = rev(srn_graph_labels)), 
                                          log10(value))) + 
  geom_point(color = '#F0E68C') + stat_summary(fun.y = "median", colour = '#EEC900', geom = "point", size = 2.5) + 
  coord_flip() + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("B") + 
  scale_x_discrete(labels = rev(srn_labels)) + 
  theme(plot.title = element_text(face = "bold"), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 6))

crc_MDA_graph <- ggplot(crc_MDA_full, aes(factor(otu, 
                                                 levels = rev(unique(crc_MDA_full$otu)), 
                                                 labels = rev(crc_graph_labels)), 
                                          log10(value))) + 
  geom_point(color = '#FFB6C1') + stat_summary(fun.y = "median", colour = '#DC143C', geom = "point", size = 2.5) + 
  coord_flip() + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("C") + 
  scale_x_discrete(labels = rev(crc_labels)) + 
  theme(plot.title = element_text(face = "bold"), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 6))


full_figS3 <- grid.arrange(adn_MDA_graph, srn_MDA_graph, crc_MDA_graph, nrow = 1)



ggsave(file = "results/figures/FigureS3.pdf", full_figS3, 
       width=12, height = 8, dpi = 300)



