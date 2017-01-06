### Figure S3
### Description of variables in IF models
## Marc Sze


#Load needed code and packages
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", "wesanderson", "knitr", "rmarkdown"))

#Read data needed
IF_MDA_data_summary <- read.csv("results/tables/IF_model_top_vars_MDA_Summary.csv", 
                                    header = T, stringsAsFactors = F)

IF_MDA_full <- read.csv("results/tables/IF_model_top_vars_MDA_full_data.csv", 
                            header = T, stringsAsFactors = F)

tax <- read.delim('data/process/final.taxonomy', sep='\t', header=T, row.names=1)

occurance_data <- read.csv("results/tables/IF_rf_wCV_imp_vars_summary.csv", header = T, stringsAsFactors = F)

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)

# Pull OTUs that are only in the MDA data
OTU_IDs <- unique(filter(IF_MDA_full, variables != "fit_result")[, "variables"])
select_tax_df <- tax_df[OTU_IDs, ]
low_tax_ID <- gsub("_unclassified", "", createTaxaLabeller(select_tax_df))

# create labels for factor values with low taxonomy
graph_labels <- c("FIT", paste(names(low_tax_ID), " (", low_tax_ID, ")", sep = ""))

# align occurance data with MDA data
rownames(occurance_data) <- occurance_data$Variable
occurance_data <- occurance_data[unique(IF_MDA_full$variables), ]


# create part A plot
IFA <- ggplot(IF_MDA_full, aes(factor(variables, 
                                              levels = rev(unique(IF_MDA_full$variables)), 
                                              labels = rev(graph_labels)), 
                                       log10(value))) + 
  geom_jitter(aes(color = variable)) + 
  coord_flip() + theme_bw() + ylab("Log10 MDA") + xlab("Variable") +  ggtitle("A") + 
  theme(plot.title = element_text(face = "bold", size = 32, hjust = 0), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 8))

# create part B plot

IFB <- ggplot(occurance_data, aes(factor(Variable, 
                                             levels = rev(unique(IF_MDA_full$variables)), 
                                             labels = rev(graph_labels)), total_appearance)) + 
  geom_bar(stat = "identity", fill = "white", color = "black", width = 0.75) + 
  geom_hline(aes(yintercept = 50, color = "red"), linetype = 2, size = 1) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110)) + 
  coord_flip() + theme_bw() + 
  ylab("Total Model Appearance Percentage") + xlab("Variable") + ggtitle("B") + 
  theme(plot.title = element_text(face = "bold", size = 32, hjust = 0), 
        axis.title = element_text(face = "bold"), 
        legend.position = "none", 
        axis.text.y = element_text(size = 8), 
        panel.grid = element_blank())


test <- grid.arrange(IFA, IFB, nrow = 2, ncol = 1)

ggsave(file = "results/figures/FigureS3.pdf", test, 
       width=8, height = 10, dpi = 300)







