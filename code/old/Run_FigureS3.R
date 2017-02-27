### Figure S3
### Description of variables in IF models
## Marc Sze


#Load needed code and packages
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", "wesanderson"))

#Read data needed
IF_MDA_data_summary <- read.csv("results/tables/reduced_IF_model_top_vars_MDA_Summary.csv", 
                                    header = T, stringsAsFactors = F)

IF_MDA_full <- read.csv("results/tables/reduced_IF_model_top_vars_MDA.csv", 
                            header = T, stringsAsFactors = F)

tax <- read.delim('data/process/final.taxonomy', sep='\t', header=T, row.names=1)

occurance_data <- read.csv("results/tables/IF_rf_wCV_imp_vars_summary.csv", header = T, stringsAsFactors = F)

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)

# Order the data file from most to least important based on mean MDA
IF_MDA_data_summary <- IF_MDA_data_summary[order(IF_MDA_data_summary$mean_MDA, decreasing = TRUE), ]

# Pull OTUs that are only in the MDA data
OTU_IDs <- unique(IF_MDA_data_summary$variable)
select_tax_df <- tax_df[OTU_IDs, ]
low_tax_ID <- gsub("2", "", gsub("_unclassified", "", createTaxaLabeller(select_tax_df)))

# create labels for factor values with low taxonomy
graph_labels <- c(paste(gsub("2", "", low_tax_ID), " (", names(low_tax_ID), ")", sep = ""))
OTU_names <- names(low_tax_ID)

test <- c()
for(i in 1:length(low_tax_ID)){
  if(is.na(low_tax_ID[i])){
    
    test <- c(test, bquote(paste("FIT", sep = "")))
  } else {
    
    test <- c(test, bquote(paste(italic(.(low_tax_ID[i])) ~ "(", .(OTU_names[i]), ")", sep = "")))
  }
  
}

test2 <- do.call(expression, test)

# align occurance data with MDA data
rownames(occurance_data) <- occurance_data$Variable
occurance_data <- occurance_data[unique(IF_MDA_full$variables), ]


# create part A plot
IFA <- ggplot(IF_MDA_full, aes(factor(variables, 
                                              levels = rev(unique(IF_MDA_data_summary$variable)), 
                                              labels = rev(graph_labels)), 
                                       log10(value))) + 
  geom_point(aes(color = variable)) + stat_summary(fun.y = "median", colour = "black", geom = "point", size = 2) + 
  coord_flip() + theme_bw() + ylab("Log10 MDA") + xlab("Variable") +  ggtitle("A") + 
  scale_x_discrete(labels = rev(test2)) + 
  theme(plot.title = element_text(face = "bold", size = 32, hjust = 0), 
        legend.position = "none", 
        axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 6))

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
        axis.text.y = element_text(size = 4.5), 
        panel.grid = element_blank())


test <- grid.arrange(IFA, IFB, nrow = 2, ncol = 1)

ggsave(file = "results/figures/FigureS3.pdf", test, 
       width=8, height = 10, dpi = 300)







