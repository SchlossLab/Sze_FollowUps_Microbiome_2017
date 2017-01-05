### Figure S2
### Description of variables in the both the lesion and IF models
## Marc Sze

#Load needed code and packages
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", "wesanderson", "knitr", "rmarkdown"))

#Read data needed
lesion_MDA_data_summary <- read.csv("results/tables/lesion_model_top_vars_MDA_Summary.csv", 
                            header = T, stringsAsFactors = F)

lesion_MDA_full <- read.csv("results/tables/lesion_model_top_vars_MDA_full_data.csv", 
                            header = T, stringsAsFactors = F)

tax <- read.delim('data/process/final.taxonomy', sep='\t', header=T, row.names=1)


# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)

ggplot(lesion_MDA_full, aes(factor(variables, 
                                   levels = rev(unique(lesion_MDA_full$variables)), 
                                   labels = rev(unique(lesion_MDA_full$variables))), 
                                   log10(value))) + 
  geom_point(aes(color = variable)) + stat_summary(fun.y = "mean", colour = "black", geom = "point", size = 3) + 
  coord_flip() + theme_bw() + ylab("Log10 MDA") + xlab("Variable") +  
  theme(legend.position = "none")




# If wanting to add a table and figure together
#grid.arrange(test1, test2, nrow = 1, ncol = 2)
#testTheme <- ttheme_minimal(core = list(fg_params = list(cex = 0.60)), 
#                            colhead = list(fg_params=list(cex = 0.60)))
#
#test2 <- tableGrob(imp_model_vars, rows = NULL, theme = testTheme)


