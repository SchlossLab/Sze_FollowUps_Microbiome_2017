### Create Figure 3 graph
### Show results for initial and follow up samples based on testing
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))


### Load in needed data tables
reduced_IF_model_roc <- read.csv("data/process/tables/reduced_IF_test_data_roc.csv", header = T)
IF_MDA_data_summary <- read.csv("data/process/tables/reduced_IF_model_top_vars_MDA_Summary.csv", 
                                header = T, stringsAsFactors = F)

IF_MDA_full <- read.csv("data/process/tables/reduced_IF_model_top_vars_MDA.csv", 
                        header = T, stringsAsFactors = F)

IF_red_follow_up_probability <- read.csv("data/process/tables/reduced_IF_follow_up_probability_summary.csv", 
                                         stringsAsFactors = F, header = T)

tax <- read.delim('data/process/final.taxonomy', sep='\t', header=T, row.names=1)

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
low_tax_ID <- gsub("_", " ", gsub("2", "", gsub("_unclassified", "", createTaxaLabeller(select_tax_df))))

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

# Create Figure

initial_model <- grid.arrange(
  
  # ROC Graph
  filter(reduced_IF_model_roc, run != "full_roc") %>% ggplot(aes(x = sensitivities, y = specificities)) + 
    geom_polygon(data = filter(reduced_IF_model_roc, run != "full_roc"), alpha = 0.5, fill = "gray") + 
    geom_line(data = filter(reduced_IF_model_roc, run != "full_roc"), 
              aes(group = run), size = 1.25, color = "gray") + 
    geom_line(data = filter(reduced_IF_model_roc, run == "full_roc"), 
              size = 1.5, color = "cornflowerblue") + 
    scale_x_continuous(trans = "reverse") + theme_bw() + geom_abline(intercept = 1, linetype = 2, size = 1) + 
    ggtitle("A") + xlab("Sensitivity") + ylab("Specificity") + 
    theme(plot.title = element_text(face = "bold"), 
          axis.title = element_text(face = "bold")), 
  
  # MDA Graph
  ggplot(IF_MDA_full, aes(factor(variables, 
                                 levels = rev(unique(IF_MDA_data_summary$variable)), 
                                 labels = rev(graph_labels)), 
                          log10(value))) + 
    geom_point(aes(color = variable)) + stat_summary(fun.y = "median", colour = "black", geom = "point", size = 2) + 
    coord_flip() + theme_bw() + ylab("Log10 MDA") + xlab("") +  ggtitle("B") + 
    scale_x_discrete(labels = rev(test2)) + 
    theme(plot.title = element_text(face = "bold"), 
          legend.position = "none", 
          axis.title = element_text(face = "bold"), 
          axis.text.y = element_text(size = 6)), 
  
  # Graph the IF model CRC data only
  filter(IF_red_follow_up_probability, diagnosis != "adenoma") %>% 
    ggplot(aes(factor(sampleType, 
                      levels = c("initial", "followup")), Yes, group = factor(EDRN))) + 
    geom_point(aes(color=factor(disease_free, 
                                levels = c("n", "y", "unknown"))), size = 2) + 
    geom_line(aes(color=factor(disease_free, 
                               levels = c("n", "y", "unknown"))), alpha = 0.6) + 
    scale_color_manual(name = "Cancer Free\n on Follow Up", 
                       label = c("No", "Yes", "Unknown"),  
                       values = wes_palette("GrandBudapest")) + 
    scale_x_discrete(breaks = c("initial", "followup"), 
                     labels = c("Initial", "Follow Up")) + 
    coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(aes(yintercept = 0.5), linetype = 2) + 
    ggtitle("C") + ylab("Initial Postive Probability") + 
    xlab("") + theme_bw() + theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = c(0.25, 0.20), 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the IF model adenoma data only
  filter(IF_red_follow_up_probability, diagnosis == "adenoma") %>%
    ggplot(aes(factor(sampleType, levels = c("initial", "followup")), 
               Yes, group = factor(EDRN))) + 
    geom_point(aes(color=factor(Dx_Bin))) + 
    geom_line(aes(color = factor(Dx_Bin))) + 
    scale_color_manual(name = "Polyp Type", values = c("cyan", "blue"), 
                       breaks = c("adenoma", "adv_adenoma"), 
                       labels = c("Adenoma", "SRN")) + 
    scale_x_discrete(
      breaks = c("initial", "followup"), 
      labels = c("Initial", "Follow Up")) + 
    coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(aes(yintercept = 0.5), linetype = 2) + 
    ggtitle("D") + ylab("Initial Postive Probability") + xlab("") + theme_bw() + 
    theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = c(0.20, 0.15), 
      plot.title = element_text(face="bold", hjust = 0)))


# Save figures and write necessary tables
ggsave(file = "results/figures/Figure4.tiff", initial_model, 
       width=8.5, height = 11, dpi = 300)






 








