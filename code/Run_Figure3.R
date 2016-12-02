### Create Figure 3
### Graph to show AUC of full data with chosen model
### Graph to show how min, middle, and max performed on 20% test data
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", 
           "wesanderson", "caret", "pROC"))

### Load in needed data tables
test_roc <- read.csv("results/tables/test_data_roc.csv", header = T)
auc_info <- read.csv("results/tables/auc_summary.csv", header = T, row.names = 1)


### Graph the ROC curve of the middle model for all data

roc_curve_graph <- ggplot(test_roc, aes(sensitivities, specificities)) + 
  geom_line(aes(group = run, color = run), size = 1.5) + scale_x_continuous(trans = "reverse") + 
  theme_bw() + xlab("Sensitivity") + ylab("Specificity") + geom_abline(intercept = 1, linetype = 2) + 
  scale_color_discrete(name = "Train Set Models", 
                       breaks = c("best_roc", "middle_roc", "worse_roc"), 
                       labels = c(paste("Best Training Model\n(cvROC = ", 
                                        round(auc_info["best", "ROC_cv"], digits = 3), ")", sep = ""), 
                                  paste("Middle Training Model\n(cvROC = ", 
                                        round(auc_info["middle", "ROC_cv"], digits = 3), ")", sep=""), 
                                  paste("Worse Training Model\n(cvROC = ", 
                                        round(auc_info["worse", "ROC_cv"], digits = 3), ")", sep=""))) + 
  theme(axis.title = element_text(face = "bold"), legend.position = c(0.75, 0.4), 
        legend.title = element_text(face="bold"), legend.text = element_text(size = 10), 
        legend.key = element_rect(size = 15), legend.key.size = unit(2.5, 'lines'))


# Save graph image as a pdf
ggsave(file = "results/figures/Figure3.pdf", roc_curve_graph, 
       width=8, height = 8, dpi = 300)



