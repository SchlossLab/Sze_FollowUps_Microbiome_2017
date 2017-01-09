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
lesion_model_roc <- read.csv("results/tables/test_data_roc.csv", header = T)
reduced_lesion_model_roc <- read.csv("results/tables/reduced_test_data_roc.csv", header = T)
lesion_auc_info <- read.csv("results/tables/auc_summary.csv", 
  header = T, row.names = 1)


IF_model_roc <- read.csv("results/tables/IF_test_data_roc.csv", header = T)
reduced_IF_model_roc <- read.csv("results/tables/reduced_IF_test_data_roc.csv", header = T)

roc_curve_graph <- grid.arrange(
  # Graph showing all variable and reduced variable lesion model
  test <- filter(lesion_model_roc, run != "middle_roc") %>% 
    ggplot(aes(x = sensitivities, y = specificities)) + 
    geom_polygon(alpha = 0.25) + 
    geom_line(aes(group = run), size = 1.25, color = "black", alpha = 0.5) + 
    geom_polygon(data = filter(reduced_lesion_model_roc, run != "full_roc"), alpha = 0.25, fill = "red") + 
    geom_line(data = filter(reduced_lesion_model_roc, run != "full_roc"), 
              aes(group = run), size = 1.25, color = "red", alpha = 0.5) + 
    geom_line(data = filter(reduced_lesion_model_roc, run == "full_roc"), 
              size = 1.5, color = "cornflowerblue") + 
    scale_x_continuous(trans = "reverse") + theme_bw() + geom_abline(intercept = 1, linetype = 2, size = 1) + 
    ggtitle("A") + theme(plot.title = element_text(face= "bold")), 
  
  # Graph showing all variable and reduced variable IF model
  ggplot(IF_model_roc, aes(x = sensitivities, y = specificities)) + 
    geom_polygon(alpha = 0.25) + 
    geom_line(aes(group = run), size = 1.25, color = "black", alpha = 0.5) + 
    geom_polygon(data = filter(reduced_IF_model_roc, run != "full_roc"), alpha = 0.25, fill = "red") + 
    geom_line(data = filter(reduced_IF_model_roc, run != "full_roc"), 
              aes(group = run), size = 1.25, color = "red", alpha = 0.5) + 
    geom_line(data = filter(reduced_IF_model_roc, run == "full_roc"), 
              size = 1.5, color = "cornflowerblue") + 
    scale_x_continuous(trans = "reverse") + theme_bw() + geom_abline(intercept = 1, linetype = 2, size = 1) + 
    ggtitle("B") + theme(plot.title = element_text(face = "bold")), 
  nrow = 1, ncol = 2)
 

# Save graph image as a pdf
ggsave(file = "results/figures/Figure3.tiff", roc_curve_graph, device = "tiff", 
       width=8, height = 8, dpi = 300)



