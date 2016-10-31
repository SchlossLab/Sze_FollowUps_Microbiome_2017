### Prediction of Follow Ups Analysis
### How do the with and without FIT models do in correctly calling initial and follow up samples
### P-value table of comparison and Figure 4 graph
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", 
  "gridExtra", "scales", "wesanderson"))


# Read in data tables
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", 
  stringsAsFactors = F, header = T)
rf_prediction_summary <- read.csv(
  'results/tables/rf_prediction_summary.csv', header = T)
aucrf_model_cutoff <- read.csv('results/tables/aucrf_model_cutoffs.csv', 
  header = T, row.names = 1)

# create table to store the data
pvalue_summary <- as.data.frame(matrix(
  nrow = 4, ncol = 2, dimnames = c(
    list(rows = c("model_compare_i", "model_compare_f", 
      "wfit_if", "wofit_if"), cols = c("adenoma", "carcinoma")))))

# Create variable vectors to cycle through during the for loop
samples <- c("initial", "followup")
model_used <- c("wfit", "wofit")

# create 2by2 table, execute fisher exact test, and populate the table 
for(i in 1:length(rownames(pvalue_summary))){
  
  for(j in 1:length(colnames(pvalue_summary))){
    
    if (i <= 2){
      
      pvalue_summary[i, "carcinoma"] <- 
        makeANDfish_2by2(rf_prediction_summary, c("wfit", "wofit"), 
          c("Yes", "No"), aucrf_model_cutoff, model = TRUE, 
          model_sample_type = paste(samples[i]), 
                       model_type = NULL, remove_sample = "adenoma")
      
      pvalue_summary[i, "adenoma"] <- 
        makeANDfish_2by2(filter(rf_prediction_summary, 
          diagnosis != "N/D"), c("wfit", "wofit"), c("Yes", "No"), 
                       aucrf_model_cutoff, model = TRUE, 
                       model_sample_type = paste(samples[i]), 
                       model_type = NULL, remove_sample = "adenocarcinoma")
      
    } else if (i>2){
      
      # carcionma only
      pvalue_summary[i, "carcinoma"] <- 
        makeANDfish_2by2(rf_prediction_summary, c("initial", "followup"), 
          c("Yes", "No"), aucrf_model_cutoff, model = FALSE, 
          model_sample_type = NULL, model_type = paste(model_used[i-2]), 
          remove_sample = "adenoma")
      
      # adenoma only
      pvalue_summary[i, "adenoma"] <- 
        makeANDfish_2by2(filter(rf_prediction_summary, diagnosis != "N/D"), 
                         c("wfit", "wofit"), c("Yes", "No"), 
                       aucrf_model_cutoff, model = FALSE, 
                       model_sample_type = NULL, 
                       model_type = paste(model_used[i-2]), 
                       remove_sample = "adenocarcinoma")
    }
  }
}

# Remove temporary variables
rm(i, j, model_used, samples)

# create table to hold wilcoxson paired tests pvalues

wilcox_pvalue_summary <- as.data.frame(matrix(
  nrow = 3, ncol = 2, dimnames = list(
    rows = c("lesion", "adenoma_only", "carcinoma_only"), 
    cols = c("wfit", "wofit"))))

# Set up variable vector
lesion_type <- c("all", "adenocarcinoma", "adenoma")
filter_diagnosis <- c("none", "N/D", "none")

for(i in 1:length(lesion_type)){
  
  wilcox_pvalue_summary[i, "wfit"] <- wilcox.test(
    filter(rf_prediction_summary, 
      model == "wfit" & 
      sample_type == "initial" & 
      diagnosis != paste(lesion_type[i]) & 
      diagnosis != paste(filter_diagnosis[i]))[, "postive_probability"], 
    filter(rf_prediction_summary, 
      model == "wfit" & 
      sample_type == "followup" & 
      diagnosis != paste(lesion_type[i]) & 
      diagnosis != paste(filter_diagnosis[i]))[, "postive_probability"], 
    paired = TRUE)$p.value
  
  wilcox_pvalue_summary[i, "wofit"] <- wilcox.test(
    filter(rf_prediction_summary, model == "wofit" & 
      sample_type == "initial" & 
      diagnosis != paste(lesion_type[i]) & 
      diagnosis != paste(filter_diagnosis[i]))[, "postive_probability"], 
    filter(rf_prediction_summary, model == "wofit" & 
      sample_type == "followup" & 
      diagnosis != paste(lesion_type[i]) & 
      diagnosis != paste(filter_diagnosis[i]))[, "postive_probability"], 
    paired = TRUE)$p.value
}


# Create Figure 4 
# Visual summery of the pvalues obtained

accuracy_plot <- grid.arrange(
  # Graph the carcinoma wFIT data only
  filter(rf_prediction_summary, 
    diagnosis != "adenoma" & model == "wfit") %>%
    ggplot(aes(factor(sample_type, levels = c("initial", "followup")), 
      postive_probability, group = factor(
        rep(filter(good_metaf, Diagnosis != "adenoma")[, "EDRN"], 2)))) + 
    geom_point(aes(color=factor(
      rep(filter(good_metaf, Diagnosis != "adenoma")[, "Disease_Free"], 2), 
                                levels = c("n", "y", "unknown")))) + 
    geom_line(linetype = 2, alpha = 0.3) + 
    scale_color_manual(name = "Cancer\nFree", 
      values = wes_palette("GrandBudapest")) + 
    scale_x_discrete(breaks = c("initial", "followup"), 
      labels = c("Initial", "Follow Up")) + 
    coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = aucrf_model_cutoff, aes(yintercept = wfit), 
      linetype = 2) + ggtitle("A") + ylab("Postive Probability") + 
    xlab("") + theme_bw() + theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = "none", 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the adenoma wFIT data only
  filter(rf_prediction_summary, 
    diagnosis == "adenoma" & model == "wfit") %>%
    ggplot(aes(factor(sample_type, levels = c("initial", "followup")), 
      postive_probability, group = factor(
        rep(filter(good_metaf, Diagnosis == "adenoma")[, "EDRN"], 2)))) + 
    geom_point(aes(color=factor(
      rep(filter(good_metaf, Diagnosis == "adenoma")[, "Dx_Bin"], 2)))) + 
    geom_line(aes(color = factor(
      rep(filter(good_metaf, Diagnosis == "adenoma")[, "Dx_Bin"], 2)))) + 
    scale_color_manual(name = "Polyp Type", values = c("cyan", "blue"), 
                       breaks = c("Adenoma", "adv Adenoma"), 
                       labels = c("Adenoma", "SRN")) + 
    scale_x_discrete(
      breaks = c("initial", "followup"), 
      labels = c("Initial", "Follow Up")) + 
    coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = aucrf_model_cutoff, 
      aes(yintercept = wfit), linetype = 2) + 
    ggtitle("C") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = "none", 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the carcinoma woFIT data only
  filter(rf_prediction_summary,  
    diagnosis != "adenoma" & model == "wofit") %>%
    ggplot(aes(factor(sample_type, levels = c("initial", "followup")), 
      postive_probability, group = factor(
        rep(filter(good_metaf, Diagnosis != "adenoma")[, "EDRN"], 2)))) + 
    geom_point(aes(color=factor(
      rep(filter(good_metaf, Diagnosis != "adenoma")[, "Disease_Free"], 2), 
                                levels = c("n", "y", "unknown")))) + 
    geom_line(linetype = 2, alpha = 0.3) + 
    scale_color_manual(
      name = "Cancer Free", 
      label = c("No", "Yes", "Unknown"), 
      values = wes_palette("GrandBudapest")) + 
    scale_x_discrete(
      breaks = c("initial", "followup"), 
      labels = c("Initial", "Follow Up")) + 
    coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = aucrf_model_cutoff, 
      aes(yintercept = wofit), linetype = 2) + 
    ggtitle("B") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = c(0.20, 0.15), 
      plot.title = element_text(face="bold", hjust = 0)), 
  
  # Graph the adenoma woFIT data only
  filter(rf_prediction_summary, 
    diagnosis == "adenoma" & model == "wofit") %>%
    ggplot(aes(factor(sample_type, 
      levels = c("initial", "followup")), postive_probability, 
    group = factor(
      rep(filter(good_metaf, Diagnosis == "adenoma")[, "EDRN"], 2)))) + 
    geom_point(aes(color=factor(
      rep(filter(good_metaf, Diagnosis == "adenoma")[, "Dx_Bin"], 2)))) + 
    geom_line(aes(color = factor(
      rep(filter(good_metaf, Diagnosis == "adenoma")[, "Dx_Bin"], 2)))) + 
    scale_color_manual(
      name = "Polyp Type", 
      values = c("cyan", "blue"), 
      breaks = c("Adenoma", "adv Adenoma"), 
      labels = c("Adenoma", "SRN")) + 
    scale_x_discrete(
      breaks = c("initial", "followup"), 
      labels = c("Initial", "Follow Up")) + 
    coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = aucrf_model_cutoff, 
      aes(yintercept = wofit), linetype = 2) + 
    ggtitle("D") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(
      axis.title = element_text(face="bold"), 
      legend.title = element_text(face="bold"), 
      legend.position = c(0.5, 0.5), 
      plot.title = element_text(face="bold", hjust = 0))
)


# Save figures and write necessary tables

ggsave(file = "results/figures/Figure4.pdf", accuracy_plot, 
       width=8.5, height = 11, dpi = 300)

write.csv(pvalue_summary, 
  "results/tables/IF_model_fisher_pvalue_summary.csv")
write.csv(wilcox_pvalue_summary, 
  "results/tables/IF_model_wicox_paired_pvalue_summary.csv")

