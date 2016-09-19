## Interesting Graphs
## Creation of exploratory graphs that may end up in manuscript
## Marc Sze


source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))

#### Read in tables for Call Comparison Graph with and without Fit
withFIT <- read.csv("results/tables/withFIT.models.datatable.csv", header = T)
noFIT <- read.csv("results/tables/noFIT.models.datatable.csv", header = T)
withFIT_cutoffs <- read.csv("results/tables/withFIT.cutoffs.csv", header = T)
noFIT_cutoffs <- read.csv("results/tables/noFIT.cutoffs.csv", header = T)

# Remove the cancer model group
withFIT <- filter(withFIT, model != 'cancer')
withFIT_cutoffs <- filter(withFIT_cutoffs, model != 'cancer')


# Create facet names for graph and re order the levels for specific columns
Names_facet <- c('threeGroups' = "Classify N, A, C", 'SRNlesion' = "Classify SRN + Cancer", 'lesion' = "Classify Lesion")

noFIT$model <- factor(noFIT$model, levels = c("threeGroups", "SRNlesion", "lesion"))
withFIT$model <- factor(withFIT$model, levels = c("threeGroups", "SRNlesion", "lesion"))

noFIT$diseaseFree <- factor(noFIT$diseaseFree, levels = c("n", "y", "unknown"))
withFIT$diseaseFree <- factor(withFIT$diseaseFree, levels = c("n", "y", "unknown"))


# Adenoma graph only
grid.arrange(
  # With Fit in model performance with Adenoma (ALL data)
  filter(withFIT, diagnosis == "adenoma" & dataset == "All") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive)) + 
    geom_jitter(aes(color=detailed_diagnosis), width = 0.3) + 
    scale_color_manual(name = "Polyp Type", values = c("cyan", "blue"), 
                       breaks = c("Adenoma", "adv Adenoma"), labels = c("Adenoma", "SRN")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(withFIT_cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Adenomas with Fit (Train on All Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")), 
  
    # Without Fit in model performance with Adenoma (ALL data)
  filter(noFIT, diagnosis == "adenoma" & dataset == "All") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive)) + 
    geom_jitter(aes(color=detailed_diagnosis), width = 0.3) + 
    scale_color_manual(name = "Polyp Type", values = c("red", "darkred"), 
                       breaks = c("Adenoma", "adv Adenoma"), labels = c("Adenoma", "SRN")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(noFIT_cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Adenomas without Fit (Train on All Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold"))
)





# Cancer graph only
grid.arrange(
  # With Fit in model performance with Adenoma (ALL data)
  filter(withFIT, diagnosis == "adenocarcinoma" & dataset == "All") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive)) + 
    geom_jitter(aes(color=diseaseFree), width = 0.3) + 
    scale_color_manual(name = "Cancer\nFree", values = wes_palette("GrandBudapest")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(withFIT_cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Cancer with Fit (Train on All Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")), 
  
  # Without Fit in model performance with Adenoma (ALL data)
  filter(noFIT, diagnosis == "adenocarcinoma" & dataset == "All") %>%
    ggplot(aes(factor(time_point, levels = c("initial", "followup")), positive)) + 
    geom_jitter(aes(color=diseaseFree), width = 0.3) + 
    scale_color_manual(name = "Cancer\nFree", values = wes_palette("GrandBudapest")) + 
    facet_wrap(~model, labeller = as_labeller(Names_facet)) + coord_cartesian(ylim = c(0, 1)) + 
    geom_hline(data = filter(noFIT_cutoffs, dataset == "All"), aes(yintercept = cutoff), linetype = 2) + 
    ggtitle("Cancer without Fit (Train on All Data)") + ylab("Postive Probability") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold"))
)








