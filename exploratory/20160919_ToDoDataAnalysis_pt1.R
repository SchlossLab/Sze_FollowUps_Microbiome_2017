## A To Do List to follow up on stemming from original exploratory analysis
## Focus strictly on Lesion, SRNLesion, and three groups (Normal, Adenoma, Cancer) classifications
## Marc Sze

# Load required dependencies and libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "vegan"))


## Need to look at differences between initial and follow up between microbiome (thetayc) and fit.
# Does the microbiome show a difference in community versus the fit?

metaI <- read.delim('data/process/initials_metadata.tsv', header=T, sep='\t')
metaF <- read.delim('data/process/followUps_metadata.txt', header=T, sep='\t') %>% mutate(lesion = factor(NA, levels=c(0,1)))
thetaCompTotal <- dissplit('data/process/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist',split=F, meta = F)

difference_table_treatment <- pickDistanceValues(
  as.character(metaF$initial), as.character(metaF$followUp), thetaCompTotal, metaF, 
  c("fit_result", "fit_followUp", "Dx_Bin", "dx")) %>% 
  mutate(fit_difference = fit_result - fit_followUp)

thetaDiff_follow_pvalue <- wilcox.test(distance ~ dx, data = difference_table_treatment)$p.value

fitDiff_follow_pvalue <- wilcox.test(fit_difference ~ dx, data = difference_table_treatment)$p.value

grid.arrange(
  # Difference from bacterial community structure between adenoma and cancer
  ggplot(difference_table_treatment, aes(factor(dx, levels = c("adenoma", "cancer")), distance, group = 1)) + 
    geom_jitter(aes(color=Dx_Bin), width = 0.3, size = 4, alpha = 0.5) + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + 
    scale_color_manual(name = "Lesion Type", values = wes_palette("GrandBudapest"), 
                       breaks = c("Adenoma", "adv Adenoma", "Cancer"), labels = c("Adenoma", "SRN", "Cancer")) + 
    coord_cartesian(ylim = c(0, 1)) + ylab("Thetayc Distance") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")) + 
    annotate("text", label = paste("P-value = ", round(thetaDiff_follow_pvalue, digits = 2)), x = 1.5, y = 0.75), 
  
  # Difference from fit between adenoma and cancer
  ggplot(difference_table_treatment, aes(factor(dx, levels = c("adenoma", "cancer")), fit_difference, group = 1)) + 
    geom_jitter(aes(color=Dx_Bin), width = 0.3, size = 4, alpha = 0.5) + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + 
    scale_color_manual(name = "Lesion Type", values = wes_palette("GrandBudapest"), 
                       breaks = c("Adenoma", "adv Adenoma", "Cancer"), labels = c("Adenoma", "SRN", "Cancer")) + 
    ylab("Change in Fit from Initial to Follow Up") + xlab("") + theme_bw() + 
    theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), 
          title = element_text(face="bold")) + 
    annotate("text", label = paste("P-value = ", format(round(fitDiff_follow_pvalue, digits = 7))), x = 1.5, y = 0)
)


## Need to look at how the follow up samples look like versus their initial samples and then visualize versus
## the normal values

breakDown_samples <- c(rep("initial", length(rownames(metaF))), rep("follow_up", length(rownames(metaF))), 
                       rep("normal", length(metaI$sample[metaI$dx == "normal"])))

theta_init_follow_dist <- pickDistanceValues(
  as.character(c(metaF$initial, metaF$followUp, metaI$sample[metaI$dx == "normal"])), 
  as.character(c(metaF$initial, metaF$followUp, metaI$sample[metaI$dx == "normal"])), 
  thetaCompTotal, withMeta = FALSE)

set.seed(050416)
adonis(as.dist(theta_init_follow_dist) ~ factor(breakDown_samples))

thetayc.mds <- metaMDS(as.dist(theta_init_follow_dist)) %>% scores() %>% as.data.frame() %>% 
  mutate(samples = factor(breakDown_samples))

ggplot(data = thetayc.mds, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(color=samples)) + theme_bw() + coord_equal() + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon")


## Need to look at how FIT measurements compare to normal individuals

normalFit <- filter(metaI, dx == "normal") %>% select(fit_result)

initial_Fit_Pvalue <- wilcox.test(normalFit$fit_result, metaF$initial)$p.value
followUp_Fit_Pvalue <- wilcox.test(normalFit$fit_result, metaF$fit_followUp)$p.value


## Need to try the opposite and look at what RF and Borutat would pull out when looking solely at those with 
## follow up

## Need to look into more detail into the differences between these models

















