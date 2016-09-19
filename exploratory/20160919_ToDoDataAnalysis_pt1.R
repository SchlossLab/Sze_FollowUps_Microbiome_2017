## A To Do List to follow up on stemming from original exploratory analysis
## Focus strictly on Lesion, SRNLesion, and three groups (Normal, Adenoma, Cancer) classifications
## Marc Sze


## Need to look at differences between initial and follow up between microbiome (thetayc) and fit.
# Does the microbiome show a difference in community versus the fit?

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


## Need to try the opposite and look at what RF and Borutat would pull out when looking solely at those with 
## follow up

## Need to look into more detail into the differences between these models

## Need to look at how FIT measurements compare to a random sampling of normal individuals
# Probably need to bootstrap this















