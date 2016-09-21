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

set.seed(050416)
thetayc.mds <- metaMDS(as.dist(theta_init_follow_dist)) %>% scores() %>% as.data.frame() %>% 
  mutate(samples = factor(breakDown_samples))

ggplot(data = thetayc.mds, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(color=samples)) + theme_bw() + coord_equal() + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon")


## Need to look at how FIT measurements compare to normal individuals

normalFit <- filter(metaI, dx == "normal") %>% select(fit_result)

initial_Fit_Pvalue <- wilcox.test(normalFit$fit_result, metaF$initial)$p.value
followUp_Fit_Pvalue <- wilcox.test(normalFit$fit_result, metaF$fit_followUp)$p.value

grid.arrange(
  # Initial Fit compared to normal Fit
  c(normalFit$fit_result, metaF$fit_result) %>% as.data.frame() %>% 
    mutate(sampleType = c(rep("normal", length(rownames(normalFit))), 
                                                           rep("initial", length(rownames(metaF))))) %>% 
    ggplot(aes(factor(sampleType, levels = c("normal", "initial")), log10(.+1), group = 1)) + 
    geom_jitter(aes(color=factor(sampleType, levels = c("normal", "initial"))), width = 0.5) + 
    scale_color_discrete(name = "Sample Type") + coord_cartesian(ylim = c(0, 4)) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + ggtitle("Fit") + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black")  + 
    ylab("Fit Result (ng/mL)") + xlab("") + theme_bw() + 
    theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), 
                                  axis.ticks.x = element_blank(), plot.title = element_text(face = "bold")) + 
    annotate("text", label = paste("P-value = ", format(round(initial_Fit_Pvalue, digits = 42), scientific = TRUE)), 
             x = 1.5, y = 3), 
  
  # Follow up microbiome sequences compared to normal
  c(normalFit$fit_result, metaF$fit_followUp) %>% as.data.frame() %>% 
    mutate(sampleType = c(rep("normal", length(rownames(normalFit))), 
                          rep("follow_up", length(rownames(metaF))))) %>% 
    ggplot(aes(factor(sampleType, levels = c("normal", "follow_up")), log10(.+1), group = 1)) + 
    geom_jitter(aes(color=factor(sampleType, levels = c("normal", "follow_up"))), width = 0.5) + 
    scale_color_discrete(name = "Sample Type") + coord_cartesian(ylim = c(0, 4)) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black") + 
    ylab("Fit Result (ng/mL)") + xlab("") + theme_bw() + 
    theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), 
                                  axis.ticks.x = element_blank()) + 
    annotate("text", label = paste("P-value = ", format(round(followUp_Fit_Pvalue, digits = 7))), x = 1.5, y = 3)
)


## Need to look at how important variables look at initial versus follow up

lesion_confirmed_vars_NOFIT <- read.csv('results/tables/lesion_confirmed_vars_NOFIT.csv', 
                                  header = TRUE, stringsAsFactors = F) %>% filter(otus != "fit_result")
SRNlesion_confirmed_vars_NOFIT <- read.csv('results/tables/SRNlesion_confirmed_vars_NOFIT.csv', 
                                     header = TRUE, stringsAsFactors = F) %>% filter(otus != "fit_result")
threeGroups_confirmed_vars_NOFIT <- read.csv('results/tables/threeGroups_confirmed_vars_NOFIT.csv', 
                                       header = TRUE, stringsAsFactors = F) %>% filter(otus != "fit_result")

shared <- read.delim('data/process/followUps.final.an.unique_list.0.03.subsample.0.03.filter.shared', 
                     header=T, sep='\t')

# Lesion Initial
pValue_test_data_lesion_NOFIT_init <- as.data.frame(cbind(c(metaF$initial, metaI$sample[metaI$dx == "normal"]), 
                  c(rep("lesion", length(rownames(metaF))), rep("normal", length(rownames(metaI[metaI$dx == "normal", ])))))) %>% 
  rename(Group = V1, sampleType = V2) %>% mutate(Group = as.integer(as.character((Group)))) %>% inner_join(shared, by = "Group") %>% 
  select(Group, sampleType, one_of(lesion_confirmed_vars_NOFIT$otus))


pValue_table_lesion_NOFIT_init <- getCorrectedPvalue(pValue_test_data_lesion_NOFIT_init, 
                                                lesion_confirmed_vars_NOFIT$otus, "sampleType") %>% 
  rename(initial = sampleType)

# Lesion Follow up
pValue_test_data_lesion_NOFIT_follow <- as.data.frame(cbind(c(metaF$followUp, metaI$sample[metaI$dx == "normal"]), 
                                                     c(rep("lesion", length(rownames(metaF))), rep("normal", length(rownames(metaI[metaI$dx == "normal", ])))))) %>% 
  rename(Group = V1, sampleType = V2) %>% mutate(Group = as.integer(as.character((Group)))) %>% inner_join(shared, by = "Group") %>% 
  select(Group, sampleType, one_of(lesion_confirmed_vars_NOFIT$otus))


pValue_table_lesion_NOFIT_follow <- getCorrectedPvalue(pValue_test_data_lesion_NOFIT_follow, 
                                                lesion_confirmed_vars_NOFIT$otus, "sampleType") %>% 
  rename(follow_up = sampleType)



# SRN lesion or Cancer Initial
pValue_test_data_SRNlesion_NOFIT_init <- as.data.frame(cbind(c(metaF$initial, metaI$sample[metaI$Dx_Bin == "Normal" | metaI$Dx_Bin == "Adenoma"]), 
                                                     c(rep("SRNlesionORCancer", length(rownames(metaF))), 
                                                       rep("normalorAdenoma", 
                                                           length(rownames(metaI[metaI$Dx_Bin == "Normal" | metaI$Dx_Bin == "Adenoma", ])))))) %>% 
  rename(Group = V1, sampleType = V2) %>% mutate(Group = as.integer(as.character((Group)))) %>% inner_join(shared, by = "Group") %>% 
  select(Group, sampleType, one_of(SRNlesion_confirmed_vars_NOFIT$otus))


pValue_table_SRNlesion_NOFIT_init <- getCorrectedPvalue(pValue_test_data_SRNlesion_NOFIT_init, 
                                                SRNlesion_confirmed_vars_NOFIT$otus, "sampleType") %>% 
  rename(initial = sampleType)

# SRN lesion or Cancer Follow up
pValue_test_data_SRNlesion_NOFIT_follow <- as.data.frame(cbind(c(metaF$followUp, 
                                                                 metaI$sample[metaI$Dx_Bin == "Normal" | metaI$Dx_Bin == "Adenoma"]), 
                                                             c(rep("SRNlesionORCancer", length(rownames(metaF))), 
                                                               rep("normalorAdenoma", 
                                                                   length(rownames(metaI[metaI$Dx_Bin == "Normal" | metaI$Dx_Bin == "Adenoma", ])))))) %>% 
  rename(Group = V1, sampleType = V2) %>% mutate(Group = as.integer(as.character((Group)))) %>% inner_join(shared, by = "Group") %>% 
  select(Group, sampleType, one_of(SRNlesion_confirmed_vars_NOFIT$otus))


pValue_table_SRNlesion_NOFIT_follow <- getCorrectedPvalue(pValue_test_data_SRNlesion_NOFIT_follow, 
                                                        SRNlesion_confirmed_vars_NOFIT$otus, "sampleType") %>% 
  rename(follow_up = sampleType)


# Graph

# Convert taxa table to a data frame with columns for each taxonomic division
tax <- read.delim('data/process/followUps.final.an.unique_list.0.03.cons.taxonomy', 
                  sep='\t', header=T, row.names=1)
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
# This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components
rm(tax)

lesion_selected_taxa <- tax_df[lesion_confirmed_vars_NOFIT[, 'otus'], ]
lesion_selected_labs <- createTaxaLabeller(lesion_selected_taxa)

pvalue_summary_table_lesion <- cbind(pValue_table_lesion_NOFIT_init, pValue_table_lesion_NOFIT_follow, lesion_selected_labs) %>% 
  rename(lowest_tax_ID = lesion_selected_labs)

write.csv(pvalue_summary_table_lesion, "results/tables/lesion_Imp_vars_NOFIT_pvalue_summary.csv")

SRNlesion_selected_taxa <- tax_df[SRNlesion_confirmed_vars_NOFIT[, 'otus'], ]
SRNlesion_selected_labs <- createTaxaLabeller(SRNlesion_selected_taxa)

pvalue_summary_table_SRNlesion <- cbind(pValue_table_SRNlesion_NOFIT_init, pValue_table_SRNlesion_NOFIT_follow, 
                                     SRNlesion_selected_labs) %>% rename(lowest_tax_ID = SRNlesion_selected_labs)

write.csv(pvalue_summary_table_SRNlesion, "results/tables/SRNlesion_Imp_vars_NOFIT_pvalue_summary.csv")


specific_data_lesion_init <- melt(pValue_test_data_lesion_NOFIT_init, id = c("Group", "sampleType"))

specific_data_SRNlesion_init <- melt(pValue_test_data_SRNlesion_NOFIT_init, id = c("Group", "sampleType"))

specific_data_lesion_follow <- melt(pValue_test_data_lesion_NOFIT_follow, id = c("Group", "sampleType"))

specific_data_SRNlesion_follow <- melt(pValue_test_data_SRNlesion_NOFIT_follow, id = c("Group", "sampleType"))

# Graph for Lesion model only
grid.arrange(
  # Initial select microbiome sequences compared to normal
  filter(specific_data_lesion_init, 
         variable %in% rownames(pvalue_summary_table_lesion)[pvalue_summary_table_lesion$initial <= 0.05]) %>% 
    ggplot(aes(factor(sampleType, levels = c("normal", "lesion")), log10(value+1), group = 1)) + 
    geom_jitter(aes(color=factor(sampleType, levels = c("normal", "lesion"))), width = 0.5) + 
    scale_color_discrete(name = "Sample Type", breaks = c("normal", "lesion"), labels = c("normal", "initial")) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + ggtitle("Lesion Model") + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black") + coord_cartesian(ylim = c(0, 3.5)) + 
    facet_wrap(~variable, labeller = as_labeller(lesion_selected_labs)) + ylab("Log10 Subsampled Sequence Reads") + 
    xlab("") + theme_bw() + theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), 
                                  axis.ticks.x = element_blank(), plot.title = element_text(face = "bold")) + 
    annotate("text", label = paste("P-value = ", format(round(pvalue_summary_table_lesion$initial[
      pvalue_summary_table_lesion$initial <= 0.05], digits = 7))), x = 1.5, y = -0.1), 
  
  # Follow up microbiome sequences compared to normal
  filter(specific_data_lesion_follow, 
         variable %in% rownames(pvalue_summary_table_lesion)[pvalue_summary_table_lesion$initial <= 0.05]) %>% 
    ggplot(aes(factor(sampleType, levels = c("normal", "lesion")), log10(value+1), group = 1)) + 
    geom_jitter(aes(color=factor(sampleType, levels = c("normal", "lesion"))), width = 0.5) + 
    scale_color_discrete(name = "Sample Type", breaks = c("normal", "lesion"), labels = c("normal", "follow up")) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black") + coord_cartesian(ylim = c(0, 3.5)) + 
    facet_wrap(~variable, labeller = as_labeller(lesion_selected_labs)) + ylab("Log10 Subsampled Sequence Reads") + 
    xlab("") + theme_bw() + theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), 
                                  axis.ticks.x = element_blank()) + 
    annotate("text", label = paste("P-value = ", format(round(pvalue_summary_table_lesion$follow_up[
      pvalue_summary_table_lesion$initial <= 0.05], digits = 7))), x = 1.5, y = -0.1)
)


# Graph for SRN Lesion model only
grid.arrange(
  # Initial select microbiome sequences compared to normal
  filter(specific_data_SRNlesion_init, 
         variable %in% rownames(pvalue_summary_table_SRNlesion)[pvalue_summary_table_SRNlesion$initial <= 0.05]) %>% 
    ggplot(aes(factor(sampleType, levels = c("normalorAdenoma", "SRNlesionORCancer")), log10(value+1), group = 1)) + 
    geom_jitter(aes(color=factor(sampleType, levels = c("normalorAdenoma", "SRNlesionORCancer"))), width = 0.5) + 
    scale_color_discrete(name = "Sample Type", breaks = c("normalorAdenoma", "SRNlesionORCancer"), 
                         labels = c("Normal or Adenoma", "initial")) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + ggtitle("SRN + Cancer Model") + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black") + coord_cartesian(ylim = c(0, 3.5)) + 
    facet_wrap(~variable, labeller = as_labeller(SRNlesion_selected_labs)) + ylab("Log10 Subsampled Sequence Reads") + 
    xlab("") + theme_bw() + theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), 
                                  axis.ticks.x = element_blank(), plot.title = element_text(face = "bold")) + 
    annotate("text", label = paste("P-value = ", format(round(pvalue_summary_table_SRNlesion$initial[
      pvalue_summary_table_SRNlesion$initial <= 0.05], digits = 7))), x = 1.5, y = -0.1), 
  
  # Follow up microbiome sequences compared to normal
  filter(specific_data_SRNlesion_follow, 
         variable %in% rownames(pvalue_summary_table_SRNlesion)[pvalue_summary_table_SRNlesion$initial <= 0.05]) %>% 
    ggplot(aes(factor(sampleType, levels = c("normalorAdenoma", "SRNlesionORCancer")), log10(value+1), group = 1)) + 
    geom_jitter(aes(color=factor(sampleType, levels = c("normalorAdenoma", "SRNlesionORCancer"))), width = 0.5) + 
    scale_color_discrete(name = "Sample Type", breaks = c("normalorAdenoma", "SRNlesionORCancer"), 
                         labels = c("Normal or Adenoma", "follow up")) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black") + coord_cartesian(ylim = c(0, 3.5)) + 
    facet_wrap(~variable, labeller = as_labeller(SRNlesion_selected_labs)) + ylab("Log10 Subsampled Sequence Reads") + 
    xlab("") + theme_bw() + theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), 
                                  axis.ticks.x = element_blank()) + 
    annotate("text", label = paste("P-value = ", format(round(pvalue_summary_table_SRNlesion$follow_up[
      pvalue_summary_table_SRNlesion$initial <= 0.05], digits = 7))), x = 1.5, y = -0.1)
)


## Need to try the opposite and look at what RF and Borutat would pull out when looking solely at those with 
## follow up

## Need to look into more detail into the differences between these models


## Need to look at using random forest on follow up versus normal for what OTUs are different














