## Run thetayc and Beta diversity analysis
## Main goal is to test whether theres is a beta diversity difference between inital and followups
## Marc Sze

# Load required dependencies and libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "vegan"))

# Load needed data
thetaCompTotal <- dissplit('data/process/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist',split=F, meta = F)
metaI <- read.csv("results/tables/mod_metadata/metaI_final.csv", stringsAsFactors = F, header = T)
metaF <- read.csv("results/tables/mod_metadata/metaF_final.csv", stringsAsFactors = F, header = T)
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)

difference_table_treatment <- pickDistanceValues(
  as.character(metaF$initial), as.character(metaF$followUp), thetaCompTotal, metaF, 
  c("fit_result", "fit_followUp", "Dx_Bin", "dx")) %>% 
  mutate(fit_difference = fit_result - fit_followUp)


pValueList <- list(thetaDiff = wilcox.test(distance ~ dx, data = difference_table_treatment)$p.value, 
                   fitDiff = wilcox.test(fit_difference ~ dx, data = difference_table_treatment)$p.value)

#Difference between initial and follow up for thetayc and fit broken down by adenoma and cancer
diff_adn_v_crc <- grid.arrange(
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
    annotate("text", label = paste("P-value = ", round(pValueList[["thetaDiff"]], digits = 2)), x = 1.5, y = 0.75), 
  
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
    annotate("text", label = paste("P-value = ", format(round(pValueList[["fitDiff"]], digits = 7))), x = 1.5, y = 0)
)


ggsave(file = "results/figures/adn_v_crc_diff.tiff", diff_adn_v_crc, width=8, height = 8, dpi = 300)


## Train Set bacterial community composition

test_set_theta_init_follow_dist <- thetaCompTotal[rownames(thetaCompTotal) %in% as.character(metaI$sample), 
                                                  colnames(thetaCompTotal) %in% as.character(metaI$sample)]

breakDown_samples <- metaI$dx[as.character(metaI$sample) %in% colnames(test_set_theta_init_follow_dist)]

set.seed(050416)
pValueList[["bdiver_Train"]] <- adonis(as.dist(test_set_theta_init_follow_dist) ~ factor(breakDown_samples))

set.seed(050416)
thetayc.mds.List <- list(bdiver_Train = metaMDS(as.dist(test_set_theta_init_follow_dist)) %>% scores() %>% as.data.frame() %>% 
  mutate(samples = factor(breakDown_samples)))

train_bdiver <- ggplot(data = thetayc.mds.List[["bdiver_Train"]], aes(x=NMDS1, y=NMDS2)) + geom_point(aes(color=samples)) + 
  theme_bw() + coord_equal() + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon") + 
  annotate("text", label = paste("PERMANOVA = ", round(pValueList[["bdiver_Train"]]$aov.tab$`Pr(>F)`[1], 
                                                       digits = 3)), x = -0.5, y = -0.3)


ggsave(file = "results/figures/Train_Bdiversity.tiff", train_bdiver, width=8, height = 8, dpi = 300)

## Follow up versus initial 

breakDown_samples <- c(rep("initial", length(rownames(metaF))), rep("follow_up", length(rownames(metaF))))

theta_init_follow_dist <- pickDistanceValues(
  as.character(c(metaF$initial, metaF$followUp)), as.character(c(metaF$initial, metaF$followUp)), 
  thetaCompTotal, withMeta = FALSE)

set.seed(050416)
pValueList[["bdiver_Test_IF"]] <- adonis(as.dist(theta_init_follow_dist) ~ factor(breakDown_samples))

set.seed(050416)
thetayc.mds.List[["bdiver_Test_IF"]] <- metaMDS(as.dist(theta_init_follow_dist)) %>% scores() %>% as.data.frame() %>% 
  mutate(samples = factor(breakDown_samples))

diff_init_follow <- ggplot(data = thetayc.mds.List[["bdiver_Test_IF"]], aes(x=NMDS1, y=NMDS2)) + geom_point(aes(color=samples)) + 
  theme_bw() + coord_equal() + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon") + 
  scale_colour_manual(values = c("red", "green")) + 
  annotate("text", label = paste("PERMANOVA = ", round(pValueList[["bdiver_Test_IF"]]$aov.tab$`Pr(>F)`[1], 
                                                       digits = 3)), x = 0.2, y = -0.5)

ggsave(file = "results/figures/init_v_follow_diff.tiff", diff_init_follow, width=8, height = 8, dpi = 300)


## Follow and Initial versus normal in train set

breakDown_samples <- c(rep("initial", length(rownames(metaF))), rep("follow_up", length(rownames(metaF))), 
                       rep("normal", length(metaI$sample[metaI$dx == "normal"])))

theta_init_follow_dist <- pickDistanceValues(
  as.character(c(metaF$initial, metaF$followUp, metaI$sample[metaI$dx == "normal"])), 
  as.character(c(metaF$initial, metaF$followUp, metaI$sample[metaI$dx == "normal"])), 
  thetaCompTotal, withMeta = FALSE)

set.seed(050416)
pValueList[["bdiver_Test_IFN"]] <- adonis(as.dist(theta_init_follow_dist) ~ factor(breakDown_samples))

set.seed(050416)
thetayc.mds.List[["bdiver_Test_IFN"]] <- metaMDS(as.dist(theta_init_follow_dist)) %>% scores() %>% as.data.frame() %>% 
  mutate(samples = factor(breakDown_samples))

diff_init_follow <- ggplot(data = thetayc.mds.List[["bdiver_Test_IFN"]] , aes(x=NMDS1, y=NMDS2)) + geom_point(aes(color=samples)) + 
  theme_bw() + coord_equal() + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon") + 
  annotate("text", label = paste("PERMANOVA = ", round(pValueList[["bdiver_Test_IFN"]]$aov.tab$`Pr(>F)`[1], 
                                                       digits = 3)), x = -0.3, y = -0.59)

ggsave(file = "results/figures/init_v_follow_norm_diff.tiff", diff_init_follow, width=8, height = 8, dpi = 300)


## Compare crc of init versus follow

crc_only_theta_init_follow_dist <- pickDistanceValues(
  as.character(c(metaF$initial[metaF$dx != "adenoma"], metaF$followUp[metaF$dx != "adenoma"])), 
  as.character(c(metaF$initial[metaF$dx != "adenoma"], metaF$followUp[metaF$dx != "adenoma"])), 
  thetaCompTotal, withMeta = FALSE)

breakDown_samples <- c(rep("initial", length(metaF$dx[metaF$dx != "adenoma"])), 
                       rep("follow_up", length(metaF$dx[metaF$dx != "adenoma"])))

set.seed(050416)
pValueList[["bdiver_Test_crc_IF"]] <- adonis(as.dist(crc_only_theta_init_follow_dist) ~ factor(breakDown_samples))

set.seed(050416)
thetayc.mds.List[["bdiver_Test_crc_IF"]] <- metaMDS(as.dist(crc_only_theta_init_follow_dist)) %>% 
  scores() %>% as.data.frame() %>% 
  mutate(samples = factor(breakDown_samples))

ggplot(data = thetayc.mds.List[["bdiver_Test_crc_IF"]], aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(color=samples)) + theme_bw() + coord_equal() + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon")



## Compare crc of init versus follow to normal 

crc_only_theta_init_follow_dist <- pickDistanceValues(
  as.character(c(metaF$initial[metaF$dx != "adenoma"], metaF$followUp[metaF$dx != "adenoma"], 
                 metaI$sample[metaI$dx == "normal"])), 
  as.character(c(metaF$initial[metaF$dx != "adenoma"], metaF$followUp[metaF$dx != "adenoma"], 
                 metaI$sample[metaI$dx == "normal"])), thetaCompTotal, withMeta = FALSE)

breakDown_samples <- c(rep("initial", length(metaF$dx[metaF$dx != "adenoma"])), 
                       rep("follow_up", length(metaF$dx[metaF$dx != "adenoma"])), 
                       rep("normal", length(metaI$sample[metaI$dx == "normal"])))

set.seed(050416)
pValueList[["bdiver_Test_crc_IFN"]] <- adonis(as.dist(crc_only_theta_init_follow_dist) ~ factor(breakDown_samples))

set.seed(050416)
thetayc.mds.List[["bdiver_Test_crc_IFN"]] <- metaMDS(as.dist(crc_only_theta_init_follow_dist)) %>% 
  scores() %>% as.data.frame() %>% 
  mutate(samples = factor(breakDown_samples))

ggplot(data = thetayc.mds.List[["bdiver_Test_crc_IFN"]], aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(color=samples)) + theme_bw() + coord_equal() + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon")


## Compare crc of init versus norm

crc_only_theta_init_follow_dist <- pickDistanceValues(
  as.character(c(metaF$initial[metaF$dx != "adenoma"], metaI$sample[metaI$dx == "normal"])), 
  as.character(c(metaF$initial[metaF$dx != "adenoma"], metaI$sample[metaI$dx == "normal"])), 
  thetaCompTotal, withMeta = FALSE)

breakDown_samples <- c(rep("initial", length(metaF$dx[metaF$dx != "adenoma"])), 
                       rep("normal", length(metaI$sample[metaI$dx == "normal"])))

set.seed(050416)
pValueList[["bdiver_Test_crc_IN"]] <- adonis(as.dist(crc_only_theta_init_follow_dist) ~ factor(breakDown_samples))

set.seed(050416)
thetayc.mds.List[["bdiver_Test_crc_IN"]] <- metaMDS(as.dist(crc_only_theta_init_follow_dist)) %>% 
  scores() %>% as.data.frame() %>% 
  mutate(samples = factor(breakDown_samples))

ggplot(data = thetayc.mds.List[["bdiver_Test_crc_IN"]], aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(color=samples)) + theme_bw() + coord_equal() + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon")


## Compare adenoma of init versus follow

polyp_only_theta_init_follow_dist <- pickDistanceValues(
  as.character(c(metaF$initial[metaF$dx == "adenoma"], metaF$followUp[metaF$dx == "adenoma"])), 
  as.character(c(metaF$initial[metaF$dx == "adenoma"], metaF$followUp[metaF$dx == "adenoma"])), 
  thetaCompTotal, withMeta = FALSE)

breakDown_samples <- c(rep("initial", length(metaF$dx[metaF$dx == "adenoma"])), 
                       rep("follow_up", length(metaF$dx[metaF$dx == "adenoma"])))

set.seed(050416)
pValueList[["bdiver_Test_adn_IF"]] <- adonis(as.dist(polyp_only_theta_init_follow_dist) ~ factor(breakDown_samples))

set.seed(050416)
thetayc.mds.List[["bdiver_Test_adn_IF"]] <- metaMDS(as.dist(polyp_only_theta_init_follow_dist), trace = 0) %>% 
  scores() %>% as.data.frame() %>%  mutate(samples = factor(breakDown_samples))

ggplot(data = thetayc.mds.List[["bdiver_Test_adn_IF"]], aes(x=NMDS1, y=NMDS2)) + geom_point(aes(color=samples)) + theme_bw() + coord_equal() + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon") + 
  annotate("text", label = paste("PERMANOVA = ", round(pValueList[["bdiver_Test_adn_IF"]]$aov.tab$`Pr(>F)`[1], 
                                                       digits = 3)), x = 0.25, y = -0.50)


## Compare adenoma init, follow versus normal

polyp_only_theta_init_follow_dist <- pickDistanceValues(
  as.character(c(metaF$initial[metaF$dx == "adenoma"], metaF$followUp[metaF$dx == "adenoma"], 
                 metaI$sample[metaI$dx == "normal"])), 
  as.character(c(metaF$initial[metaF$dx == "adenoma"], metaF$followUp[metaF$dx == "adenoma"], 
                 metaI$sample[metaI$dx == "normal"])), thetaCompTotal, withMeta = FALSE)

breakDown_samples <- c(rep("initial", length(metaF$dx[metaF$dx == "adenoma"])), 
                       rep("follow_up", length(metaF$dx[metaF$dx == "adenoma"])), 
                       rep("normal", length(metaI$sample[metaI$dx == "normal"])))

set.seed(050416)
pValueList[["bdiver_Test_adn_IFN"]] <- adonis(as.dist(polyp_only_theta_init_follow_dist) ~ factor(breakDown_samples))

set.seed(050416)
thetayc.mds.List[["bdiver_Test_adn_IFN"]] <- metaMDS(as.dist(polyp_only_theta_init_follow_dist)) %>% scores() %>% 
  as.data.frame() %>% mutate(samples = factor(breakDown_samples))

ggplot(data = thetayc.mds.List[["bdiver_Test_adn_IFN"]], aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(color=samples)) + theme_bw() + coord_equal() + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon")


## Compare adenoma init versus normal

polyp_only_theta_init_follow_dist <- pickDistanceValues(
  as.character(c(metaF$initial[metaF$dx == "adenoma"], metaI$sample[metaI$dx == "normal"])), 
  as.character(c(metaF$initial[metaF$dx == "adenoma"], metaI$sample[metaI$dx == "normal"])), 
  thetaCompTotal, withMeta = FALSE)

breakDown_samples <- c(rep("initial", length(metaF$dx[metaF$dx == "adenoma"])), 
                       rep("normal", length(metaI$sample[metaI$dx == "normal"])))

set.seed(050416)
pValueList[["bdiver_Test_adn_IN"]] <- adonis(as.dist(polyp_only_theta_init_follow_dist) ~ factor(breakDown_samples))

set.seed(050416)
thetayc.mds.List[["bdiver_Test_adn_IN"]] <- metaMDS(as.dist(polyp_only_theta_init_follow_dist)) %>% scores() %>% 
  as.data.frame() %>% mutate(samples = factor(breakDown_samples))

ggplot(data = thetayc.mds.List[["bdiver_Test_adn_IN"]], aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(color=samples)) + theme_bw() + coord_equal() + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon")


## Compare adenoma final versus normal

polyp_only_theta_init_follow_dist <- pickDistanceValues(
  as.character(c(metaF$followUp[metaF$dx == "adenoma"], metaI$sample[metaI$dx == "normal"])), 
  as.character(c(metaF$followUp[metaF$dx == "adenoma"], metaI$sample[metaI$dx == "normal"])), 
  thetaCompTotal, withMeta = FALSE)

breakDown_samples <- c(rep("follow_up", length(metaF$dx[metaF$dx == "adenoma"])), 
                       rep("normal", length(metaI$sample[metaI$dx == "normal"])))

set.seed(050416)
pValueList[["bdiver_Test_adn_FN"]] <- adonis(as.dist(polyp_only_theta_init_follow_dist) ~ factor(breakDown_samples))

set.seed(050416)
thetayc.mds.List[["bdiver_Test_adn_FN"]] <- metaMDS(as.dist(polyp_only_theta_init_follow_dist)) %>% scores() %>% 
  as.data.frame() %>% mutate(samples = factor(breakDown_samples))

ggplot(data = thetayc.mds.List[["bdiver_Test_adn_FN"]], aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(color=samples)) + theme_bw() + coord_equal() + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon")








