### Create an NMDS with all followup samples in it
### Are there any differences between Adenoma and cancer in the initial and follow up
## Marc Sze


# Load needed libraries and created functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "vegan", "knitr"))


# Load needed data
thetaCompTotal <- dissplit('data/process/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist',split=F, meta = F)
metaI <- read.csv("results/tables/mod_metadata/metaI_final.csv", stringsAsFactors = F, header = T)
metaF <- read.csv("results/tables/mod_metadata/metaF_final.csv", stringsAsFactors = F, header = T)
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)


# Load in needed data from distance matrix in a specific order  
all_followUp_combined <- pickDistanceValues(
  as.character(c(metaF$initial[metaF$dx == "adenoma"], 
                 metaF$initial[metaF$dx == "cancer"], 
                 metaF$followUp[metaF$dx == "adenoma"], 
                 metaF$followUp[metaF$dx == "cancer"])), 
  as.character(c(metaF$initial[metaF$dx == "adenoma"], 
                 metaF$initial[metaF$dx == "cancer"], 
                 metaF$followUp[metaF$dx == "adenoma"], 
                 metaF$followUp[metaF$dx == "cancer"])), thetaCompTotal, withMeta = FALSE)

# Create vector with group labels that match custom order of created distance matrix
breakDown_samples <- c(rep("initialA", length(metaF$dx[metaF$dx == "adenoma"])),
                       rep("initialC", length(metaF$dx[metaF$dx == "cancer"])), 
                       rep("follow_upA", length(metaF$dx[metaF$dx == "adenoma"])), 
                       rep("follow_upC", length(metaF$dx[metaF$dx == "cancer"])))


# Run the PERMANOVA
set.seed(050416)
pValueAdonis <- adonis(as.dist(all_followUp_combined) ~ factor(breakDown_samples))

# Create the NMDS axis 1 and 2
set.seed(050416)
thetayc.mds <- metaMDS(as.dist(all_followUp_combined), trace = 0) %>% 
                           scores() %>% as.data.frame() %>%  mutate(samples = factor(breakDown_samples))

# Graph the NMDS
ggplot(data = thetayc.mds, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(color=samples)) + 
  theme_bw() + coord_equal() + ggtitle("All Follow Ups Combined") + 
  stat_ellipse(aes(group = samples, color = samples, fill = samples), alpha = 0.25, geom = "polygon") + 
  annotate("text", label = paste("PERMANOVA = ", round(pValueAdonis$aov.tab$`Pr(>F)`[1], 
                                                       digits = 3)), x = 0.20, y = -0.53, size = 2.5)

