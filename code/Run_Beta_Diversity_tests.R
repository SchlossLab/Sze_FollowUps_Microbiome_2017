### NMDS and PERMANOVA Anlysis of Test (follow up samples)
### How does the overall community composition change between inital and follow up
## Marc Sze


# Load needed libraries and created functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", 
  "scales", "wesanderson", "vegan", "knitr"))


# Load needed data
thetaCompTotal <- dissplit('data/process/final.thetayc.0.03.lt.ave.dist',
  split=F, meta = F)
metaF <- read.csv("results/tables/mod_metadata/metaF_final.csv", 
  stringsAsFactors = F, header = T)


# Load in needed data from distance matrix in a specific order for both Adenoma and Carcinoma  
all_followUp_combined <- pickDistanceValues(
  as.character(c(metaF$initial[metaF$dx == "adenoma"], 
                 metaF$initial[metaF$dx == "cancer"], 
                 metaF$followUp[metaF$dx == "adenoma"], 
                 metaF$followUp[metaF$dx == "cancer"])), 
  as.character(c(metaF$initial[metaF$dx == "adenoma"], 
                 metaF$initial[metaF$dx == "cancer"], 
                 metaF$followUp[metaF$dx == "adenoma"], 
                 metaF$followUp[metaF$dx == "cancer"])), 
  thetaCompTotal, withMeta = FALSE)

# Create vector with group labels that match custom order of created distance matrix
breakDown_samples <- c(rep("initialA", 
  length(metaF$dx[metaF$dx == "adenoma"])),
                       rep("initialC", 
                        length(metaF$dx[metaF$dx == "cancer"])), 
                       rep("follow_upA", 
                        length(metaF$dx[metaF$dx == "adenoma"])), 
                       rep("follow_upC", 
                        length(metaF$dx[metaF$dx == "cancer"])))


# Create table to store the relevant PERMANOVA and F.Model values

beta_diver_summary <- as.data.frame(
  matrix(nrow = 3, ncol = 3, dimnames = list(
  rows = c("All", "Adenoma", "Carcinoma"), 
  cols = c("F_model", "R2", "PERMANOVA"))))


# Run the PERMANOVA to see if there is any difference between any of the groups together
set.seed(050416)
All_followups <- 
adonis(as.dist(all_followUp_combined) ~ factor(breakDown_samples))

# Populate 'All' component of table
beta_diver_summary['All', ] <- c(All_followups$aov.tab$F.Model[1], 
  All_followups$aov.tab$R2[1], All_followups$aov.tab$`Pr(>F)`[1])


# Load in needed Adenoma Data
polyp_only_theta_init_follow_dist <- pickDistanceValues(
  as.character(c(metaF$initial[metaF$dx == "adenoma"], 
    metaF$followUp[metaF$dx == "adenoma"])), 
  as.character(c(metaF$initial[metaF$dx == "adenoma"], 
    metaF$followUp[metaF$dx == "adenoma"])), 
  thetaCompTotal, withMeta = FALSE)

# Re populate the breakDown_samples vector to match Adenoma Data
breakDown_samples <- c(
  rep("initial", length(metaF$dx[metaF$dx == "adenoma"])), 
  rep("follow_up", length(metaF$dx[metaF$dx == "adenoma"])))

# Run adonis on adenoma samples only for initial and follow up
set.seed(050416)
adenoma_followups <- adonis(
  as.dist(polyp_only_theta_init_follow_dist) ~ factor(breakDown_samples))

# Populate 'Adenoma' component of table
beta_diver_summary['Adenoma', ] <- c(
  adenoma_followups$aov.tab$F.Model[1], adenoma_followups$aov.tab$R2[1], 
  adenoma_followups$aov.tab$`Pr(>F)`[1])

# Create NMDS axis of Adenoma only from thetayc distances
set.seed(050416)
thetayc.mds.List <- list(
  bdiver_Test_adn_IF = metaMDS(as.dist(polyp_only_theta_init_follow_dist), 
    trace = 0) %>% scores() %>% as.data.frame() %>%  
  mutate(samples = factor(breakDown_samples)))

# Load in needed Carcinoma data
crc_only_theta_init_follow_dist <- pickDistanceValues(
  as.character(c(metaF$initial[metaF$dx != "adenoma"], 
    metaF$followUp[metaF$dx != "adenoma"])), 
  as.character(c(metaF$initial[metaF$dx != "adenoma"], 
    metaF$followUp[metaF$dx != "adenoma"])), 
  thetaCompTotal, withMeta = FALSE)

# Re populate the breakDown_samples vector to match Carcinoma Data
breakDown_samples <- c(
  rep("initial", length(metaF$dx[metaF$dx != "adenoma"])), 
  rep("follow_up", length(metaF$dx[metaF$dx != "adenoma"])))

# Run adonis on carcinoma samples only for initial and follow up
set.seed(050416)
carcinoma_followups <- adonis(
  as.dist(crc_only_theta_init_follow_dist) ~ factor(breakDown_samples))

# Populate 'Carcinoma' component of table
beta_diver_summary['Carcinoma', ] <- c(
  carcinoma_followups$aov.tab$F.Model[1], carcinoma_followups$aov.tab$R2[1], 
  carcinoma_followups$aov.tab$`Pr(>F)`[1])

# Create NMDS axis of Carcinoma only from thetayc distances
set.seed(050416)
thetayc.mds.List[["bdiver_Test_crc_IF"]] <- metaMDS(
  as.dist(crc_only_theta_init_follow_dist), trace = 0) %>% 
  scores() %>% as.data.frame() %>% 
  mutate(samples = factor(breakDown_samples))

# Save table for later use
write.csv(beta_diver_summary, "results/tables/beta_diver_summary.csv")


# Plot the NMDS for Adenoma Only and Carcinoma Only
bdiver_init_follow_sep <- grid.arrange(
  # Adenoma Only initial and follow up
  ggplot(data = thetayc.mds.List[["bdiver_Test_adn_IF"]], 
    aes(x=NMDS1, y=NMDS2)) + geom_point(aes(color=samples)) + 
  theme_bw() + coord_equal() + ggtitle("A") + 
    stat_ellipse(aes(group = samples, color = samples, fill = samples), 
      alpha = 0.25, geom = "polygon") + 
    theme(plot.title = element_text(face = "bold", hjust = 0)) + 
    annotate("text", label = paste("PERMANOVA = ", 
      round(adenoma_followups$aov.tab$`Pr(>F)`[1], digits = 3)), 
    x = 0.20, y = -0.50, size = 2.5), 
  # Carcinoma Only initial and follow up
  ggplot(data = thetayc.mds.List[["bdiver_Test_crc_IF"]], 
    aes(x=NMDS1, y=NMDS2)) + 
    geom_point(aes(color=samples)) + theme_bw() + 
    coord_equal() + ggtitle("B") + 
    stat_ellipse(aes(group = samples, color = samples, fill = samples), 
      alpha = 0.25, geom = "polygon") + 
    theme(plot.title = element_text(face = "bold", hjust = 0)) + 
    annotate("text", label = paste("PERMANOVA = ", 
      round(carcinoma_followups$aov.tab$`Pr(>F)`[1], digits = 3)), 
             x = 0.25, y = -0.50, size = 2.5), ncol = 1, nrow = 2)


ggsave(file = "results/figures/Figure2.pdf", bdiver_init_follow_sep, 
  width=8, height = 8, dpi = 300)










