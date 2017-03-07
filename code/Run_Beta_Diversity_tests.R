### NMDS and PERMANOVA Anlysis of Test (follow up samples)
### How does the overall community composition change between inital and follow up
## Marc Sze


# Load needed libraries and created functions
source('code/functions.R')

loadLibs(c("dplyr", "vegan"))


# Load needed data
thetaCompTotal <- read.dist('data/process/final.thetayc.0.03.lt.ave.dist')
metaF <- read.csv("data/process/mod_metadata/metaF_final.csv", 
  stringsAsFactors = F, header = T)


# Load in needed data from distance matrix in a specific order for both Adenoma and Carcinoma  
dataList <- list(
  #choose all follow up distances
  all_followUp_dist = pickDistanceValues(
    as.character(c(metaF$initial[metaF$dx == "adenoma"], 
                   metaF$initial[metaF$dx == "cancer"], 
                   metaF$followUp[metaF$dx == "adenoma"], 
                   metaF$followUp[metaF$dx == "cancer"])), 
    as.character(c(metaF$initial[metaF$dx == "adenoma"], 
                   metaF$initial[metaF$dx == "cancer"], 
                   metaF$followUp[metaF$dx == "adenoma"], 
                   metaF$followUp[metaF$dx == "cancer"])), 
    thetaCompTotal, withMeta = FALSE), 
  #choose only adenoma distances
  polyp_dist = pickDistanceValues(
    as.character(c(metaF$initial[metaF$dx == "adenoma"], 
                   metaF$followUp[metaF$dx == "adenoma"])), 
    as.character(c(metaF$initial[metaF$dx == "adenoma"], 
                   metaF$followUp[metaF$dx == "adenoma"])), 
    thetaCompTotal, withMeta = FALSE), 
  #choose only crc distances
  crc_dist = pickDistanceValues(
    as.character(c(metaF$initial[metaF$dx != "adenoma"], 
                   metaF$followUp[metaF$dx != "adenoma"])), 
    as.character(c(metaF$initial[metaF$dx != "adenoma"], 
                   metaF$followUp[metaF$dx != "adenoma"])), 
    thetaCompTotal, withMeta = FALSE))
  
# Create vector with group labels that match custom order of created distance matrix
breakDown_samples <- list(
  #call vector for all follow up samples
  all = c(rep("initialA", length(metaF$dx[metaF$dx == "adenoma"])), 
          rep("initialC", length(metaF$dx[metaF$dx == "cancer"])), 
          rep("follow_upA", length(metaF$dx[metaF$dx == "adenoma"])), 
          rep("follow_upC", length(metaF$dx[metaF$dx == "cancer"]))), 
  #call vector for adenoma follow up samples
  adn = c(rep("initial", length(metaF$dx[metaF$dx == "adenoma"])), 
          rep("follow_up", length(metaF$dx[metaF$dx == "adenoma"]))), 
  #call vector for crc follow up samples
  crc = c(rep("initial", length(metaF$dx[metaF$dx != "adenoma"])), 
          rep("follow_up", length(metaF$dx[metaF$dx != "adenoma"]))))

# Create table to store the relevant PERMANOVA, R2, and F.Model values
beta_diver_summary <- as.data.frame(
  matrix(nrow = 3, ncol = 3, dimnames = list(
  rows = c("All", "Adenoma", "Carcinoma"), 
  cols = c("F_model", "R2", "PERMANOVA"))))

# Create list to store the relevant MDS data for graphs
thetayc.mds.list <- list(all = c(), adn = c(), crc = c())

# Analyze the bacterial community structure using PERMANOVA
for(i in 1:length(dataList)){
  
  # Run the PERMANOVA to see if there is any difference between any of the groups together
  set.seed(050416)
  tempAnalysisData <- 
    adonis(as.dist(dataList[[i]]) ~ factor(breakDown_samples[[i]]))
  
  # Populate 'All' component of table
  beta_diver_summary[i, ] <- c(tempAnalysisData$aov.tab$F.Model[1], 
                               tempAnalysisData$aov.tab$R2[1], tempAnalysisData$aov.tab$`Pr(>F)`[1])
  
  # Create NMDS axis values from thetayc distances and store for graphing
  set.seed(050416)
  thetayc.mds.list[[i]] <- metaMDS(as.dist(dataList[[i]]), trymax = 3000, trace = 0) %>% 
      scores() %>% as.data.frame() %>% mutate(samples = factor(breakDown_samples[[i]]))
  
}

# Save tables for later use
write.csv(beta_diver_summary, "data/process/tables/beta_diver_summary.csv")
write.csv(thetayc.mds.list[["adn"]], "data/process/tables/thetayc_adn_IF.csv", row.names = F)
write.csv(thetayc.mds.list[["crc"]], "data/process/tables/thetayc_crc_IF.csv", row.names = F)


