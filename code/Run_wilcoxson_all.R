### All OTU Paired Wilcoxson Test
### Are there differences in any OTU in initial and follow up
### Uses the full lesion model data
## Marc Sze


#Load needed libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "scales", "wesanderson", "ggplot2"))

#Load needed data sets
lesion_data <- read.csv("results/tables/full_test_data.csv", header = T, row.names = 1)










