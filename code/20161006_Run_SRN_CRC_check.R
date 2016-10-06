### Assess if we gain anything from the 
### Base analysis without metadata or OTU groupings
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))


# Load in data tables needed
withFit_data <- read.csv("results/tables/withFIT.models.datatable.csv", header = T, stringsAsFactors = F)
woFit_data <- read.csv("results/tables/noFIT.models.datatable.csv", header = T, stringsAsFactors = F)
withFit_cutoffs <- read.csv("results/tables/withFIT.cutoffs.csv", header = T, stringsAsFactors = F)
woFit_cutoffs <- read.csv("results/tables/noFIT.cutoffs.csv", header = T, stringsAsFactors = F)


length(rownames(filter(withFit_data, L1 == "SRNlesion_rf_opt", time_point == "initial", detailed_diagnosis == "Cancer", 
               positive > withFit_cutoffs$cutoff[1]))) 









