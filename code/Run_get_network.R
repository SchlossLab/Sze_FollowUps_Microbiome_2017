### Common OTUs Network
### Create data to be used in Network
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "reshape2", "scales"))

### Load in needed data tables
common_pred_data <- read.csv("results/tables/common_pred_data.csv", header = T, stringsAsFactors = F)

### Create data table
lesion_pos <- filter(common_pred_data, lesion == "Yes") %>% select(-lesion)
lesion_neg <- filter(common_pred_data, lesion == "No") %>% select(-lesion)

comparisons <- colnames(select(common_pred_data, -lesion))

### Create co-occurance percentages
neg_occurance <- get_co_occurance(lesion_neg, comparisons)
pos_occurance <- get_co_occurance(lesion_pos, comparisons)







