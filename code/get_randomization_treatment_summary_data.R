### Pull summary stats from Randomization
### Treatment models
### Marc Sze

# Load needed libraries 
library(dplyr)


# Read in needed data
treatment_groups <- c("adn", "SRN", "CRC")

treat_rand <- sapply(treatment_groups, 
                     function(x) read.csv(paste("data/process/tables/", x, 
                                                "_randomization_treatment_ROC_model_summary.csv", sep = ""), 
                                          header = T, stringsAsFactors = F), simplify = F)

treat_actual <- sapply(c("adn", "srn", "crc"), 
                       function(x) read.csv(paste("data/process/tables/", x, 
                                                  "_treatment_ROC_model_summary.csv", sep = ""), 
                                            header = T, stringsAsFactors = F), simplify = F)


# Function to generate model stats
get_summary_stats <- function(i, column_of_int, dataList){
  
  mod_column <- ifelse(dataList[[i]][, column_of_int] < 0.5, 
                       1-dataList[[i]][, column_of_int], dataList[[i]][, column_of_int])
  
  average <- mean(mod_column)
  std_dv <- sd(mod_column)
  
  final_info <- c(average = average, std_dv = std_dv, dis_group = i)
  
  return(final_info)
  
}


# Function to compare the outputted AUCs from random to actual
get_pvalues <- function(i_rand, i_act, rand_dataList = treat_rand, 
                        act_dataList = treat_actual, 
                        column_of_int = "ROC"){
  
  random_set <- ifelse(rand_dataList[[i_rand]][, column_of_int] < 0.5, 
                       1-rand_dataList[[i_rand]][, column_of_int], 
                       rand_dataList[[i_rand]][, column_of_int])
  
  actual_set <- ifelse(act_dataList[[i_act]][, column_of_int] < 0.5, 
                       1-act_dataList[[i_act]][, column_of_int], 
                       act_dataList[[i_act]][, column_of_int])
  
  final_data <- c(pvalue = t.test(random_set, actual_set)$p.value, study = i_act)
  
  return(final_data)
  
}


# generate the data table with the summary stats
summary_data_rand <- t(as.data.frame.list(lapply(treatment_groups, 
                       function(x) get_summary_stats(x, "ROC", treat_rand))))

summary_data_act <- t(as.data.frame.list(lapply(c("adn", "srn", "crc"), 
                                                 function(x) get_summary_stats(x, "ROC", treat_actual))))


# Generate the data tables with the pvalues
pvalue_summary <- t(mapply(get_pvalues, treatment_groups, c("adn", "srn", "crc")))


# Write out the actual data
write.csv(summary_data_rand, "data/process/tables/randomization_treatment_AUC_summary.csv", row.names = F)

write.csv(summary_data_act, "data/process/tables/actual_treatment_AUC_summary.csv", row.names = F)

write.csv(pvalue_summary, "data/process/tables/treatment_comparison_AUC_summary.csv", row.names = F)


















