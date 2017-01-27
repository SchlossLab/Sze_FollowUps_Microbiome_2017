### List of Stored Functions
### CRC Follow up project
### Marc Sze


#This function loads and installs libraries if they are not already loaded by the user 
loadLibs <- function(deps){
  for (dep in deps){
    if (dep %in% installed.packages()[,"Package"] == FALSE){
      install.packages(as.character(dep), quiet=TRUE);
    }
    library(dep, verbose=FALSE, character.only=TRUE)
  }
}


# Run alpha diversity tests
# Order of paired samples has to already be set with the exact order matching 
# e.g. if you have 10 samples that are pairs 1 and 6 have to be the same person
get_alpha_pvalues <- function(data_table, alpha = TRUE, 
                              table_names = c("sobs", "shannon", "evenness"), 
                              numComp = length(table_names), multi = "BH"){
  
  #Set up output table
  alpha_pvalue_table <- as.data.frame(matrix(ncol = 2, nrow = numComp, dimnames = list(rows = table_names, 
                                                                                       cols = c("pvalue", "BH_adj_pvalue"))))
  #Get p-values of paired test
  alpha_pvalue_table[, "pvalue"] <- apply(data_table[, c("sobs", "shannon", "shannoneven")], 2, 
                                          function(x){
                                            wilcox.test(x[data_table$sampleType == "initial"], 
                                                        x[data_table$sampleType == "followups"], 
                                                        paired = TRUE)$p.value})
  
  # Get BH corrected p-values
  alpha_pvalue_table$BH_adj_pvalue <- p.adjust(alpha_pvalue_table$pvalue, method = multi)
  
  return(alpha_pvalue_table)
}



# This function creates a square distance matrix and splits it based 
# on specified criteria for init and follow. 
dissplit <- function(file, metafile, init='initial', follow='followUp', split=T, meta = T){
  source('code/read.dist.R')
  dist <- read.dist(file)
  intra <- c()
  # Complete only if split is wanted
  if (meta == T){
    
    for(i in 1:nrow(metafile)){
      intra[i] <- as.numeric(dist[
        as.character(metafile[i, init]),as.character(metafile[i, follow])])
    }
    intra_ade <- intra[metaF$dx=='adenoma']
    intra_canc <- intra[metaF$dx=='cancer']
    
    ###Pull out only interindividual distances and convert to vector
    inter <- dist[as.character(metaF$initial),as.character(metaF$followUp)]
    
    #Remove those that are in the intra and replace with NA
    diag(inter) <- NA
    #Check that it only took out 67 measurements
    length(which(is.na(inter)))
    
    inter <- as.numeric(as.vector(unlist(inter)))
    inter <- inter[!is.na(inter)]
    
    combined <- list(intra_ade, intra_canc, inter)
    names(combined) <- c("intra_ade", "intra_canc", "inter")
    
  }
  
 if (split == T){
    return(combined)
  } else {
    return(dist)
  }
  
  
}


# Create a function that will pick distance values between samples based on user defined
# two vectors and outputs a data frame.
pickDistanceValues <- function(vec1, vec2, distanceTable, metadata, group, withMeta = TRUE){
  
  tempDistance <- c()
  
  if (withMeta == TRUE){
    
    for (i in 1:length(vec1)){
      
      tempDistance <- c(tempDistance, distanceTable[vec1[i], vec2[i]])
      
    }
    
    tempTable <- cbind(as.data.frame(as.numeric(tempDistance)), metadata[, group])
    colnames(tempTable) <- c("distance", group)
    
  } else{
    
    tempTable <- distanceTable[vec1, vec2]
    
  }
  
  
  
  return(tempTable) 
}


# Function to obtain paired wilcoxson tests on probabilities at initial and follow up
# from a given model
# dataTable must have columns named as Dx_Bin and sampleType for it to work
# probability table should have a "No" and "Yes" column for the probability call of the model
getProb_PairedWilcox <- function(dataTable, 
                                 rown = c("lesion", "all_adenoma", "carcinoma_only", "SRN_only"), 
                                 # next two variables are used to set up filtering criteria
                                 lesion_group = c("all", "cancer", "adenoma", "cancer"), 
                                 extra_specifics = c("none", "none", "adv_adenoma", "adenoma")){
  
  # create table to hold wilcoxson paired tests pvalues
  
  wilcox_pvalue_summary <- as.data.frame(matrix(
    nrow = 4, ncol = 1, dimnames = list(
      rows = rown, 
      cols = "Pvalue")))
  
  # Set up variable vector
  lesion_type <- lesion_group
  filter_diagnosis <- extra_specifics
  
  for(i in 1:length(lesion_type)){
    
    wilcox_pvalue_summary[i, "Pvalue"] <- wilcox.test(
      filter(dataTable, sampleType == "initial" & 
               Dx_Bin != paste(lesion_type[i]) & 
               Dx_Bin != paste(filter_diagnosis[i]))[, "Yes"], 
      filter(dataTable, sampleType == "followup" & 
               Dx_Bin != paste(lesion_type[i]) & 
               Dx_Bin != paste(filter_diagnosis[i]))[, "Yes"], 
      paired = TRUE)$p.value
    
  }
  
  # Add Benjamini-Hochberg correction
  wilcox_pvalue_summary <- cbind(
    wilcox_pvalue_summary, 
    BH_correction = p.adjust(wilcox_pvalue_summary$Pvalue, 
                             method = "BH")) 
  
  return(wilcox_pvalue_summary)
}


# Create a label key for the facet wrap function for ggplot2
createTaxaLabeller <- function(taxaTable){
  # collects all entries from the first data list that is not unclassified
  dataList <- apply(taxaTable, 1, function(x) x[x != "unclassified"])
  # creates a vector of the lowest ID taxonomy
  if (class(dataList) == "list"){

    tempCall <- unname(unlist(lapply(dataList, function(x) x[length(x)])))
  
  } else{

    tempCall <- apply(dataList, 2, function(x) x[length(x[x != "Bacteria"])+1])
  }
  
  # assigns names to the vector that are the OTU labels
  names(tempCall) <- rownames(taxaTable)
  # returns that vector
  return(tempCall)
}


# Create function to obtain confusion summary data
# Of importance is the Mcnemar P-value for comparison of actual versus predicted similarity
# needs caret to work properly and needs a meta data table with a column called Disease_Free
# Disease_Free needs to be in a "y"/"n" format
get_confusion_data <- function(dataTable, metaData){
  
  loadLibs(c("caret", "dplyr"))
  #add needed columns for the testing
  dataTable <- mutate(
    dataTable, 
    predict_call = factor(ifelse(Yes > 0.5, "Yes", "No"))) %>% 
    mutate(
      initial_call = factor(c(rep("Yes", length(rownames(metaData))), 
                              ifelse(metaData$Disease_Free == "n", "Yes", "No"))))
  
  #run testing on all the data
  confusion_all <- confusionMatrix(dataTable$predict_call, dataTable$initial_call, 
                                   positive = "Yes")
  
  #run testing based on initial only
  confusion_initial <- confusionMatrix(
    dataTable$predict_call[1:length(rownames(metaData))], 
    dataTable$initial_call[1:length(rownames(metaData))], 
    positive = "Yes")
  
  #run testing based on follow up only
  confusion_follow <- confusionMatrix(
    dataTable$predict_call[(length(rownames(metaData))+1):
                                         length(rownames(dataTable))], 
    dataTable$initial_call[(length(rownames(metaData))+1):
                                         length(rownames(dataTable))], positive = "Yes")
  
  #aggregate all the tests together into a single summary file to be read out
  confusion_summary <- cbind(
    all = c(confusion_all$overall, confusion_all$byClass), 
    initial = c(confusion_initial$overall, confusion_initial$byClass), 
    followup = c(confusion_follow$overall, confusion_follow$byClass))
}


# This function creates a confusion table
# Needs to have a metadata data with the column "Disease_Free"
# This column needs to be in the format of "n"/"y"
make_confusionTable <- function(dataTable, metaData, n = 1, m = 66){
  
  #add needed columns for the testing
  dataTable <- mutate(
    dataTable, 
    predict_call = factor(ifelse(Yes > 0.5, "Yes", "No"))) %>% 
    mutate(
      initial_call = factor(c(rep("Yes", length(rownames(good_metaf))), 
                              ifelse(metaData$Disease_Free == "n", "Yes", "No"))))
  
  temp_list <- confusionMatrix(
    dataTable$predict_call[n:m], 
    dataTable$initial_call[n:m], 
    positive = "Yes")
  
  c_table <- matrix(temp_list$table, nrow = 2, ncol = 2, 
                            dimnames = list(nrow = c("pred_no", "pred_yes"), 
                                            ncol = c("ref_no", "ref_yes")))
}

















