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
get_alpha_pvalues <- function(data_table, table_names = c("sobs", "shannon", "evenness"), 
                              numComp = length(table_names)){
  
  #Set up output table
  alpha_pvalue_table <- as.data.frame(matrix(ncol = 1, nrow = numComp, 
                                             dimnames = list(rows = table_names, 
                                                             cols = "pvalue")))
  #Get p-values of paired test
  alpha_pvalue_table[, "pvalue"] <- apply(data_table[, c("sobs", "shannon", "shannoneven")], 2, 
                                          function(x){
                                            wilcox.test(x[data_table$sampleType == "initial"], 
                                                        x[data_table$sampleType == "followups"], 
                                                        paired = TRUE)$p.value})
  
  return(alpha_pvalue_table)
}


# This loads in a mothur triangle matrix and converts it to a square matrix
read.dist <- function(file, input='lt', make.square=T, diag=0){
  if(input=='lt'){
    stuff <- scan(file, what='character') #gets all of the elements of the file
    n <- as.numeric(stuff[1]) # n = number of groups/samples in file
    stuff <- stuff[-1] # removes number of groups from list of stuff
    m <- data.frame(matrix(NA, nrow=n, ncol=n) ) #makes empty matrix based on number of groups
    diag(m) <- diag #fills in diagonal with specified value
    
    c <- 1 # c keeps track of postion in stuff vector
    for(i in 1:n){
      group <- stuff[c] #get group name
      colnames(m)[i] <- group
      rownames(m)[i] <- group
      
      if(i > 1){
        m[i,1:(i-1)] <- stuff[(c+1):(c+i-1)] #fills in matrix with values from stuff
      }
      c<-c+i #this because math
    }
    if(make.square){
      m[upper.tri(m)] <- t(m)[upper.tri(m)] #fills in upper triangle
    }
  }
  
  if(input=='square'){ #reads in square matrix
    m<-read.table(file, skip=1, row.names=1)
    colnames(m) <- rownames(m)
  }
  return(m)
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
                                 not_group = c("all", "cancer", "adenoma", "cancer"), 
                                 extra_specifics = c("none", "none", "adv_adenoma", "adenoma")){
  
  # create table to hold wilcoxson paired tests pvalues
  
  wilcox_pvalue_summary <- as.data.frame(matrix(
    nrow = length(rown), ncol = 1, dimnames = list(
      rows = rown, 
      cols = "Pvalue")))
  
  # Set up variable vector
  lesion_type <- not_group
  filter_diagnosis <- extra_specifics
  
  for(i in 1:length(rown)){
    
    wilcox_pvalue_summary[i, "Pvalue"] <- wilcox.test(
      filter(dataTable, sampleType == "initial" & 
               Dx_Bin != paste(lesion_type[i]) & 
               Dx_Bin == paste(filter_diagnosis[i]))[, "Yes"], 
      filter(dataTable, sampleType == "followup" & 
               Dx_Bin != paste(lesion_type[i]) & 
               Dx_Bin == paste(filter_diagnosis[i]))[, "Yes"], 
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

# Creates names for the get_tax_level_shared function
get_tax_substring <- function(tax, tax_level){
  substring <- unlist(strsplit(tax, ";"))[tax_level]
  paste(substring, collapse='.')
}
get_tax_name <- function(tax_file, tax_level){
  
  
  tax_data <- read.table(file=tax_file, header=T, stringsAsFactors=F)
  taxonomy <- tax_data$Taxonomy
  taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
  taxonomy <- gsub('"', '', taxonomy)
  
  tax_substring <- sapply(taxonomy, get_tax_substring, tax_level)
  
  names(tax_substring) <- tax_data$OTU
  
  tax_substring
}

# Get the total number based on shared file and tax file.
get_tax_level_shared <- function(shared_file, tax_file, tax_level){
  
  shared_otus <- read.delim(file=shared_file, header=T, stringsAsFactors=F, row.names=2)[,-c(1,2)]
  is_present <- apply(shared_otus, 2, sum) > 0
  shared <- shared_otus[,is_present]
  
  taxonomy <- get_tax_name(tax_file, tax_level)
  taxonomy <- taxonomy[colnames(shared)]
  unique_taxa <- levels(as.factor(taxonomy))
  
  shared_tax_level <- NULL
  
  for(ut in unique_taxa){
    otus <- names(taxonomy[taxonomy %in% ut])
    sub_shared <- shared_otus[,colnames(shared_otus) %in% otus]
    
    if(is.null(dim(sub_shared))){
      shared_tax_level <- cbind(shared_tax_level, sub_shared)
    } else {
      tax_level_count <- apply(sub_shared, 1, sum)
      shared_tax_level <- cbind(shared_tax_level, tax_level_count)
    }
  }
  colnames(shared_tax_level) <- unique_taxa
  rownames(shared_tax_level) <- rownames(shared)
  return(shared_tax_level)
}




# Create function to obtain confusion summary data
# Of importance is the Mcnemar P-value for comparison of actual versus predicted similarity
# needs caret to work properly and needs a meta data table with a column called Disease_Free
# Disease_Free needs to be in a "y"/"n" format
get_confusion_data <- function(dataTable, metaData, 
                               to_filter1 = "cancer", to_filter2 = "cancer", columns = "Dx_Bin", DF_data = "No"){
  
  loadLibs(c("caret", "dplyr"))
  
  tempMetaData <- metaData[metaData[, columns] != to_filter1 & metaData[, columns] != to_filter2, ]
  dataTable <- dataTable[dataTable[, columns] != to_filter1 & dataTable[, columns] != to_filter2, ]
  n <- length(rownames(tempMetaData))
  
 
  #add needed columns for the testing
  if(DF_data == "No"){
    
    dataTable <- mutate(
      dataTable, 
      predict_call = factor(ifelse(Yes > 0.5, "Yes", "No"))) %>% 
      mutate(
        initial_call = factor(c(rep("Yes", n), rep("No", n))))
  } else{
    
    dataTable <- mutate(
      dataTable, 
      predict_call = factor(ifelse(Yes > 0.5, "Yes", "No"))) %>% 
      mutate(
        initial_call = factor(c(rep("Yes", n), 
                                ifelse(tempMetaData$Disease_Free == "n", "Yes", "No"))))
  }
  
  
  #run testing on all the data
  confusion_all <- confusionMatrix(dataTable$predict_call, dataTable$initial_call, 
                                   positive = "Yes")
  
  #run testing based on initial only
  confusion_initial <- confusionMatrix(
    dataTable$predict_call[1:length(rownames(tempMetaData))], 
    dataTable$initial_call[1:length(rownames(tempMetaData))], 
    positive = "Yes")
  
  #run testing based on follow up only
  confusion_follow <- confusionMatrix(
    dataTable$predict_call[(length(rownames(tempMetaData))+1):
                                         length(rownames(dataTable))], 
    dataTable$initial_call[(length(rownames(tempMetaData))+1):
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
make_confusionTable <- function(dataTable, metaData, 
                                column = "Dx_Bin", to_filter1 = "cancer", 
                                to_filter2 = "cancer", DF_data = "No", sampling = "sampleType", samples = "initial"){
  
  
  tempMetaData <- metaData[metaData[, "Dx_Bin"] != to_filter1 & metaData[, "Dx_Bin"] != to_filter2, ]
  dataTable <- dataTable[dataTable[, column] != to_filter1 & dataTable[, column] != to_filter2 & 
                           dataTable[, sampling] == samples, ]
  n <- length(rownames(tempMetaData))
  
  #add needed columns for the testing
  if(DF_data == "No"){
    
    dataTable <- mutate(
      dataTable, 
      predict_call = factor(ifelse(Yes > 0.5, "Yes", "No"))) %>% 
      mutate(
        initial_call = factor(ifelse(dataTable[, sampling] == "initial", "Yes", "No"), levels = c("No", "Yes")))
  } else{
    
    dataTable <- mutate(
      dataTable, 
      predict_call = factor(ifelse(Yes > 0.5, "Yes", "No"))) %>% 
      mutate(
        initial_call = factor(c(rep("Yes", n), 
                                ifelse(tempMetaData$Disease_Free == "n", "Yes", "No"))))
  }
  
  
   temp_list <- confusionMatrix(
    dataTable$predict_call, 
    dataTable$initial_call, 
    positive = "Yes")
  
  c_table <- matrix(temp_list$table, nrow = 2, ncol = 2, 
                            dimnames = list(nrow = c("pred_no", "pred_yes"), 
                                            ncol = c("ref_no", "ref_yes")))
}

# Function to track co-occurances of specific OTUs
# Sums all times in a sample each OTU is found together 
# Normalizes by the OTU pair that has the least number of positive samples
get_co_occurance <- function(comparisonTable, variables){
  
  dataTable <- c()
  
  for(i in 1:length(variables)){
    
    coldata <- c()
    maximums <- c()
    
    for(j in 1:length(variables)){
      tempdata <- c()
      
      
      for(k in 1:length(rownames(comparisonTable))){
        y = 0
        
        if(comparisonTable[k, j] > 0 & comparisonTable[k, i] > 0){
          
          y = y + 1
          
        }
        tempdata <- c(tempdata, y)
        
      }
      
      maximums <- c(maximums, ifelse(
        length(comparisonTable[, i][comparisonTable[, i] > 0]) > length(comparisonTable[, j][comparisonTable[, j] > 0]), 
        length(comparisonTable[, j][comparisonTable[, j] > 0]), 
        length(comparisonTable[, i][comparisonTable[, i] > 0])))
      
      coldata <- rbind(coldata, sum(tempdata))
    }
    
    
    dataTable <- cbind(dataTable, as.numeric(format(coldata/maximums[i], digits = 3)))
    
  }
  
  dataTable[upper.tri(dataTable)] <- t(dataTable)[upper.tri(dataTable)]
  
  colnames(dataTable) <- c(1:length(variables))
  rownames(dataTable) <- c(1:length(variables))
  
  return(dataTable)
}


# Function to get list of co-occurances above a cutoff
get_list_coocur <- function(dataTable, cutoff = 0.5){
  
  x = 1
  test <- list()
  
  for(i in 1:length(rownames(dataTable))){
    
    OTUs <- c()
    
    if(i != length(rownames(dataTable))){
      
      for(j in (1+x):length(rownames(dataTable))){
        
        if(dataTable[j, i] > cutoff){
          
          OTUs <- c(OTUs, as.numeric(rownames(dataTable)[j]))
        }
        
      }
      
      test[[colnames(dataTable)[i]]] <- OTUs
      x = x + 1
    }
    
  }
  return(test)
}


# Generate a total count for the total number of connections
get_connection_totals <- function(dataList){
  
  x = 0
  for(i in 1:length(dataList)){
    
    for(j in 1:length(dataList[[i]])){
      
      x = x + 1
    }
  }
  return(x)
}


















