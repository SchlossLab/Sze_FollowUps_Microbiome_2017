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





# This function obtains the Wilcoxson Rank Sum and generates the bonferroni corrected P-values based on supplied groups
getCorrectedPvalue <- function(testTable, variableList, group, wilcox = T, kruskal = F){
  # Load needed library and initialize table
  loadLibs("dplyr")
  finalTable <- c()
  # Iterate through apply function to analyze p-value
  for (i in 1:length(group)){
    
    if (wilcox == T){
      
      # Runs that actual wilcoxson test across every column
      tempData <- apply(select(testTable, one_of(variableList)), 2, function(x) 
        wilcox.test(x~testTable[, group[i]])$p.value) %>% p.adjust(method = "bonferroni")
      # Aggregate the table together
      finalTable <- cbind(finalTable, tempData)
    }
    
    if (kruskal == T){
      
      # Runs that actual wilcoxson test across every column
      tempData <- apply(select(testTable, one_of(variableList)), 2, function(x) 
        kruskal.test(x~testTable[, group[i]])$p.value) %>% p.adjust(method = "bonferroni")
      # Aggregate the table together
      finalTable <- cbind(finalTable, tempData)
    }
    
  }
  # Change column names and return final table as a data frame
  colnames(finalTable) <- group
  return(as.data.frame(finalTable))
}


# This function creates a table from different ROC lists and creates a ggplot 
# useable table of sensitivities and specificities
makeSensSpecTable <- function(rocNameList, variableList, modelList){
  # Load needed library if not already installed and loaded and initialize empty table
  loadLibs("dplyr")
  sens_specif_table <- c()
  # loop to iterate through each of the listed ROC lists
  for (i in 1:length(modelList)){
    # Gets senstivities and specificities and converts them to a data frame
    tempfile <- as.data.frame.list(rocNameList[[i]][variableList])
    # create a column that specifies model used
    tempfile2 <- mutate(tempfile, model = rep(modelList[i], length(rownames(tempfile))))
    # Add this to existing data frame
    sens_specif_table <- rbind(sens_specif_table, tempfile2)
  }
  # Align factor levels with order in modeList
  sens_specif_table$model <- factor(sens_specif_table$model, levels = modelList)
  # Output final table
  return(sens_specif_table)
}


# This function compares all possible ROC curves generated against each other and can 
# initialize a multiple comparison correction if necessary.
getROCPvalue <- function(rocNameList, modelList, totalModels, multi = F){
  # Load needed library and set up initial conditions
  loadLibs("pROC")
  dataTable <- c()
  x <- 1
  # Set up initial loop to run through the different models
  for (i in 1:totalModels){
    # Make sure that it only compares groups that have not had any comparisons completed yet
    if (i != totalModels){
      # Create an empty data table
      tempVector <- c()
      # Initialize second loop for the comparison group
      for (j in (i+1):totalModels){
        # Run the actual DeLong's test for correlated ROC curves and grab p-value
        tempVector <- c(tempVector, unname(unlist(roc.test(rocNameList[[i]], rocNameList[[j]])['p.value'])))
      }
      # Initialize statement for adding NA's if comparing against self
      if (i != 1){
        # Runs the NA addition and adds to counter
        tempVector <- c(rep(NA, totalModels - (totalModels - x)), tempVector)
        x <- x + 1
      }
      # Combine newly generated vector to existing stored data table
      dataTable <- cbind(dataTable, tempVector)
    }
  }
  # Add row names and column names for easier interpretation
  rownames(dataTable) <- modelList[2:totalModels]
  colnames(dataTable) <- modelList[1:(totalModels - 1)]
  # If bonferroni correction is wanted then it will enter here before returning the data table
  if (multi == T){
    # Get the total number of comparisons that were made
    totalComparisons <- unname(table(is.na(dataTable))["FALSE"])
    # Create a for loop bounded by the total comparisons made
    for (k in 1:length(colnames(dataTable))){
      # Run the p adjust base function over each column
      dataTable[, k] <- p.adjust(dataTable[, k], method = "BH", n = totalComparisons)
    }
  }
  
  # Output final data table with p-values for each comparison
  return(as.data.frame(dataTable))
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


# Get similarities of Important variables between models for a VennDiagram
# Dereplicated since I decided not to got this route and 
# Did not see the utility in having this information.
compare4ModelVariables <- function(variableList){
  
  loadLibs("VennDiagram")
  
  tempList <- list(n12 = c(), n13 = c(), n14 = c(), n23 = c(), n24 = c(), n34 = c(), 
                   n123 = c(), n124 = c(), n134 = c(), n234 = c(), n1234 = c())
  
  tempSimOTUList <- list(n12 = c(), n13 = c(), n14 = c(), n23 = c(), n24 = c(), n34 = c())
  
  x = 1
  
  for (i in 1:length(variableList)){
    
    tempValues <- c()
    if (i != length(variableList)){
      
      for (j in (i+1):length(variableList)){
        
        tempList[[x]] <- length(calculate.overlap(variableList[c(i,j)])[[3]])
        tempSimOTUList[[x]] <- calculate.overlap(variableList[c(i,j)])[[3]]
    
        x <- x + 1
      }
      
    }

  }
  tempList[['n123']] <- length(calculate.overlap(list(variableList[[3]], tempSimOTUList[['n12']])))
  tempList[['n124']] <- length(calculate.overlap(list(variableList[[4]], tempSimOTUList[[1]])))
  tempList[['n134']] <- length(calculate.overlap(list(variableList[[4]], tempSimOTUList[[2]])))
  tempList[['n234']] <- length(calculate.overlap(list(variableList[[4]], tempSimOTUList[[4]])))
  tempList[['n123']] <- length(calculate.overlap(list(tempSimOTUList[['n13']], tempSimOTUList[['n24']])))
  
  
  return(tempList)
}



# Function that creates a 2 x 2 table looking at the number of positives and negatives predicted
makePosNegTable <- function(dataTable, variable = "positive", cutoffTable, rfDone = "threeway_rf_opt", model = "threeGroups", 
                            cutdata_selection = "All", cutoffCall = "cutoff", cutModel_call = "model", datasetCall = "dataset", 
                            filterDiagnosis = "adenoma", SampleTime = "initial", x = 2, y = 2){
  # Load required library
  loadLibs("dplyr")
  # create the cutoff value to be used
  test <- cutoffTable[, cutoffCall][withFIT_cutoffs[, cutModel_call] == paste(model) & 
                                          withFIT_cutoffs[, datasetCall] == paste(cutdata_selection)]
  # create an empty matrix to enter results into and name the rows and columns
  test2 <- as.data.frame(matrix(nrow = x, ncol = y))
  rownames(test2) <- c(paste(filterDiagnosis, "Group"), paste("Not", filterDiagnosis, "Group"))
  colnames(test2) <- c("PredictPos", "PredictNeg")
  
  # Look at first group predicted positives
  test2[paste(filterDiagnosis, "Group"), "PredictPos"] <- length(rownames(filter(
    dataTable[dataTable[, variable] > test, ]  %>% filter(L1 == rfDone[1])  %>% filter(diagnosis == filterDiagnosis[1]), 
    time_point == SampleTime[1])))
  
  # Look at first group predicted negatives
  test2[paste(filterDiagnosis, "Group"), "PredictNeg"] <- length(rownames(filter(
    dataTable[dataTable[, variable] < test, ]  %>% filter(L1 == rfDone[1])  %>% filter(diagnosis == filterDiagnosis[1]), 
    time_point == SampleTime[1])))
  
  # Look at second group predicted positives
  test2[paste("Not", filterDiagnosis, "Group"), "PredictPos"] <- length(rownames(filter(
    dataTable[dataTable[, variable] > test, ]  %>% filter(L1 == rfDone[1])  %>% filter(diagnosis != filterDiagnosis[1]), 
    time_point == SampleTime[1])))
  
  # Look at second groups predicted negatives
  test2[paste("Not", filterDiagnosis, "Group"), "PredictNeg"] <- length(rownames(filter(
    dataTable[dataTable[, variable] < test, ]  %>% filter(L1 == rfDone[1])  %>% filter(diagnosis != filterDiagnosis[1]), 
    time_point == SampleTime[1])))
  
  return(test2)
  
}


# Test function for foreach loop for modularizing AUCRF
testdata <- function(dataList, group){
  
  loadLibs(c("dplyr", "AUCRF", "pROC"))
  
  #DataTable <- dataList[[group]] 
  
  assign("dataList", dataList, envir = globalenv())
  assign("group", group, envir = globalenv())
  
  set.seed(050416)
  assign("AUCdata", AUCRF(as.formula(paste(group, " ~ .", sep = "")), 
                          dataList[[group]], pdel=0.05, ntree=500, ranking='MDA'), envir = globalenv())
  
  test_rf_opt <- AUCdata$RFopt
  
  assign("model_train_probs", predict(test_rf_opt, type='prob')[, 2], envir = globalenv())
  
  model_train_roc <- roc(as.formula(paste("dataList", "[[", "group", "]]", "$", group, "~ model_train_probs", sep = "")))
  
  set.seed(050416)
  modelAUCRF_wCV <- AUCRFcv(AUCdata, 10, 20) # if works need to set this back to 10, 20
  
  modelAUCRF_wCV_selected <- as.data.frame(modelAUCRF_wCV$Psel[modelAUCRF_wCV$Psel > 0.5])
  colnames(modelAUCRF_wCV_selected) <- "Probs_present"
  
  write.csv(modelAUCRF_wCV_selected, 
            paste("data/process/old/",group,"_AUCRFwCV.csv", sep = ""))
  #if works need to set back to data/process/old
  
  
  
  finalList <- list(model = AUCdata, modelRFopt = test_rf_opt, probs = model_train_probs, train_roc = model_train_roc, 
                    wCV = modelAUCRF_wCV, selectedOTUs = modelAUCRF_wCV_selected)
  

  rm(group, dataList, model_train_probs, AUCdata, envir = globalenv())
  
  return(finalList)
  
  
}


# A function to convert into a categorical based on a defined cutoff
makeCat <- function(dataTable, cutoff, ignorecols = 2){
  
  tempData <- as.data.frame(matrix(nrow = length(rownames(dataTable)), ncol = length(colnames(dataTable))))
  
  for(i in (ignorecols+1):length(colnames(dataTable))){
    tempVector <- rep(NA, length(dataTable[, i]))
    tempVector[dataTable[, i] <= cutoff] <- 0
    tempVector[dataTable[, i] > cutoff] <- 1
    tempData[, i] <- tempVector
  }
  
  for(j in 1:ignorecols){
    
    tempData[, j] <- dataTable[, j] 
  }
  
  colnames(tempData) <- colnames(dataTable)
  
  return(tempData)
}


# Function to count the number of appearances that are greater than 0
getCount <- function(vectorOfinterest){
  
  count = 0
  for (i in 1:length(vectorOfinterest)){
    if(vectorOfinterest[i] > 0){
      
      count = count + 1

    } else {
      
      count = count
    }
    
  }
  return(count)
}



# Function to create a count table
# Create so that initial and follow ups are created from one table and not overwriting themselves
# Add place for disease group to be in a column

createCountTable <- function(data, cutoff, L1_analysis, time_point_analyzed, diagnosis_used){
  
  testTable <- as.data.frame(matrix(nrow = length(diagnosis_used)*length(time_point_analyzed), ncol = 4))
  colnames(testTable) <- c("timePoint", "Predicted", "GroupTotal", "DiseaseGroup")
  
  x <- 1
  
  for(i in 1:length(time_point_analyzed)){
 
    for(j in 1:length(diagnosis_used)){
      
      testTable$Predicted[x] <- 
        length(rownames(filter(withFit_data, L1 == paste(L1_analysis), time_point == paste(time_point_analyzed[i]), 
                               detailed_diagnosis == paste(diagnosis_used[j]), positive > cutoff)))
      
      testTable$timePoint[x] <- time_point_analyzed[i]
      
      testTable$GroupTotal[x] <- 
        length(rownames(filter(withFit_data, L1 == paste(L1_analysis), time_point == paste(time_point_analyzed[i]), 
                               detailed_diagnosis == paste(diagnosis_used[j])))) 
      
      testTable$DiseaseGroup[x] <- diagnosis_used[j]
      
      x <- x + 1
     
    }
  }
  
  return(testTable)
}




# obtain wilcoxson p-values from initial and follow up data tables and correct for multiple comparisons
get_abund_pvalues <- function(initialData, followupData, multi = T, correction = "bonferroni", 
                              pairedTest = TRUE, start = 2){
  
  pvalues <- c()
  
  for(i in start:length(colnames(initialData))){
    
    pvalues <- c(pvalues, 
                 wilcox.test(initialData[, i], followupData[, i], paired = pairedTest)$p.value)
  }
  
  if(multi == T){
    
    padjust_vals <- p.adjust(pvalues, method = paste(correction))
    
    temptable <- as.data.frame(cbind(colnames(initialData)[start:length(colnames(initialData))], pvalues, padjust_vals))
    colnames(temptable) <- c("otus", "pvalues", "adj_pvalues")
    temptable$pvalues <- as.numeric(as.character(temptable$pvalues))
    temptable$adj_pvalues <- as.numeric(as.character(temptable$adj_pvalues))
    temptable$otus <- as.character(temptable$otus)
    
    
  }else{
    
    temptable <- as.data.frame(cbind(colnames(initialData)[start:length(colnames(initialData))], pvalues))
    colnames(temptable) <- c("otus", "pvalues")
    as.data.frame(temptable)
    temptable$pvalues <- as.numeric(as.character(temptable$pvalues))
    temptable$otus <- as.character(temptable$otus)
  }
  
 return(temptable)
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

















