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
      dataTable[, k] <- p.adjust(dataTable[, k], method = "bonferroni", n = totalComparisons)
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
    
    tempCall <- apply(dataList, 2, function(x) x[length(x[x != "unclassified"])])
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













create.table <- function(initial_data, follow_data, otus, variable_list){
  otuData <- unname(unlist(c(
    initial_data[, otus], follow_data[, otus])))
  
  temp_file <- c()
  for(i in 1:length(otus)){
    temp_file <- c(temp_file, rep(otus[i], length(rownames(initial_c_shared))))
    
    if(i == length(otus)){
      OTUcat <- rep(temp_file, 
                   length(otuData)/length(temp_file))
    }
  }
  
}


test_operation <- function(x){
  AUCRF(lesion~., data=x, pdel=0.05, ntree=500, ranking='MDA')
}








