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
dissplit <- function(file, metafile, init='initial', follow='followUp'){
  source('code/read.dist.R')
  dist <- read.dist(file)
  intra <- c()
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
  
  return(combined)
}


# This function obtains the Wilcoxson Rank Sum and generates the bonferroni corrected P-values based on supplied groups
getCorrectedPvalue <- function(testTable, variableList, group){
  # Load needed library and initialize table
  loadLibs("dplyr")
  finalTable <- c()
  # Iterate through apply function to analyze p-value
  for (i in 1:length(group)){
    # Runs that actual wilcoxson test across every column
    tempData <- apply(select(testTable, one_of(variableList)), 2, function(x) 
      wilcox.test(x~testTable[, group[i]])$p.value) %>% p.adjust(method = "bonferroni")
    # Aggregate the table together
    finalTable <- cbind(finalTable, tempData)
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
  tempCall <- unname(unlist(lapply(dataList, function(x) x[length(x)])))
  # assigns names to the vector that are the OTU labels
  names(tempCall) <- rownames(taxaTable)
  # returns that vector
  return(tempCall)
}



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








