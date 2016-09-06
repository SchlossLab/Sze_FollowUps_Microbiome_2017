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








