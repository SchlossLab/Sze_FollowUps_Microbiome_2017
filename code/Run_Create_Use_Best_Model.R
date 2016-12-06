### Create actual model based on training and testing
### Use the full data set now to create model
### Workflow based on Jenna Wiens suggested workflow
## Marc Sze

#Load needed libraries
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "caret","scales", "wesanderson", "randomForest"))


# Read in necessary data frames

test_data <- read.csv("results/tables/full_test_data.csv", header = TRUE)

# Create Random FOrest model




