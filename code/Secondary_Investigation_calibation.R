### Calibration plots of the full models
### Investigate model fit on follow up data
## Marc Sze

### Hold in storage in case someone is really interested in the model performances


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", 
           "wesanderson", "caret", "randomForest"))

# Run full model on test data calibration test
full_test <- read.csv("results/tables/lesion_test_data_probs_summary.csv", header = T, 
                      row.names = 1, stringsAsFactors = F)

full_test2 <- as.data.frame(cbind(obs = ifelse(full_test$lesion == "Yes", 1, 0), probs = full_test$No))
full_test2$obs <- factor(full_test2$obs)

full_cali <- calibration(obs ~ probs, data = if_test2, cuts = 10)

plot(calibration(obs ~ probs, data = full_test2, cuts = 10))


# Run reduced full model on test data calibration test
red_test <- read.csv("results/tables/reduced_lesion_test_data_probs_summary.csv", 
                     header = T, stringsAsFactors = F)

red_test2 <- as.data.frame(cbind(obs = ifelse(red_test$lesion == "Yes", 1, 0), probs = red_test$No))
red_test2$obs <- factor(red_test2$obs)

red_cali <- calibration(obs ~ probs, data = red_test2, cuts = 10)

plot(calibration(obs ~ probs, data = red_test2, cuts = 10))


# Run full model on follow up calibration test
test <- read.csv("results/tables/follow_up_probability_summary.csv", header = T, stringsAsFactors = F)

test2 <- as.data.frame(cbind(obs = c(rep(1, 66), rep(1, 1), rep(0, 65)), probs = test$No))
test2$obs <- factor(test2$obs)

plot(calibration(obs ~ probs, data = test2, cuts = 10))

# Run IF model on follow up calibration test
if_test <- read.csv("results/tables/IF_follow_up_probability_summary.csv", header = T, stringsAsFactors = F)

if_test2 <- as.data.frame(cbind(obs = c(rep(1, 66), rep(1, 1), rep(0, 65)), probs = if_test$No))
if_test2$obs <- factor(if_test2$obs)

if_cali <- calibration(obs ~ probs, data = if_test2, cuts = 2)

plot(calibration(obs ~ probs, data = if_test2, cuts = 10))

















