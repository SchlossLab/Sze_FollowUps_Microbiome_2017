### Run the Random Forest on initial versus follow up
### Base analysis without metadata or OTU groupings
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))


# Read in data tables
good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)
shared <- read.delim('data/process/followUps.final.an.unique_list.0.03.subsample.0.03.filter.shared', 
                     header=T, sep='\t')

# Create a vector with the sampleIDs in desired order and required metadata
sampleIDs <- c(good_metaf$initial, good_metaf$followUp)

crc_sampleIDs <- c(good_metaf$initial[good_metaf$Diagnosis != "adenoma"], 
                   good_metaf$followUp[good_metaf$Diagnosis != "adenoma"])

adn_sampleIDs <- c(good_metaf$initial[good_metaf$Diagnosis == "adenoma"], 
                   good_metaf$followUp[good_metaf$Diagnosis == "adenoma"])


# Modify the shared file to only include initial and follow up samples in correct order
rownames(shared) <- as.character(shared$Group)
mod_shared <- shared[as.character(sampleIDs), ]

mod_crc_shared <- shared[as.character(crc_sampleIDs), ]

mod_adn_shared <- shared[as.character(adn_sampleIDs), ]

train_data <- list(All = )

# Create a training data set
train_set <- list(All = mutate(mod_shared, lesion = as.factor(c(rep(0, 67), rep(1, 67)))) %>% 
                    select(lesion, contains("Otu0")), 
                  cancerOnly = mutate(mod_crc_shared, lesion = as.factor(c(rep(0, 26), rep(1, 26)))) %>% 
                    select(lesion, contains("Otu0")), 
                  adnOnly = mutate(mod_adn_shared, lesion = as.factor(c(rep(0, 41), rep(1, 41)))) %>% 
                    select(lesion, contains("Otu0")))

# Set up initial list to store the run data
orig_RF_run <- list(All = c(), cancerOnly = c(), adnOnly = c())
orig_rf_opt <- list(All = c(), cancerOnly = c(), adnOnly = c())
orig_probs <- list(All = c(), cancerOnly = c(), adnOnly = c())
orig_roc <- list(All = c(), cancerOnly = c(), adnOnly = c())

# Run the actual AUCRF for each different group
for(i in 1:length(train_set)){
  
  set.seed(050416)
  orig_RF_run[[i]] <- AUCRF(as.formula(paste(colnames(train_set[[i]])[1], " ~ .", sep = "")), 
                            data=train_set[[i]], pdel=0.05, ntree=500, ranking='MDA')
  orig_rf_opt[[i]] <- orig_RF_run[[i]]$RFopt
  orig_probs[[i]] <- predict(orig_rf_opt[[i]], type = 'prob')[, 2]
  orig_roc[[i]] <- roc(train_set[[i]][, 1] ~ orig_probs[[i]])
  
}

# Remove and load data needed to compare original test run with initial/follow up RF runs
rm(list = setdiff(ls(), "orig_RF_run"))

IF_model <- orig_RF_run
rm(orig_RF_run)

load("exploratory/RFwoFit.RData")
rm(selected_train, sens_specif_table, data, corr_pvalue_ROC_table, orig_probs, impfactor_Data_List, orig_roc, train,  
   selected_probs, selected_RF_run, variableList, modelList, selected_rf_opt, selected_roc, rocNameList, confirmed_vars, orig_rf_opt)

train_model <- orig_RF_run
rm(orig_RF_run)

## Need to think of a way to statistically anayze this.

# Get OTU list from the training model and models created with both adenoma and crc (All), 
# only adenoma (adnOnly), and only crc (crcOnly)
train_otus <- train_model$lesion$Xopt
all_IF_otus <- IF_model$All$Xopt
crc_IF_otus <- IF_model$cancerOnly$Xopt
adn_IF_otus <- IF_model$adnOnly$Xopt


# Find number of overlapping otus and keep those that overlap.
all_ID_otus <- train_otus[train_otus %in% all_IF_otus]
crc_ID_otus <- train_otus[train_otus %in% crc_IF_otus]
adn_ID_otus <- train_otus[train_otus %in% adn_IF_otus]

### Create a data table to graph these specific OTUs

# Read in and create a useable taxonomy table
tax <- read.delim('data/process/followUps.final.an.unique_list.0.03.cons.taxonomy', 
                  sep='\t', header=T, row.names=1)

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)
# This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components


# create needed labels for Boruta picked important variables for each model

all_selected_taxa <- tax_df[all_ID_otus, ]
all_selected_labs <- createTaxaLabeller(all_selected_taxa)


# Create needed initial and follow up data tables
initial <- inner_join(good_metaf, shared, by = c("initial" = "Group")) %>%
  select(lesion, one_of(all_ID_otus))

followups <- inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
  select(lesionf, one_of(all_ID_otus)) %>% rename(lesion = lesionf)

#creating the data tables to be used
select_otu_data <- melt(
  rbind(mutate(initial, sampleType = "initial", EDRN = good_metaf$EDRN, Diagnosis = good_metaf$Diagnosis), 
        mutate(followups, sampleType = "follow_up", EDRN = good_metaf$EDRN, Diagnosis = good_metaf$Diagnosis)), 
  id=c(paste("lesion"), "EDRN", "Diagnosis", "sampleType")) %>% 
  mutate(Disease_Free = rep(good_metaf$Disease_Free, length(rownames(.))/67))

filter(select_otu_data, Diagnosis != "adenoma") %>% 
ggplot(aes(factor(sampleType, levels = c("initial", "follow_up")), 
                            log10(value+1.1), group = factor(EDRN))) + 
  geom_line(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
  geom_point(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
  facet_wrap(~variable, labeller = as_labeller(all_selected_labs)) + coord_cartesian(ylim = c(0, 4)) + 
  theme_bw() + ylab("Log Total Sequences") + xlab("") + ggtitle("Cancer") + 
  scale_colour_manual(name = "Cancer Free", labels = c("No", "Yes", "Unknown"), values = c("darkred", "orange", "red")) + 
  theme(legend.position=c(0.85,0.05), plot.margin = unit(c(1, 1, 1, 1), "lines"), 
        plot.title = element_text(size=20, face="bold"), legend.title = element_text(face="bold"))

###Re-Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson"))

# Run paired wilcoxson statistical tests

ALL_adjusted_pvalues_abund <- get_abund_pvalues(initial, followups)

crc_adjusted_pvalues_abund <- get_abund_pvalues(
  inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
    filter(Diagnosis != "adenoma") %>% select(lesion, one_of(crc_ID_otus)), 
  inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
    filter(Diagnosis != "adenoma") %>% 
    select(lesionf, one_of(crc_ID_otus)) %>% rename(lesion = lesionf))

adn_adjusted_pvalues_abund <- get_abund_pvalues(
  inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
    filter(Diagnosis == "adenoma") %>% select(lesion, one_of(adn_ID_otus)), 
  inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
    filter(Diagnosis == "adenoma") %>% 
    select(lesionf, one_of(adn_ID_otus)) %>% rename(lesion = lesionf))


