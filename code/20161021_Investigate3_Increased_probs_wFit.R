### Investigate Reason what causes 3 samples in with Fit model postive probability to increase
## Marc Sze


# Load in needed libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("pROC","randomForest","AUCRF", "Boruta", "dplyr", "tidyr", "ggplot2", "reshape2", 
           "gridExtra", "scales", "wesanderson", "vegan", "knitr"))


# Read in needed data
withFit_cutoffs <- read.csv("results/tables/withFIT.cutoffs.csv", header = T, stringsAsFactors = F)
woFit_cutoffs <- read.csv("results/tables/noFIT.cutoffs.csv", header = T, stringsAsFactors = F)

withFit_data <- read.csv("results/tables/withFIT.models.datatable.csv", header = T, stringsAsFactors = F)
woFit_data <- read.csv("results/tables/noFIT.models.datatable.csv", header = T, stringsAsFactors = F)

good_metaf <- read.csv("results/tables/mod_metadata/good_metaf_final.csv", stringsAsFactors = F, header = T)



# Generate a vector which contains all samples EDRN that increased in with fit model
higher_prob_samples <- withFit_data[withFit_data$positive[withFit_data$time_point == "followup"] > 
                       withFit_data$positive[withFit_data$time_point == "initial"] & withFit_data$model == "lesion" & 
    withFit_data$dataset == "All" & withFit_data$diagnosis != "adenoma" & 
      withFit_data$time_point == "initial", "EDRN"]

# Subset data for the wofit model to only contain these three samples
woFit_select3 <- filter(woFit_data, EDRN %in% higher_prob_samples, model == "lesion", dataset == "All")
  # don't seem to decrease by much or increase by much (pretty much stay the same with respect to probabilities)

# Subset metadata to observe what that looks like
select3_metadata <- filter(good_metaf, EDRN %in% higher_prob_samples)
  # 2/3 FITs increased in value (However, not the one with Cancer still present (EDRN 13710))
  # This inidivual increases in Parvimonas in the non-fit model but this is not in the with fit model.

# Read in shared and taxonomy data to look at what isn't changing in these three samples

shared <- read.delim('data/process/followUps.final.an.unique_list.0.03.subsample.0.03.filter.shared', 
                     header=T, sep='\t')

tax <- read.delim('data/process/followUps.final.an.unique_list.0.03.cons.taxonomy', 
                  sep='\t', header=T, row.names=1)

wfit_pvalue_summary <- read.csv("results/tables/lesion_OTU_WF_Pvalue_summary.csv", header = T, stringsAsFactors = F)
wofit_pvalue_summary <- read.csv("results/tables/lesion_OTU_WoF_Pvalue_summary.csv", header = T, stringsAsFactors = F)

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
rm(tax)
# This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components

# create needed labels for Boruta picked important variables for each model

lesion_selected_taxa <- tax_df[wfit_pvalue_summary[, 'otus'], ]
lesion_selected_labs <- createTaxaLabeller(lesion_selected_taxa)
lesion_selected_labs <- c(fit_result = "fit", lesion_selected_labs)



# Create needed initial and follow up data tables
initial <- inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
  filter(EDRN %in% higher_prob_samples) %>% 
  select(lesion, fit_result, one_of(wfit_pvalue_summary[, 'otus']))

followups <- inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
  filter(EDRN %in% higher_prob_samples) %>% 
  select(lesionf, fit_followUp, one_of(wfit_pvalue_summary[, 'otus'])) %>% 
  rename(lesion = lesionf) %>% rename(fit_result = fit_followUp)

#creating the data tables to be used
select_otu_data <- melt(
  rbind(mutate(initial, sampleType = "initial", 
                 EDRN = filter(good_metaf, EDRN %in% higher_prob_samples)[, "EDRN"]), 
        mutate(followups, sampleType = "follow_up", 
               EDRN = filter(good_metaf, EDRN %in% higher_prob_samples)[, "EDRN"])), 
  id=c(paste("lesion"), "sampleType", "EDRN")) %>% 
  mutate(Disease_Free = 
           rep(c("n", "y", "y"), length(rownames(.))/3))

ggplot(select_otu_data, aes(factor(sampleType, levels = c("initial", "follow_up")), 
             log10(value+1.1), group = factor(EDRN))) + 
  geom_line(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
  geom_point(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
  facet_wrap(~variable, labeller = as_labeller(lesion_selected_labs)) + coord_cartesian(ylim = c(0, 4)) + 
  theme_bw() + ylab("Log Total Sequences") + xlab("") + ggtitle("Cancer") + 
  scale_colour_manual(name = "Cancer Free", labels = c("No", "Yes", "Unknown"), values = c("darkred", "orange", "red")) + 
  theme(legend.position=c(0.85,0.05), plot.margin = unit(c(1, 1, 1, 1), "lines"), 
        plot.title = element_text(size=20, face="bold"), legend.title = element_text(face="bold"))

    # Based on graph can say that common to all three samples are that they are losing a Lachnospiraceae OTU
    # and gaining two different Bacteroides OTUs

## Try the same thing but wofit
# create needed labels for Boruta picked important variables for each model

lesion_selected_taxa <- tax_df[wofit_pvalue_summary[, 'otus'], ]
lesion_selected_labs <- createTaxaLabeller(lesion_selected_taxa)

# Create needed initial and follow up data tables
initial <- inner_join(good_metaf, shared, by = c("initial" = "Group")) %>% 
  filter(EDRN %in% higher_prob_samples) %>% 
  select(lesion, one_of(wofit_pvalue_summary[, 'otus']))

followups <- inner_join(good_metaf, shared, by = c("followUp" = "Group")) %>% 
  filter(EDRN %in% higher_prob_samples) %>% 
  select(lesionf, one_of(wofit_pvalue_summary[, 'otus'])) %>% 
  rename(lesion = lesionf)

#creating the data tables to be used
select_otu_data <- melt(
  rbind(mutate(initial, sampleType = "initial", 
               EDRN = filter(good_metaf, EDRN %in% higher_prob_samples)[, "EDRN"]), 
        mutate(followups, sampleType = "follow_up", 
               EDRN = filter(good_metaf, EDRN %in% higher_prob_samples)[, "EDRN"])), 
  id=c(paste("lesion"), "sampleType", "EDRN")) %>% 
  mutate(Disease_Free = rep(c("n", "y", "y"), length(rownames(.))/3))

ggplot(select_otu_data, aes(factor(sampleType, levels = c("initial", "follow_up")), 
                            log10(value+1.1), group = factor(EDRN))) + 
  geom_line(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
  geom_point(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
  facet_wrap(~variable, labeller = as_labeller(lesion_selected_labs)) + coord_cartesian(ylim = c(0, 4)) + 
  theme_bw() + ylab("Log Total Sequences") + xlab("") + ggtitle("Cancer") + 
  scale_colour_manual(name = "Cancer Free", labels = c("No", "Yes", "Unknown"), values = c("darkred", "orange", "red")) + 
  theme(legend.position=c(0.85,0.05), plot.margin = unit(c(1, 1, 1, 1), "lines"), 
        plot.title = element_text(size=20, face="bold"), legend.title = element_text(face="bold"))

  # This holds for the without fit model with the addition of an elevated Blautia OTU, another
  # Lachnospiraceae OTU, and a decrease in a Lactococcus OTU













