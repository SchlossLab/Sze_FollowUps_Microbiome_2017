## Explore Follow up data
## Marc Sze

###Load needed Libraries and functions
source('code/functions.R')
source('code/graphFunctions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "wesanderson", "vegan"))

### Read in necessary data 

tax <- read.delim('data/process/followUps.final.an.unique_list.0.03.cons.taxonomy', 
                  sep='\t', header=T, row.names=1)
shared <- read.delim('data/process/followUps.final.an.unique_list.0.03.subsample.0.03.filter.shared', 
                     header=T, sep='\t')

lesion_confirmed_vars <- read.csv('results/tables/lesion_confirmed_vars.csv', 
                                  header = TRUE, stringsAsFactors = F) %>% filter(otus != "fit_result")
SRNlesion_confirmed_vars <- read.csv('results/tables/SRNlesion_confirmed_vars.csv', 
                                     header = TRUE, stringsAsFactors = F) %>% filter(otus != "fit_result")
threeGroups_confirmed_vars <- read.csv('results/tables/threeGroups_confirmed_vars.csv', 
                                       header = TRUE, stringsAsFactors = F) %>% filter(otus != "fit_result")

# Convert taxa table to a data frame with columns for each taxonomic division
tax_df <- data.frame(do.call('rbind', strsplit(as.character(tax$Taxonomy), ';')))
rownames(tax_df) <- rownames(tax)
colnames(tax_df) <- c("Domain", "Phyla", "Order", "Class", "Family", "Genus")
tax_df <- as.data.frame(apply(tax_df, 2, function(x) gsub("\\(\\d{2}\\d?\\)", "", x)))
  # This uses the R regex to remove numerals and the brackets the \\ is used to escape specific function components
rm(tax)

# Create label lists for graphs to be used later
SRNlesion_selected_taxa <- tax_df[SRNlesion_confirmed_vars[, 'otus'], ]
SRNlesion_selected_labs <- createTaxaLabeller(SRNlesion_selected_taxa)

lesion_selected_taxa <- tax_df[lesion_confirmed_vars[, 'otus'], ]
lesion_selected_labs <- createTaxaLabeller(lesion_selected_taxa)

threeGroups_selected_taxa <- tax_df[threeGroups_confirmed_vars[, 'otus'], ]
threeGroups_selected_labs <- createTaxaLabeller(threeGroups_selected_taxa)


### Organize tables for train and test sets
metaF <- read.delim('data/process/followUps_metadata.txt', header=T, sep='\t')
metaI <- read.delim('data/process/initials_metadata.tsv', header=T, sep='\t')

### Need to amend and separate Adenoma and CRC
outcomes_f <- read.csv('data/process/followUp_outcome_data.csv', header = T, 
                       stringsAsFactors = F)

### Need to amend outcomes file so it matches the rest
outcomes_f <- outcomes_f[match(metaF$EDRN, outcomes_f$EDRN), ]

###Pull out the thetayc distance between the initial and followup sample within the
###same person for all and split by adenoma or cancer

thetaCompTotal <- dissplit(
  'data/process/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist', metaF)

intra_ade <- as.vector(unlist(thetaCompTotal['intra_ade']))
intra_canc <- as.vector(unlist(thetaCompTotal['intra_canc']))
inter <- as.vector(unlist(thetaCompTotal['inter']))
rm(thetaCompTotal)

###Pull out the share OTUs between the initial and followup sample within the
###same person for all and split by adenoma or cancer

sobsCompTotal <- dissplit(
  'data/process/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist', metaF)

sobs_intra_ade <- as.vector(unlist(sobsCompTotal['intra_ade']))
sobs_intra_canc <- as.vector(unlist(sobsCompTotal['intra_canc']))
sobs_inter <- as.vector(unlist(sobsCompTotal['inter']))
rm(sobsCompTotal)

# Get OTU abundances for initial samples that have followups
initial <- merge(metaF, shared, by.x='initial', by.y='Group')
initial <- initial[,c('fit_result',
                      colnames(initial)[grep('Otu[0123456789]', 
                                             colnames(initial))])]

# Get OTU abundances for follow ups
followups <- merge(metaF, shared, by.x='followUp', by.y='Group')
followups <- followups[,c('fit_result',
                          colnames(followups)[grep('Otu[0123456789]', 
                                                   colnames(followups))])]


cancer_out_f <- outcomes_f[which(outcomes_f$Diagnosis == "adenocarcinoma" | 
                                   outcomes_f$Diagnosis == "N/D"), ]

meta_cancer <- merge(metaF, cancer_out_f, by.x='EDRN', by.y='EDRN')

# Create a both treatment group for radiation and chemo
meta_cancer$both_treatments[which(
  (meta_cancer$chemo_received == "no" & meta_cancer$radiation_received == "yes") | 
    (meta_cancer$chemo_received == "yes" & meta_cancer$radiation_received == "no") | 
    (meta_cancer$chemo_received == "no" & meta_cancer$radiation_received == "no"))] <- "no"

meta_cancer$both_treatments[which(
  meta_cancer$chemo_received == "yes" & meta_cancer$radiation_received == "yes")] <- "yes"


# Run a PERMANOVA on thetayc distance matrix based on treatment breakdown

thetaCompTotal <- dissplit('data/process/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist',split=F, meta = F)

cancerOnlyTheta_follow <- filter(thetaCompTotal, rownames(thetaCompTotal) %in% as.character(meta_cancer$followUp))
cancerOnlyTheta_follow <- cancerOnlyTheta_follow[, colnames(cancerOnlyTheta_follow) %in% as.character(meta_cancer$followUp)]
rownames(cancerOnlyTheta_follow) <- colnames(cancerOnlyTheta_follow)

rownames(meta_cancer) <- as.character(meta_cancer$followUp)

meta_cancer <- meta_cancer[rownames(cancerOnlyTheta_follow), ]

set.seed(050416)
adonis(as.dist(cancerOnlyTheta_follow) ~ meta_cancer$chemo_received)

set.seed(050416)
adonis(as.dist(cancerOnlyTheta_follow) ~ meta_cancer$radiation_received)

set.seed(050416)
adonis(as.dist(cancerOnlyTheta_follow) ~ meta_cancer$both_treatments)


# Run and test change in thetayc and if that is different 
difference_table_treatment <- pickDistanceValues(vec1, vec2, thetaCompTotal, meta_cancer, 
                                                 c("chemo_received", "radiation_received", "both_treatments"))

t.test(distance ~ chemo_received, data = difference_table_treatment)
t.test(distance ~ radiation_received, data = difference_table_treatment)
t.test(distance ~ both_treatments, data = difference_table_treatment)


# Create statistics test for each OTU and create a summary table (LESION).
pValue_test_data_lesion <- rename(meta_cancer, Group = followUp) %>% inner_join(shared, by = "Group") %>% 
  select(Group, chemo_received, radiation_received, number_of_meds, both_treatments, 
         one_of(lesion_confirmed_vars$otus))

pValue_table_chemo_rad_lesion <- getCorrectedPvalue(
  pValue_test_data_lesion, lesion_confirmed_vars$otus, c("chemo_received", "radiation_received", "both_treatments"))
      # No significant differences between groups based on lesion important OTUs

# Create statistics test for each OTU and create a summary table (SRNLESION).
pValue_test_data_SRNlesion <- rename(meta_cancer, Group = followUp) %>% inner_join(shared, by = "Group") %>% 
  select(Group, chemo_received, radiation_received, number_of_meds, both_treatments, 
         one_of(SRNlesion_confirmed_vars$otus))

pValue_table_chemo_rad_SRNlesion <- getCorrectedPvalue(
  pValue_test_data_SRNlesion, SRNlesion_confirmed_vars$otus, c("chemo_received", "radiation_received", "both_treatments"))
    # No significant differences between groups based on lesion important OTUs


# Create statistics test for each OTU and create a summary table (THREEGROUPS).
pValue_test_data_threeGroups <- rename(meta_cancer, Group = followUp) %>% inner_join(shared, by = "Group") %>% 
  select(Group, chemo_received, radiation_received, number_of_meds, both_treatments, 
         one_of(threeGroups_confirmed_vars$otus))

pValue_table_chemo_rad_threeGroups <- getCorrectedPvalue(
  pValue_test_data_threeGroups, threeGroups_confirmed_vars$otus, c("chemo_received", "radiation_received", "both_treatments"))
# No significant differences between groups based on lesion important OTUs


# Make data table of graph for each model tested
specific_data_lesion <- melt(pValue_test_data_lesion, 
                      id = c("Group", "chemo_received", "radiation_received", "both_treatments", "number_of_meds"))

specific_data_SRNlesion <- melt(pValue_test_data_SRNlesion, 
                             id = c("Group", "chemo_received", "radiation_received", "both_treatments", "number_of_meds"))

specific_data_threeGroups <- melt(pValue_test_data_threeGroups, 
                                id = c("Group", "chemo_received", "radiation_received", "both_treatments", "number_of_meds"))

# Visualize breakdown by Chemo Recieved
    # mean_cl_boot obtains population confidence limits via nonparametric bootstraping for mean without assuming normality
    # e.g. smean.cl.boot(x, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)
    # smean.cl..... and mean.cl... are equivalent.
ggplot(specific_data_lesion, aes(factor(chemo_received), log10(value + 1), group = variable)) + 
  geom_jitter(aes(color=chemo_received), width = 0.5) + stat_summary(fun.y = mean, colour = "black", geom = "line") + 
  stat_summary(fun.data = "mean_cl_boot", colour = "black") + 
  facet_wrap(~variable, labeller = as_labeller(lesion_selected_labs)) + ylab("Subsampled Sequence Reads") + 
  xlab("") + theme_bw() + theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), 
                                axis.ticks.x = element_blank())

#Visualize breakdown by Both Received versus only one or none
ggplot(specific_data_lesion, aes(factor(both_treatments), log10(value + 1), group=variable)) + 
  geom_jitter(aes(color=both_treatments), width = 0.5) + stat_summary(fun.y = mean, colour = "black", geom = "line") + 
  stat_summary(fun.data = "mean_cl_boot", colour = "black") + 
  facet_wrap(~variable, labeller = as_labeller(lesion_selected_labs)) + ylab("Subsampled Sequence Reads") + 
  xlab("") + theme_bw() + theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), 
                                axis.ticks.x = element_blank())

# can add position=position_dodge(width=0.5) in the geom_jitter call instead of
# width but lose the variation at the 0 point
