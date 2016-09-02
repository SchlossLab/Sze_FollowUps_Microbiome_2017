## Explore Follow up data
## Marc Sze

# Load packages
library(ggplot2)
library(gridExtra)

### Read in necessary data 

tax <- read.delim('data/followUps.final.an.unique_list.0.03.cons.taxonomy', 
                  sep='\t', header=T, row.names=1)
shared <- read.delim('data/followUps.final.an.unique_list.0.03.subsample.0.03.filter.shared', 
                     header=T, sep='\t')

### Organize tables for train and test sets
metaF <- read.delim('data/followUps_metadata.txt', header=T, sep='\t')
metaF$lesion <- factor(NA, levels=c(0,1))
metaF$lesion[metaF$dx!='normal'] <-1
metaI <- read.delim('data/initials_metadata.tsv', header=T, sep='\t')
metaI$lesion <- factor(NA, levels=c(0,1))
metaI$lesion[metaI$dx=='normal'] <-0
metaI$lesion[metaI$dx!='normal'] <-1

### Need to amend and separate Adenoma and CRC
outcomes_f <- read.csv('data/followUp_outcome_data.csv', header = T, 
                       stringsAsFactors = F)

### Need to amend outcomes file so it matches the rest
outcomes_f <- outcomes_f[match(metaF$EDRN, outcomes_f$EDRN), ]

###Pull out the thetayc distance between the initial and followup sample within the
###same person for all and split by adenoma or cancer

source('code/functions.R')

thetaCompTotal <- dissplit(
  'data/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist', metaF)

intra_ade <- as.vector(unlist(thetaCompTotal['intra_ade']))
intra_canc <- as.vector(unlist(thetaCompTotal['intra_canc']))
inter <- as.vector(unlist(thetaCompTotal['inter']))
rm(thetaCompTotal)

###Pull out the share OTUs between the initial and followup sample within the
###same person for all and split by adenoma or cancer

sobsCompTotal <- dissplit(
  'data/followUps.final.an.unique_list.thetayc.0.03.lt.ave.dist', metaF)

sobs_intra_ade <- as.vector(unlist(sobsCompTotal['intra_ade']))
sobs_intra_canc <- as.vector(unlist(sobsCompTotal['intra_canc']))
sobs_inter <- as.vector(unlist(sobsCompTotal['inter']))
rm(sobsCompTotal)

# Get OTU abundances for initial samples that have followups
initial <- merge(metaF, shared, by.x='initial', by.y='Group')
initial <- initial[,c('lesion','fit_result',
                      colnames(initial)[grep('Otu[0123456789]', 
                                             colnames(initial))])]

# Get OTU abundances for follow ups
followups <- merge(metaF, shared, by.x='followUp', by.y='Group')
followups <- followups[,c('lesion','fit_result',
                          colnames(followups)[grep('Otu[0123456789]', 
                                                   colnames(followups))])]


cancer_out_f <- outcomes_f[which(outcomes_f$Diagnosis == "adenocarcinoma" | 
                                   outcomes_f$Diagnosis == "N/D"), ]

meta_cancer <- merge(metaFConly, cancer_out_f, by.x='EDRN', by.y='EDRN')


### Grab specific OTUs and fit results that were selected as important features for classification
### Main goal is to compare how different treatments may affect the follow up microbiome samples


taxa_int <- c('Otu000034' = "Collinsella", 
              'Otu000059' = "Clostridiales(1)", 
              'Otu000062' = "Odoribacter", 
              'Otu000084' = "Bilophila", 
              'Otu000103' = "Lachnospiraceae", 
              'Otu000126' = "P.asaccharolytica",
              'Otu000145' = "Enterorhabdus", 
              'Otu000147' = "Prevotella", 
              'Otu000205' = "F.nucleatum", 
              'Otu000217' = "Coriobacteriaceae", 
              'Otu000333' = "Fusobacterium", 
              'Otu000381' = "Firmicutes", 
              'Otu000397' = "P. micra", 
              'Otu000563' = "Porphyromonas", 
              'Otu000566' = "P.stomatis",
              'Otu000671' = "Dialister", 
              'Otu000902' = "Clostridiales(2)")


select_taxa <- unname(taxa_int)
select_otus <- names(taxa_int)

meta_cancer$both_treatments[which(
  (meta_cancer$chemo_received == "no" & meta_cancer$radiation_received == "yes") | 
    (meta_cancer$chemo_received == "yes" & meta_cancer$radiation_received == "no") | 
    (meta_cancer$chemo_received == "no" & meta_cancer$radiation_received == "no"))] <- "no"

meta_cancer$both_treatments[which(
  meta_cancer$chemo_received == "yes" & meta_cancer$radiation_received == "yes")] <- "yes"


# Make data table of graph 
specific_data <- rename(meta_cancer, Group = followUp) %>% inner_join(shared, by = "Group") %>% 
  select(Group, chemo_received, radiation_received, number_of_meds, both_treatments, one_of(select_otus)) %>% 
  melt(id = c("Group", "chemo_received", "radiation_received", "both_treatments", "number_of_meds"))

# Visualize breakdown by Chemo Recieved
    # mean_cl_boot obtains population confidence limits via nonparametric bootstraping for mean without assuming normality
    # e.g. smean.cl.boot(x, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)
    # smean.cl..... and mean.cl... are equivalent.
ggplot(specific_data, aes(factor(chemo_received), log10(value + 1), group = variable)) + 
  geom_jitter(aes(color=chemo_received), width = 0.5) + stat_summary(fun.y = mean, colour = "black", geom = "line") + 
  stat_summary(fun.data = "mean_cl_boot", colour = "black") + 
  facet_wrap(~variable, labeller = as_labeller(taxa_int)) + ylab("Subsampled Sequence Reads") + 
  xlab("") + theme_bw() + theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), 
                                axis.ticks.x = element_blank())

# Visualize breakdown by Radiation Received
ggplot(specific_data, aes(factor(radiation_received), log10(value + 1), group = variable)) + 
  geom_jitter(aes(color=radiation_received), width = 0.5) + stat_summary(fun.y = mean, colour = "black", geom = "line") + 
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +
  facet_wrap(~variable, labeller = as_labeller(taxa_int)) + ylab("Subsampled Sequence Reads") + 
  xlab("") + theme_bw() + theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), 
                                axis.ticks.x = element_blank())

#Visualize breakdown by Both Received versus only one or none
ggplot(specific_data, aes(factor(both_treatments), log10(value + 1), group=variable)) + 
  geom_jitter(aes(color=both_treatments), width = 0.5) + stat_summary(fun.y = mean, colour = "black", geom = "line") + 
  stat_summary(fun.data = "mean_cl_boot", colour = "black") + 
  facet_wrap(~variable, labeller = as_labeller(taxa_int)) + ylab("Subsampled Sequence Reads") + 
  xlab("") + theme_bw() + theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), 
                                axis.ticks.x = element_blank())

# can add position=position_dodge(width=0.5) in the geom_jitter call instead of
# width but lose the variation at the 0 point
