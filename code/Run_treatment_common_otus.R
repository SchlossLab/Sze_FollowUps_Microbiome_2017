### Create data on treatment affects on common OTUs from each model
### Could these treatments influence findings of surgical removal
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "gridExtra"))

# Read in data tables
chemo_rad_stats_summary <- read.csv('data/process/tables/probs_chemo_rad_pvalue_summary.csv', 
                              header = T, stringsAsFactors = F) %>% rename(comparison = X)

good_metaf <- read.csv('data/process/mod_metadata/good_metaf_final.csv', header = T, stringsAsFactors = F)

common_vars_summary <- read.csv('data/process/tables/pvalue_IF_lesion_common_imp_vars.csv', 
                                header = T, row.names = 1, stringsAsFactors = F)

chemo_rad_summary <- read.csv('data/process/tables/chemo_rad_summary.csv', header = T, stringsAsFactors = F)

shared <- read.delim('data/process/final.shared', header = T, stringsAsFactors = F) %>% 
  select(-label, -numOtus) %>% mutate(Group = as.character(Group))

# Convert to relative abundance


# Shrink shared file down to only the samples initial then follow ups
samples_to_keep <- as.character(c(good_metaf$initial, good_metaf$followUp))
shared <- filter(shared, Group %in% samples_to_keep) %>% 
  slice(match(samples_to_keep, Group))

# Select only OTUs that were common
shared <- shared %>% select(one_of(common_vars_summary$otu))

# First graph of the panel: Chemo received

ggplot(chemo_rad_summary, 
       aes(factor(chemo, levels = c("yes", "no"), labels = c("Yes", "No")), 
           red_IF, group = 1)) + 
  geom_jitter(aes(color = chemo), width = 0.15, size = 4) + 
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1) + 
  stat_summary(fun.y = mean, colour = "black", geom = "line") + 
  ylab("Reduction in Initial Sample Positive Probability") + 
  xlab("Chemotherapy Received") + 
  coord_cartesian(ylim = c(0.4, 0.8)) + theme_bw() + 
  theme(legend.position = "none", 
        axis.title = element_text(face="bold", hjust = 0.5))


# Second graph of the panel: Rads received

ggplot(chemo_rad_summary, 
       aes(factor(rads, levels = c("yes", "no"), labels = c("Yes", "No")), 
           red_IF, group = 1)) + 
  geom_jitter(aes(color = rads), width = 0.15, size = 4) + 
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1) + 
  stat_summary(fun.y = mean, colour = "black", geom = "line") + 
  ylab("Reduction in Initial Sample Positive Probability") + 
  xlab("Radiation Received") + 
  coord_cartesian(ylim = c(0.4, 0.8)) + theme_bw() + 
  theme(legend.position = "none", 
        axis.title = element_text(face="bold", hjust = 0.5))





























