### Create Figure S3 graph
### Show results for initial and follow up samples based on 3 specific OTUs
## Marc Sze


#Load needed libraries
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2"))


# Read in necessary data frames
shared <- read.delim('data/process/final.shared', 
                     header=T, sep='\t') %>% select(Group, contains("Otu0"))
         
         
good_metaf <- read.csv('data/process/mod_metadata/good_metaf_final.csv', 
                       header = T, stringsAsFactors = F)

reduced_model <- read.csv('data/process/tables/reduced_crc_model_top_vars_MDA_Summary.csv', 
                          header = T, stringsAsFactors = F)

crc_tax <- read.csv('data/process/tables/crc_rf_otu_tax.csv', 
                    header = T, stringsAsFactors = F) %>% rename(otu = X)

# select out only necessary groups
test <- shared %>% slice(match(c(good_metaf$initial, good_metaf$followUp), Group))

# Get Relative Abundance
total_sub_Seqs <- test %>% select(contains("Otu0")) %>% rowSums()

test2 <- as.data.frame(apply(select(test, -Group), 2, function(x) x/total_sub_Seqs)) %>% 
  mutate(Group = test$Group)

# Select out for only three specific OTUs 
otus <- filter(crc_tax, Genus == "Porphyromonas" | Genus == "Parvimonas" | Genus == "Fusobacterium")[, "otu"]

# Filter shared by these specific otus
test2 <- test2 %>% select(Group, one_of(otus))

# Create modified data table
data_analysis <- data_frame(Group = c(good_metaf$initial, good_metaf$followUp), 
                            EDRN = rep(good_metaf$EDRN, 2), 
                            sampleType = c(rep("initial", length(good_metaf$EDRN)), rep("followup", length(good_metaf$EDRN))), 
                            Dx_Bin = rep(good_metaf$Dx_Bin, 2), 
                            dx = rep(good_metaf$dx, 2)) %>% 
  inner_join(test2, by = "Group") %>% 
  gather(key = otu, value = rel.abund, -EDRN, -Group, -sampleType, -Dx_Bin, -dx)



# Pull OTUs that are only in the MDA data
select_tax_df <- (crc_tax %>% slice(match(otus, otu)))[, "Genus"]
otu_num <- as.numeric(gsub("Otu", "", otus))

# create labels for factor values with low taxonomy
test <- c()
for(i in 1:length(select_tax_df)){
  
  test <- c(test, bquote(paste(italic(.(select_tax_df[i])) ~ "(OTU", .(otu_num[i]), ")", sep = "")))
  
}

label_names <- do.call(expression, test)
names(label_names) <- otus



oral_bugs <- filter(data_analysis, Dx_Bin == "cancer") %>% 
  ggplot(aes(factor(sampleType, levels = c("initial", "followup"), 
                    labels = c("Pre", "Post")), rel.abund, group = EDRN)) + 
  geom_point() + geom_line() + ylab("Relative Abundance") + xlab("") + 
  facet_wrap(~factor(otu, levels = otus, labels = label_names), labeller = label_parsed, scales = "free_y") + theme_bw()


ggsave(file = "results/figures/FigureS3.pdf", oral_bugs, 
       width=11, height = 6, dpi = 300)











