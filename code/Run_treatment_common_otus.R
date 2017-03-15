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
total_seqs <- rowSums(select(shared, -Group))

shared <- cbind(Group = shared$Group, 
                as.data.frame(apply(select(shared, -Group), 2, 
                                    function(x) x/total_seqs)))

# Shrink shared file down to only the samples initial then follow ups
samples_to_keep <- as.character(c(good_metaf$initial, good_metaf$followUp))
shared <- filter(shared, Group %in% samples_to_keep) %>% 
  slice(match(samples_to_keep, Group))

# Select only OTUs that were common
shared <- shared %>% select(one_of(common_vars_summary$otu)) %>% 
  mutate(sampleType = c(rep("initial", length(good_metaf$initial)), rep("followups", length(good_metaf$followUp))))

# Calculate P-values and write out to csv
otus <- unique(common_vars_summary$otu)

pvalue_table <- as.data.frame(matrix(nrow = length(common_lowest_IDs), 
                                     ncol = 4, dimnames = list(
                                       rown = otus, 
                                       coln = c("Pvalue_chemo", "bh_chemo", "Pvalue_rads", "bh_rads"))))

for(i in 1:length(common_lowest_IDs)){
  
  pvalue_table[i, "Pvalue_chemo"] <- wilcox.test(filter(common_data, chemo == "yes" & otu == otus[i])[, "change"], 
              filter(common_data, chemo == "no" & otu == otus[i])[, "change"])$p.value
  
  pvalue_table[i, "Pvalue_rads"] <- wilcox.test(filter(common_data, rads == "yes" & otu == otus[i])[, "change"], 
                                                 filter(common_data, rads == "no" & otu == otus[i])[, "change"])$p.value
}

pvalue_table <- pvalue_table %>% mutate(bh_chemo = p.adjust(Pvalue_chemo, method = "BH"), 
                                        bh_rads = p.adjust(Pvalue_rads, method = "BH"), 
                                        lowest_ID = common_lowest_IDs, 
                                        otu = otus)

write.csv(pvalue_table, "data/process/tables/chemo_rads_treatment_pvalue_summary.csv", row.names = F)

# Create lowest ID labels
common_lowest_IDs <- common_vars_summary$lowest_ID
table_labels <- c()
for(i in 1:length(common_lowest_IDs)){
  
  table_labels <- c(table_labels, rep(common_lowest_IDs[i], length(good_metaf$EDRN)))
}

# Create data table to be used for graphing
common_data <- ((filter(shared, sampleType == "initial") %>% select(-sampleType)) -
  (filter(shared, sampleType == "followups") %>% select(-sampleType))) %>% 
  mutate(chemo = good_metaf$chemo_received, rads = good_metaf$radiation_received, 
         Disease_Free = good_metaf$Disease_Free, Dx_Bin = good_metaf$Dx_Bin) %>% 
  gather(key = otu, value = change, one_of(common_vars_summary$otu)) %>% 
  mutate(lowest_ID = table_labels)

# Create graph labels
test <- c()
for(i in 1:length(common_vars_summary$otu)){
  
  test <- c(test, bquote(paste(
    italic(.(gsub("_", " ", common_vars_summary$lowest_ID)[i])) ~ "(", .(common_vars_summary$otu[i]), ")", sep = "")))
}

test2 <- do.call(expression, test)


final_graph <- grid.arrange(
  #IF Probability change after Chemo received
  ggplot(chemo_rad_summary, 
         aes(factor(chemo, levels = c("yes", "no"), labels = c("Yes", "No")), 
             red_IF, group = 1)) + 
    geom_jitter(aes(color = chemo), width = 0.15, size = 4) + 
    scale_color_manual(values = c("#7466ff", "#55b64e")) + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + 
    ylab("Reduction in Initial Sample Positive Probability") + 
    xlab("Chemotherapy Received") + ggtitle("A") + 
    coord_cartesian(ylim = c(0.4, 0.8)) + theme_bw() + 
    theme(plot.title = element_text(face = "bold", hjust = 0), 
          legend.position = "none", 
          axis.title = element_text(face="bold", hjust = 0.5, size = 8)), 
  
  # IF Probability change after Rads received
  ggplot(chemo_rad_summary, 
         aes(factor(rads, levels = c("yes", "no"), labels = c("Yes", "No")), 
             red_IF, group = 1)) + 
    geom_jitter(aes(color = rads), width = 0.15, size = 4) + 
    scale_color_manual(values = c("#7466ff", "#ff6700")) + 
    stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1) + 
    stat_summary(fun.y = mean, colour = "black", geom = "line") + 
    ylab("Reduction in Initial Sample Positive Probability") + 
    xlab("Radiation Received") + ggtitle("B") + 
    coord_cartesian(ylim = c(0.4, 0.8)) + theme_bw() + 
    theme(plot.title = element_text(face = "bold", hjust = 0), 
          legend.position = "none", 
          axis.title = element_text(face="bold", hjust = 0.5, size = 8)), 
  
  # Create graph for the common OTUs with chemotherapy
  ggplot(common_data, aes(otu, change, color = chemo)) + geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = c(0.1, -0.1)) + theme_bw() + 
    scale_x_discrete(labels = rev(test2)) + ylab("Follow Up Relative Abundance Change") + 
    xlab("") + ggtitle("C") + 
    scale_color_manual(values = c("#7466ff", "#55b64e"), name="Chemotherapy", 
                       breaks = c("yes", "no"), 
                       labels = c("Yes", "No")) + 
    theme(plot.title = element_text(face = "bold", hjust = 0), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.text.x = element_text(angle = 90, size = 5, hjust = 1, vjust = 0.5), 
          axis.title = element_text(face="bold", hjust = 0.5, size = 8), 
          legend.title = element_text(face="bold", size = 7), 
          legend.text = element_text(size = 7)) + 
    annotate("text", label = paste("*"), x = 31, y = 0.01, size = 6), 
  
  # Create graph for the common OTUs with radiation
  ggplot(common_data, aes(otu, change, color = rads)) + geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = c(0.1, -0.1)) + theme_bw() + 
    scale_x_discrete(labels = rev(test2)) + ylab("Follow Up Relative Abundance Change") + 
    xlab("") + ggtitle("D") + 
    scale_color_manual(values = c("#7466ff", "#ff6700"), name="Radiation", 
                       breaks = c("yes", "no"), 
                       labels = c("Yes", "No")) + 
    theme(plot.title = element_text(face = "bold", hjust = 0), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.text.x = element_text(angle = 90, size = 5, hjust = 1, vjust = 0.5), 
          axis.title = element_text(face="bold", hjust = 0.5, size = 8), 
          legend.title = element_text(face="bold", size = 7), 
          legend.text = element_text(size = 7)), layout_matrix = cbind(c(1, 3, 4), c(2, 3, 4)))





ggsave(file = "results/figures/common_treatment.pdf", final_graph, 
       width=10, height = 10, dpi = 300)













