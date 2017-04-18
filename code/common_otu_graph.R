### Plot for FHCRC
### Common OTUs
### Marc Sze

source('code/functions.R')
source('code/old/common_all_models.R')

loadLibs("gridExtra")

# Load needed data
common_vars_summary <- read.csv('data/process/tables/pvalue_adn_srn_crc_common_imp_vars.csv', 
                                header = T, row.names = 1, stringsAsFactors = F)

adn_varsMDA <- read.csv('data/process/tables/adn_reduced_model_top_vars_MDA_Summary.csv', 
                        header = T, stringsAsFactors = F)

srn_varsMDA <- read.csv('data/process/tables/srn_reduced_model_top_vars_MDA_Summary.csv', 
                        header = T, stringsAsFactors = F)

crc_varsMDA <- read.csv('data/process/tables/reduced_crc_model_top_vars_MDA_Summary.csv', 
                        header = T, stringsAsFactors = F)

# Create selection vector
common_vars <- common_vars_summary$otu

# Create MDA vectors
adn_rank <- (filter(adn_varsMDA, otu %in% common_vars) %>% slice(match(common_vars, otu)))[, "rank"]
srn_rank <- (filter(srn_varsMDA, otu %in% common_vars) %>% slice(match(common_vars, otu)))[, "rank"]
crc_rank <- (filter(crc_varsMDA, otu %in% common_vars) %>% slice(match(common_vars, otu)))[, "rank"]

# Create a combined table for graphing
combined_table <- select(common_vars_summary, otu) %>% 
  mutate(Tax = gsub("_unclassified", "", common_vars_summary$lowest_ID), 
         adn_rank = adn_rank, srn_rank = srn_rank, crc_rank = crc_rank) %>% 
  gather(key = model, value = rank, adn_rank, srn_rank, crc_rank)


otu_num <- as.numeric(gsub("Otu", "", common_vars))

test <- c()
for(i in 1:length(common_vars)){
  
  test <- c(test, 
            bquote(paste(italic(.(gsub("_", " ", combined_table$Tax[i]))) ~ "(OTU", .(otu_num[i]), ")", sep = "")))
}

common_labels <- do.call(expression, test)

# Create Graph
test_graph <- ggplot(combined_table, 
                     aes(factor(otu, levels = common_vars), rank)) + 
  geom_point(aes(color = factor(model, levels = c("adn_rank", "srn_rank", "crc_rank"), 
                             labels = c("Adenoma", "Advanced\nAdenoma", "Carcinoma"))), size = 4, alpha = 0.75) + 
  scale_y_continuous(trans = "reverse", breaks = c(1, 10, 20, 30, 40, 50, 60, 70)) + 
  scale_x_discrete(labels = common_labels) + 
  scale_color_manual(name = "Model", values = c('#228B22', '#FFD700', '#DC143C')) + 
  theme_bw() +  coord_flip(ylim = c(1, 70)) + 
  xlab("") + ylab("Importance Rank in Model") + ggtitle("B") + 
  theme(legend.title = element_text(face = "bold"), 
        axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", hjust = -0.63, size = 20),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray"))


combined_graph <- grid.arrange(venn_graph, test_graph)


ggsave(file = "results/figures/Figure4.pdf", combined_graph, 
       width=6, height = 9, dpi = 300)
