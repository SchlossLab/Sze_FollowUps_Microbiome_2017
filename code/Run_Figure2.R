### Create Figure 6 graph
### Show that CRC specific OTUs may be the difference between Adenoma and CRC
## Marc Sze

# Load needed functions
source('code/functions.R')

# Load needed packages
loadLibs(c("dplyr", "ggplot2", "gridExtra", "scales", "wesanderson"))

# Load needed data
graph_data <- read.csv("results/tables/adn_crc_maybe_diff.csv", header = T, stringsAsFactors = F)

# Add column to jitter values and italicize names
graph_data <- mutate(graph_data, 
    j_values = jitter(value, 3), 
    labs = c(
      rep(paste("paste(italic('Porphyromonas asaccharolytica '))", '(', unique(graph_data$otu)[1], ')', sep = ''), 
          length(rownames(graph_data))/length(unique(graph_data$tax_id))), 
      rep(paste("paste(italic('Fusobacterium nucleatum '))", '(', unique(graph_data$otu)[2], ')', sep = ''), 
          length(rownames(graph_data))/length(unique(graph_data$tax_id))), 
    rep(paste("paste(italic('Parvimonas micra '))", '(', unique(graph_data$otu)[3], ')', sep = ''), 
          length(rownames(graph_data))/length(unique(graph_data$tax_id))), 
    rep(paste("paste(italic('Peptostreptococcus stomatis '))", '(', unique(graph_data$otu)[4], ')', sep = ''), 
        length(rownames(graph_data))/length(unique(graph_data$tax_id)))))


# Create the figure

crc_specific <- grid.arrange(  
  # Cancer OTUs graph
  filter(graph_data, Dx_Bin == "cancer") %>% 
    ggplot(aes(factor(sampleType, levels = c("initial", "followup")), 
               j_values*100, group = factor(EDRN))) + 
    geom_line(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    geom_point(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    facet_wrap(~labs, labeller = label_parsed, scales = "free_y") + 
    theme_bw() + ylab("% Relative Abundance") + xlab("") + ggtitle("A") +  
    scale_colour_manual(name = "Cancer Free", 
                        label = c("No", "Yes", "Unknown"),  
                        values = wes_palette("GrandBudapest")) + 
    scale_x_discrete(
      breaks = c("initial", "followup"), 
      labels = c("Initial", "Follow Up")) + 
    theme(legend.title = element_text(face="bold", size = 8), 
          legend.text = element_text(size = 6), 
          legend.position = c(0.38, 0.82), 
          plot.margin = unit(c(1, 1, 1, 1), "lines"), 
          plot.title = element_text(size=20, face="bold"), 
          strip.text.x = element_text(size = 8)), 
  
  # Adenoma OTUs graph
  filter(graph_data, Dx_Bin != "cancer") %>% 
    ggplot(aes(factor(sampleType, levels = c("initial", "followup")), 
               j_values*100, group = factor(EDRN))) + 
    geom_line(aes(color = factor(Dx_Bin))) + geom_point(aes(color = factor(Dx_Bin))) +  
    facet_wrap(~labs, labeller = label_parsed, scales = "free_y") + 
    theme_bw() + ylab("% Relative Abundance") + xlab("") + ggtitle("B") + 
    scale_colour_manual(name = "Polyp Type", values = c("cyan", "blue"), 
                        breaks = c("adenoma", "adv_adenoma"), 
                        labels = c("Adenoma", "SRN")) +  
    scale_x_discrete(
      breaks = c("initial", "followup"), 
      labels = c("Initial", "Follow Up")) + 
    theme(legend.title = element_text(face="bold", size = 8), 
          legend.text = element_text(size = 6), 
          legend.position = c(0.1, 0.85), 
          plot.margin = unit(c(1, 1, 1, 1), "lines"), 
          plot.title = element_text(size=20, face="bold"), 
          strip.text.x = element_text(size = 8))

  )


ggsave(file = "results/figures/Figure2.pdf", crc_specific, 
       width=8.5, height = 11, dpi = 300)


# If wanting to use a log scale instead
# scale_y_log10(breaks = c(10, 1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001), na.value = 0.00000001)



