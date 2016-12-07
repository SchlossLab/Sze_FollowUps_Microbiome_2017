### Create Figure 6 graph
### Show that CRC specific OTUs may be the difference between Adenoma and CRC
## Marc Sze

# Load needed functions
source('code/functions.R')
source('code/graphFunctions.R')

# Load needed packages
loadLibs(c("dplyr", "ggplot2", "reshape2", "gridExtra", "scales", "wesanderson", "knitr", "rmarkdown"))

# Load needed data
graph_data <- read.csv("results/tables/adn_crc_maybe_diff.csv", header = T, stringsAsFactors = F)

graph_data <- mutate(graph_data, j_values = jitter(value))

labs <- c(
  paste("Fusobacterium Nucleatum ", "(", unique(graph_data$otu)[1], ")", sep = ""), 
  paste("Parvimonas Micra ", "(", unique(graph_data$otu)[2], ")", sep = ""), 
  paste("Peptostreptococcus Assacharolytica ", "(", unique(graph_data$otu)[3], ")", sep = ""), 
  paste("Porphyromonas Stomatis ", "(", unique(graph_data$otu)[4], ")", sep = ""))

names(labs) <- unique(graph_data$otu)

# Create the figure

crc_specific <- grid.arrange(  
  # Cancer OTUs graph
  filter(graph_data, Dx_Bin == "cancer") %>% 
    ggplot(aes(factor(sampleType, levels = c("initial", "followup")), 
               j_values, group = factor(EDRN))) + 
    geom_line(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    geom_point(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    facet_wrap(~otu, labeller = as_labeller(labs)) + coord_cartesian(ylim = c(0, 0.105)) + 
    theme_bw() + ylab("Relative Abundance") + xlab("") + ggtitle("A") + 
    scale_colour_manual(name = "Cancer Free", 
                        label = c("No", "Yes", "Unknown"),  
                        values = wes_palette("GrandBudapest")) + 
    scale_x_discrete(
      breaks = c("initial", "followup"), 
      labels = c("Initial", "Follow Up")) + 
    theme(legend.title = element_text(face="bold"), 
          legend.position = c(0.40, 0.30), 
          plot.margin = unit(c(1, 1, 1, 1), "lines"), 
          plot.title = element_text(size=20, face="bold")), 
  
  # Adenoma OTUs graph
  filter(graph_data, Dx_Bin != "cancer") %>% 
    ggplot(aes(factor(sampleType, levels = c("initial", "followup")), 
               j_values, group = factor(EDRN))) + 
    geom_line(aes(color = factor(Dx_Bin))) + geom_point(aes(color = factor(Dx_Bin))) +  
    facet_wrap(~otu, labeller = as_labeller(labs)) + coord_cartesian(ylim = c(0, 0.002)) + 
    theme_bw() + ylab("Relative Abundance") + xlab("") + ggtitle("B") + 
    scale_colour_manual(name = "Polyp Type", values = c("cyan", "blue"), 
                        breaks = c("adenoma", "adv_adenoma"), 
                        labels = c("Adenoma", "SRN")) +  
    scale_x_discrete(
      breaks = c("initial", "followup"), 
      labels = c("Initial", "Follow Up")) + 
    theme(legend.title = element_text(face="bold"), 
          legend.position = c(0.4, 0.30), 
          plot.margin = unit(c(1, 1, 1, 1), "lines"), 
          plot.title = element_text(size=20, face="bold"))

  )


ggsave(file = "results/figures/Figure6.pdf", crc_specific, 
       width=8.5, height = 11, dpi = 300)







