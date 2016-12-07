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


labs <- c(
  paste("Fusobacterium Nucleatum ", "(", unique(graph_data$otu)[1], ")", sep = ""), 
  paste("Parvimonas Micra ", "(", unique(graph_data$otu)[2], ")", sep = ""), 
  paste("Peptostreptococcus Assacharolytica ", "(", unique(graph_data$otu)[3], ")", sep = ""), 
  paste("Porphyromonas Stomatis ", "(", unique(graph_data$otu)[4], ")", sep = ""))

names(labs) <- unique(graph_data$otu)

# Create the figure

grid.arrange(  
  # Adenoma OTUs graph
  filter(graph_data, Diagnosis != "cancer") %>% 
    ggplot(aes(factor(sampleType, levels = c("initial", "follow_up")), 
               log10(jitteredValues + 1.1), group = factor(EDRN))) + 
    geom_line(aes(color = factor(Disease_Free))) + geom_point(aes(color = factor(Disease_Free))) +  
    facet_wrap(~variable, labeller = as_labeller(lesion_selected_labs)) + coord_cartesian(ylim = c(0, 4)) + 
    theme_bw() + ylab("Log Total Sequences") + xlab("") + ggtitle("Adenoma") + 
    scale_colour_manual(values = "blue") +  
    theme(legend.position="none", plot.margin = unit(c(1, 1, 1, 1), "lines"), 
          plot.title = element_text(size=20, face="bold")), 
  # Cancer OTUs graph
  filter(select_otu_data, Diagnosis == "adenocarcinoma" | Diagnosis == "N/D") %>% 
    ggplot(aes(factor(sampleType, levels = c("initial", "follow_up")), 
               log10(jitteredValues+1.1), group = factor(EDRN))) + 
    geom_line(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    geom_point(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    facet_wrap(~variable, labeller = as_labeller(lesion_selected_labs)) + coord_cartesian(ylim = c(0, 4)) + 
    theme_bw() + ylab("Log Total Sequences") + xlab("") + ggtitle("Cancer") + 
    scale_colour_manual(name = "Cancer Free", labels = c("No", "Yes", "Unknown"), values = c("darkred", "orange", "red")) + 
    theme(legend.position=c(0.9,0.8), plot.margin = unit(c(1, 1, 1, 1), "lines"), 
          plot.title = element_text(size=20, face="bold"), legend.title = element_text(face="bold")))










