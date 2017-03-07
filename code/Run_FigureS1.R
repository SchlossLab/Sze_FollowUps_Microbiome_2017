### Figure S1
### Distribution of P-values from Paired Wilcoxson Analysis
## Marc Sze

#Load needed code and packages
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "reshape2", "gridExtra", "scales", "wesanderson", "knitr", "rmarkdown"))

#Read data needed
paired_table <- read.csv("data/process/tables/OTU_paired_wilcoxson_test.csv", header = T, stringsAsFactors = F)

label_names <- c("Lesion", "Adenoma Only", "CRC Only")
names(label_names) <- c("ALL", "adn", "crc")

pvalue_distribution_paired <- ggplot(paired_table, aes(x=BH_corrected)) + geom_histogram() + 
  geom_vline(aes(xintercept = 0.05, color = "red"), linetype = 2, size = 1) + 
  facet_wrap(~factor(analysis, levels = c("ALL", "adn", "crc")), 
             labeller = as_labeller(label_names)) + theme_bw() + 
  xlab("Benjamini-Hochberg Corrected P-values") + ylab("Counts") + 
  theme(legend.position = "none")

ggsave(file = "results/figures/FigureS1.pdf", pvalue_distribution_paired, 
       width=12, height = 6, dpi = 300)
