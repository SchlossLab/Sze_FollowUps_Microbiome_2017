### Running Antibiotic comparisons 
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "ggplot2", "gridExtra"))


# Read in data tables
good_metaf <- read.csv(
  "data/process/mod_metadata/metaF_final.csv", 
  stringsAsFactors = F, header = T) %>% 
  mutate(lesion_follow = ifelse(Disease_Free == "n", 1, 0))

# Load diversity data
alpha_summary <- read.delim("data/process/final.groups.ave-std.summary", stringsAsFactors = F)
thetaCompTotal <- read.dist('data/process/final.thetayc.0.03.lt.ave.dist')

# Create variables that are to be used 
initial_samples <- good_metaf$initial
followup_samples <- good_metaf$followUp

# Create alpha table for comparisons
working_alpha <- alpha_summary %>% filter(method == "ave") %>% 
  slice(match(c(initial_samples, followup_samples), group)) %>% 
  select(group, sobs, shannon, shannoneven)

# Modify meta data so that it is in a correct format
mod_metaf <- good_metaf %>% slice(match(initial_samples, initial)) %>% 
  bind_rows(slice(good_metaf, match(followup_samples, followUp))) %>% 
  mutate(sampleType = c(rep("initial", 67), rep("followup", 67))) %>% 
  select(-followUp) %>% rename(samples = initial)

# Select specific beta function
get_beta <- function(i, j, table_name){
  
  data_table <- get(table_name)
  beta_values <- data_table[i, j]
  
  return(beta_values)
}

# create the beta data and join with original metadata
working_beta <- as.data.frame(initial_samples) %>% 
  mutate(distance = as.numeric(
    mapply(get_beta, as.character(initial_samples), as.character(followup_samples), "thetaCompTotal"))) %>% 
  rename(initial = initial_samples) %>% 
  inner_join(good_metaf, by = "initial")



## Make comparisons of surgery for adenoma and advanced adenoma within groups
surg_non_crc_summary <- working_beta %>% filter(Dx_Bin == "adenoma") %>% group_by(Surgery) %>% 
  summarise(median_value = median(distance), 
            q25 = quantile(distance)["25%"], q75 = quantile(distance)["75%"]) %>% 
  mutate(dx = c("adenoma", "adenoma")) %>% 
  bind_rows(working_beta %>% filter(Dx_Bin == "adv_adenoma") %>% group_by(Surgery) %>% 
              summarise(median_value = median(distance), 
                        q25 = quantile(distance)["25%"], q75 = quantile(distance)["75%"]) %>% 
              mutate(dx = c("adv_adenoma", "adv_adenoma")))

non_crc_surg_pvalues <- data_frame(pvalues = c(wilcox.test(filter(working_beta, Dx_Bin == "adenoma" & Surgery == "Y")[, "distance"], 
            filter(working_beta, Dx_Bin == "adenoma" & Surgery == "N")[, "distance"])$p.value, 
            wilcox.test(filter(working_beta, Dx_Bin == "adv_adenoma" & Surgery == "Y")[, "distance"], 
            filter(working_beta, Dx_Bin == "adv_adenoma" & Surgery == "N")[, "distance"])$p.value)) %>% 
  mutate(BH = p.adjust(pvalues, method = "BH"))



## Make comparisons of only surgery changes based on diagnosis
surg_only_comparisons <- data_frame(
  pvalues = c(
    wilcox.test(filter(working_beta, Dx_Bin == "adenoma" & Surgery == "Y")[, "distance"], 
                filter(working_beta, Dx_Bin == "cancer" & Surgery == "Y")[, "distance"])$p.value, 
    wilcox.test(filter(working_beta, Dx_Bin == "adv_adenoma" & Surgery == "Y")[, "distance"], 
                filter(working_beta, Dx_Bin == "cancer" & Surgery == "Y")[, "distance"])$p.value, 
    wilcox.test(filter(working_beta, Dx_Bin == "adenoma" & Surgery == "Y")[, "distance"], 
                filter(working_beta, Dx_Bin == "adv_adenoma" & Surgery == "Y")[, "distance"])$p.value)) %>% 
  mutate(BH = p.adjust(pvalues, method = "BH"), comparison = c("adnVcrc", "srnVcrc", "adnVsrn"))

# Write out needed tables
write.csv(surg_only_comparisons, "data/process/tables/surgery_only_comparisons.csv", 
          row.names = F)

write.csv(surg_non_crc_summary, "data/process/tables/surgVnonsurg_noncrc.csv", 
          row.names = F)

write.csv(non_crc_surg_pvalues, "data/process/tables/surgVnonsurg_noncrc_pvalues.csv", 
          row.names = F)


# Want to graph both sections and add a supplemental figure 

surgery_only_plot <- working_beta %>% filter(Surgery == "Y") %>% 
ggplot(aes(factor(Dx_Bin, levels = c("adenoma", "adv_adenoma", "cancer"), 
                  labels = c("Adenoma", "Advanced\nAdenoma", "Carcinoma")), distance, group = 1)) + 
  geom_jitter(aes(color=Dx_Bin), width = 0.3, size = 4, alpha = 0.5) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) + 
  scale_color_manual(name = "Lesion Type", 
                     values = c('#228B22', '#FFD700', '#DC143C'), 
                     breaks = c("adenoma", "adv_adenoma", "cancer"), 
                     labels = c("Adenoma", "Advanced\nAdenoma", "Carcinoma")) + 
  coord_cartesian(ylim = c(0, 1)) + ylab(expression(paste(theta[YC], " Distance"))) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), 
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),limits = c(-0.06, 1.06)) + 
  xlab("") + theme_bw() + ggtitle("A") + 
  theme(axis.title = element_text(face="bold", hjust = 0.5), 
        axis.text.x = element_text(size = 12), 
        legend.title = element_text(face="bold"), 
        legend.position = "none", 
        plot.title = element_text(face="bold", hjust = -0.065, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


adenoma_surg_diff <- working_beta %>% filter(Dx_Bin == "adenoma") %>% 
  ggplot(aes(factor(Surgery, levels = c("Y", "N"), labels = c("Yes", "No")), distance)) + 
  geom_jitter(aes(color = Surgery), width = 0.3, size = 4, alpha = 0.5) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) + 
  scale_color_manual(name = "Surgery Received", 
                     values = c('#00EE76', '#228B22'), 
                     breaks = c("Y", "N"), 
                     labels = c("Yes", "No")) + 
  coord_cartesian(ylim = c(0, 1)) + ylab(expression(paste(theta[YC], " Distance"))) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), 
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),limits = c(-0.06, 1.06)) + 
  xlab("") + theme_bw() + ggtitle("B") + 
  theme(axis.title = element_text(face="bold", hjust = 0.5), 
        axis.text.x = element_text(size = 12), 
        legend.title = element_text(face="bold"), 
        legend.position = "none", 
        plot.title = element_text(face="bold", hjust = -0.065, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



adv_adenoma_surg_diff <- working_beta %>% filter(Dx_Bin == "adv_adenoma") %>% 
  ggplot(aes(factor(Surgery, levels = c("Y", "N"), labels = c("Yes", "No")), distance)) + 
  geom_jitter(aes(color = Surgery), width = 0.3, size = 4, alpha = 0.5) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) + 
  scale_color_manual(name = "Surgery Received", 
                     values = c('#CDAD00', '#8B7500'), 
                     breaks = c("Y", "N"), 
                     labels = c("Yes", "No")) + 
  coord_cartesian(ylim = c(0, 1)) + ylab(expression(paste(theta[YC], " Distance"))) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), 
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),limits = c(-0.06, 1.06)) + 
  xlab("") + theme_bw() + ggtitle("C") + 
  theme(axis.title = element_text(face="bold", hjust = 0.5), 
        axis.text.x = element_text(size = 12), 
        legend.title = element_text(face="bold"), 
        legend.position = "none", 
        plot.title = element_text(face="bold", hjust = -0.065, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


test <- grid.arrange(surgery_only_plot, adenoma_surg_diff, adv_adenoma_surg_diff, 
                     layout_matrix = cbind(c(1, 2), c(1, 3)))


ggsave(file = "results/figures/FigureS1.pdf", test, 
       width=10, height = 8.5, dpi = 300)
