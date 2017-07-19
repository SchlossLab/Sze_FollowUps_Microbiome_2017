### Running Antibiotic comparisons 
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr"))


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
working_beta %>% filter(Dx_Bin == "adenoma") %>% group_by(Surgery) %>% 
  summarise(median_value = median(distance), 
            q25 = quantile(distance)["25%"], q75 = quantile(distance)["75%"])

wilcox.test(filter(working_beta, Dx_Bin == "adenoma" & Surgery == "Y")[, "distance"], 
            filter(working_beta, Dx_Bin == "adenoma" & Surgery == "N")[, "distance"])

working_beta %>% filter(Dx_Bin == "adv_adenoma") %>% group_by(Surgery) %>% 
  summarise(median_value = median(distance), 
            q25 = quantile(distance)["25%"], q75 = quantile(distance)["75%"])

wilcox.test(filter(working_beta, Dx_Bin == "adv_adenoma" & Surgery == "Y")[, "distance"], 
            filter(working_beta, Dx_Bin == "adv_adenoma" & Surgery == "N")[, "distance"])


## Make comparisons of only surgery changes based on diagnosis

wilcox.test(filter(working_beta, Dx_Bin == "adenoma" & Surgery == "Y")[, "distance"], 
            filter(working_beta, Dx_Bin == "cancer" & Surgery == "Y")[, "distance"])

wilcox.test(filter(working_beta, Dx_Bin == "adv_adenoma" & Surgery == "Y")[, "distance"], 
            filter(working_beta, Dx_Bin == "cancer" & Surgery == "Y")[, "distance"])

wilcox.test(filter(working_beta, Dx_Bin == "adenoma" & Surgery == "Y")[, "distance"], 
            filter(working_beta, Dx_Bin == "adv_adenoma" & Surgery == "Y")[, "distance"])


# Want to graph both sections and add a supplemental figure 



