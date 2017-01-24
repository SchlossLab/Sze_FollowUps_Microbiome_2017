### List of Stored Functions for graph creation
### CRC Follow up project
### Marc Sze


# This function creates a visual representation of initial and follow up samples important variables for a chosen model
# These are at initial and follow up sampling.
follow_Abundance_Fit_Graph <- function(taxaFile, labelFile, followUpFile, initialData, followData, metadata, 
                                       jitterAmount = 0.08, n, modelOFinterest, rowOFinterest = "Disease_Free", 
                                       positiveTitle = "Cancer", fit = FALSE){
  
  loadLibs(c("dplyr", "ggplot2", "reshape2", "wesanderson"))
  
  jitter_values <- jitter(rep(0, n), amount = 0.08)
  
  #creating the data tables to be used
  select_otu_data <- melt(
    rbind(select(initialData, starts_with(paste(modelOFinterest)), one_of(rownames(taxaFile))) %>% 
            mutate(sampleType = "initial", EDRN = metadata$EDRN, Diagnosis = metadata$Diagnosis), 
          select(followData, starts_with(paste(modelOFinterest)), one_of(rownames(taxaFile))) %>% 
            mutate(sampleType = "follow_up", EDRN = metadata$EDRN, Diagnosis = metadata$Diagnosis)), 
    id=c(paste(modelOFinterest), "sampleType", "Diagnosis", "EDRN")) %>% 
    mutate(jitteredValues = as.numeric(value) + rep(jitter_values, length(rownames(.))/n)) %>% 
    mutate(Disease_Free = rep(metadata[, rowOFinterest], length(rownames(.))/n))
  
  
  fit_results <- select(metadata, starts_with(modelOFinterest), EDRN, Diagnosis, fit_result, fit_followUp) %>% 
    rename(initial = fit_result, follow_up = fit_followUp) %>% 
    melt(id=c(paste(modelOFinterest), "EDRN", "Diagnosis")) %>% 
    mutate(jitteredValues = as.numeric(value) + rep(jitter_values, 2)) %>% 
    mutate(Disease_Free = rep(metadata[, rowOFinterest], 2))
  
  # Create graphs to be used
  
  # Adenoma OTUs graph
  adenomaOTUs <- filter(select_otu_data, Diagnosis == "adenoma") %>% 
    ggplot(aes(factor(sampleType, levels = c("initial", "follow_up")), 
               log10(jitteredValues + 1.1), group = factor(EDRN))) + 
    geom_line(aes(color = factor(Disease_Free))) + geom_point(aes(color = factor(Disease_Free))) +  
    facet_wrap(~variable, labeller = as_labeller(labelFile)) + coord_cartesian(ylim = c(0, 4)) + 
    theme_bw() + ylab("Log Total Sequences") + xlab("") + ggtitle("Adenoma") + 
    scale_colour_manual(values = "blue") +  
    theme(legend.position="none", plot.margin = unit(c(1, 1, 1, 1), "lines"), 
          plot.title = element_text(size=20, face="bold"))
  
  # Cancer OTUs graph
  cancerOTUs <- filter(select_otu_data, Diagnosis == "adenocarcinoma" | Diagnosis == "N/D") %>% 
    ggplot(aes(factor(sampleType, levels = c("initial", "follow_up")), 
               log10(jitteredValues+1.1), group = factor(EDRN))) + 
    geom_line(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    geom_point(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    facet_wrap(~variable, labeller = as_labeller(labelFile)) + coord_cartesian(ylim = c(0, 4)) + 
    theme_bw() + ylab("Log Total Sequences") + xlab("") + ggtitle(paste(positiveTitle)) + 
    scale_colour_manual(name = "Cancer Free", labels = c("No", "Yes", "Unknown"), values = c("darkred", "orange", "red")) + 
    theme(legend.position=c(0.9,0.8), plot.margin = unit(c(1, 1, 1, 1), "lines"), 
          plot.title = element_text(size=20, face="bold"), legend.title = element_text(face="bold"))
  
  # Adenoma Fit graph
  adenomaFIT <- filter(fit_results, Diagnosis == "adenoma") %>% 
    ggplot(aes(factor(variable, levels = c("initial", "follow_up")), 
               log10(jitteredValues+1.1), group = factor(EDRN))) + 
    geom_line(aes(color = factor(Disease_Free))) + geom_point(aes(color = factor(Disease_Free))) + 
    theme_bw() + ylab("Log FIT Result") + xlab("") + coord_cartesian(ylim = c(0, 3.5)) + 
    scale_colour_manual(values = "blue") + theme(legend.position="none")
  
  # Cancer Fit graph
  cancerFIT <- filter(fit_results, Diagnosis == "adenocarcinoma" | Diagnosis == "N/D") %>% 
    ggplot(aes(factor(variable, levels = c("initial", "follow_up")), 
               log10(jitteredValues+1.1), group = factor(EDRN))) + 
    geom_line(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    geom_point(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    theme_bw() + ylab("Log FIT Result") + xlab("") + coord_cartesian(ylim = c(0, 3.5)) + 
    scale_colour_manual(values = c("darkred", "pink", "red")) + theme(legend.position="none")
  
  
  if(fit != FALSE){
    
    return(list(adenoma_OTUs = adenomaOTUs, cancer_OTUs = cancerOTUs, adenoma_fit = adenomaFIT, cancer_fit = cancerFIT))
    
  } else{
    
    return(list(adenoma_OTUs = adenomaOTUs, cancer_OTUs = cancerOTUs)) 
  }

}



# This function creates a visual representation of initial and follow up samples important variables for a chosen model
# These are at initial and follow up sampling.
follow_Abundance_FitGroup_Graph <- function(taxaFile, labelFile, followUpFile, initialData, followData, metadata, 
                                       jitterAmount = 0.08, n, modelOFinterest, rowOFinterest = "Disease_Free", 
                                       positiveTitle = "Cancer"){
  
  loadLibs(c("dplyr", "ggplot2", "reshape2", "wesanderson"))
  
  jitter_values <- jitter(rep(0, n), amount = 0.08)
  
  #creating the data tables to be used
  select_otu_data <- melt(
    rbind(select(initialData, starts_with(paste(modelOFinterest)), one_of(rownames(taxaFile))) %>% 
            mutate(sampleType = "initial", EDRN = metadata$EDRN, Diagnosis = metadata$Diagnosis), 
          select(followData, starts_with(paste(modelOFinterest)), one_of(rownames(taxaFile))) %>% 
            mutate(sampleType = "follow_up", EDRN = metadata$EDRN, Diagnosis = metadata$Diagnosis)), 
    id=c(paste(modelOFinterest), "sampleType", "Diagnosis", "EDRN")) %>% 
    mutate(jitteredValues = as.numeric(value) + rep(jitter_values, length(rownames(.))/n)) %>% 
    mutate(Disease_Free = rep(metadata[, rowOFinterest], length(rownames(.))/n))
  
  
  fit_results <- select(metadata, starts_with(modelOFinterest), EDRN, Diagnosis, fit_init_positive, fit_follow_positive) %>% 
    rename(initial = fit_init_positive, follow_up = fit_follow_positive) %>% 
    melt(id=c(paste(modelOFinterest), "EDRN", "Diagnosis")) %>% 
    mutate(jitteredValues = as.numeric(value) + rep(jitter_values, 2)) %>% 
    mutate(Disease_Free = rep(metadata[, rowOFinterest], 2))
  
  # Create graphs to be used
  
  # Adenoma OTUs graph
  adenomaOTUs <- filter(select_otu_data, Diagnosis == "adenoma") %>% 
    ggplot(aes(factor(sampleType, levels = c("initial", "follow_up")), 
               log10(jitteredValues + 1.1), group = factor(EDRN))) + 
    geom_line(aes(color = factor(Disease_Free))) + geom_point(aes(color = factor(Disease_Free))) +  
    facet_wrap(~variable, labeller = as_labeller(labelFile)) + coord_cartesian(ylim = c(0, 3)) + 
    theme_bw() + ylab("Log Total Sequences") + xlab("") + ggtitle("Adenoma") + 
    scale_colour_manual(values = "blue") +  
    theme(legend.position="none", plot.title = element_text(size=20, face="bold"), 
          strip.text.x = element_text(size = 5), axis.text.x = element_text(size = 6))
  
  # Cancer OTUs graph
  cancerOTUs <- filter(select_otu_data, Diagnosis == "adenocarcinoma" | Diagnosis == "N/D") %>% 
    ggplot(aes(factor(sampleType, levels = c("initial", "follow_up")), 
               log10(jitteredValues+1.1), group = factor(EDRN))) + 
    geom_line(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    geom_point(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    facet_wrap(~variable, labeller = as_labeller(labelFile)) + coord_cartesian(ylim = c(0, 3)) + 
    theme_bw() + ylab("Log Total Sequences") + xlab("") + ggtitle(paste(positiveTitle)) + 
    scale_colour_manual(values = c("darkred", "pink", "red")) + 
    theme(legend.position="none", plot.title = element_text(size=20, face="bold"), 
          strip.text.x = element_text(size = 5), axis.text.x = element_text(size = 6))
  
  # Adenoma Fit graph
  adenomaFIT <- filter(fit_results, Diagnosis == "adenoma") %>% 
    ggplot(aes(factor(variable, levels = c("initial", "follow_up")), jitteredValues, group = factor(EDRN))) + 
    geom_line(aes(color = factor(Disease_Free))) + geom_point(aes(color = factor(Disease_Free))) + 
    theme_bw() + ylab("FIT Positive Result") + xlab("") + coord_cartesian(ylim = c(0, 1.2)) + 
    scale_colour_manual(values = "blue") + scale_y_continuous(breaks = c(0, 1), labels = c("No", "Yes")) + 
    theme(legend.position="none")
  
  # Cancer Fit graph
  cancerFIT <- filter(fit_results, Diagnosis == "adenocarcinoma" | Diagnosis == "N/D") %>% 
    ggplot(aes(factor(variable, levels = c("initial", "follow_up")), jitteredValues, group = factor(EDRN))) + 
    geom_line(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    geom_point(aes(color = factor(Disease_Free, levels = c("n", "y", "unknown")))) + 
    theme_bw() + ylab("FIT Positive Result") + xlab("") + coord_cartesian(ylim = c(0, 1.2)) + 
    scale_colour_manual(values = c("darkred", "pink", "red")) + 
    scale_y_continuous(breaks = c(0, 1), labels = c("No", "Yes")) + theme(legend.position="none")
  
  return(list(adenoma_OTUs = adenomaOTUs, cancer_OTUs = cancerOTUs, adenoma_fit = adenomaFIT, cancer_fit = cancerFIT))
}









