# Follow Up Samples
Niel  
November 20, 2015  







![](Baxter_followUps_2016_files/figure-html/dists-1.png) 

**Figure 1. Inter- vs Intra-personal Beta Diversity**  
The first thing I tested was whether the microbiomes of patients were more similar to themselves over time than to others. Surprisingly, there was no difference between intra- and interindividual thetaYC distances. Distances within cancer patients seems to be slightly less than within adenomas, but the difference isn't significant.  There's no difference when splitting interindividual distances by diagnosis. There's also no difference when comparing distances between initial samples to distances between followups.


![](Baxter_followUps_2016_files/figure-html/combined_model-1.png) 

**Figure 2A. Combined Model (FIT + OTUs) applied to initial and follow up samples.**  
Next, I generated a random forest model for distinguishing normal from lesion using all of the initial samples that didn't have follow ups (n=422). The optimal model used 70 OTUs and FIT. I then applied the model to the other initial samples and their follow ups. The scatter plot below shows how the model classified each individual's initial and follow up samples (n=67). The model detected 26 of 26 cancers and 25 of 41 adenomas in the initial samples. Of the 51 individuals who tested positive with the initial sample, 36 remained positive with the follow up sample. 10.4% of individuls were negative for both the initial and followups (all adenomas). Surprisingly, 13.4% of individuals were negative at the initial sampling, but became positive at the follow up.

![](Baxter_followUps_2016_files/figure-html/fit_dif-1.png) 

**Figure 2B.** This next figure shows the change in FIT result from the initial to the follow up samples. For the adenomas, it's fairly random. Some go from positive to negative, but just as many go form negative to positive. For the cancers, however, a large portion of them go from high FIT results to zero, while none of them go from zero to higher results. As a result, there is a significant reduction in FIT results for cancers, but not for adenomas. 

![](Baxter_followUps_2016_files/figure-html/comb_dif-1.png) 

**Figure 2C.** This next figure shows how the model results changed over time. Based on the previous figure, you'd expect outcome of the model to change more for the cancer samples than for the adenomas. Indeed, there is a significant reduction in the probability of lesion for the cancers but not the adenomas. This is most likely due to the reduction in FIT results indiviuduals that had cancers.


![](Baxter_followUps_2016_files/figure-html/otu_model-1.png) 

**Figure 3A. Otucome of Microbiome model (no FIT) on the initial and followup samples.** I did the same thing as figure 2, but left out the FIT results from the model. In this case the optimal used only 13 OTUs, yet it detected lesions nearly as well as the combined model. It detected detected 18 of 26 cancers (worse than the combined model) and 30 of 41 adenomas (better than the combined model) in the initial samples. Of the 48 initial samples that tested positive, 21 remained positive. In summary, compared to the combined model, the microbiome only model was better at detected adenomas, worse at detecting cancers, and individuals were more likely to appear to "normal" at the follow up.

![](Baxter_followUps_2016_files/figure-html/dif_otu-1.png) 

**Figure 3B.** Unlike the combined model, there was significant reducation in probability of lesion for both patients with adenomas and those with cancer when using the microbiome only model.


![](Baxter_followUps_2016_files/figure-html/time-1.png) 

**Figure 4.** I thought perhaps a longer time between sampling would increase the likelihood that an individual would appear "normal". However, there was no correlation between the change in model outcome and the time between samples. 

![](Baxter_followUps_2016_files/figure-html/unnamed-chunk-1-1.png) 


**Figure 6. The usual suspects.** Finally I wanted to see whether the OTUs most commonly associated with cancer were less abundant in the follow up samples. The figure above shows 6 of the 14 OTUs used in the microbiota only model. Most of them should look very familiar. You can see that for most of them there is a significant reduction in relative abundance from the initial sample to the follow up (p-values based on one-sided paired wilcoxon test). Interestngly, for Fusobacterium there is not a significant reduction unless you only compare cancer samples (p=0.046). 

**Wild Speculation:**
To me, these finding are consistent with the idea that these oral pathogens are enriched on tumors, because they like the environment. When the tumors are removed their abundance goes down. They may accelerate tumorigenesis when they colonize the tumor, but first it's the tumor and/or inflammation that enables them to colonize in the first place. It's also possible, as we've  hypothesized in lab meetings, that Fuso binds first during the earliest stages (or even prior to) tumorigenesis.  Its presence then enables the other pathogens to adhere and exacerbate disease (similar to the progression of periodontitis). That would explain why people often find Fuso associated with adenomas, while these other oral pathogens are more indicative of cancers.





