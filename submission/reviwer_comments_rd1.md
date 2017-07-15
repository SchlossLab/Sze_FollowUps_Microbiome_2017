# Response to Reviewers

## Reviewer #1

**Main points**

* Post-treatment samples were collected 188-546 days after lesion treatment, I wonder if that is not too long to evaluate the effect of the treatment. Many other factors may interfere such as an adaptation in their nutritional habits or changes in their medical chronic treatments, for instance.



* I think "Treatment" is not accurate enough, a treatment is a "fixed" factor in a design, almost any intervention can be consider as a treatment. It can be nutritional, pharmacological, surgical, .. I would specify it more, maybe something like "anti-cancer treatment". 



* The consumption of antibiotics is not reported in the study. Cancer patients received antibiotics for several reasons, pending the clinical situation. It is hard to decipher what is due to the antibiotic treatment versus the surgical intervention and the chemotherapy. At least, antibiotic regimen should be reported. Statement in lines 220-222 is not enough.  Antibiotics can have long-term effects on the gut microbiota composition by completely eradicating sensitive microorganisms.



* If the point is actually to search for microbial biomarkers of recurrence and survival (as stated in the conclusion of the abstract and the manuscript), it would be relevant to provide recurrence and survival data, at least for the 26 patients with carcinoma.



* Line 161-163. "We focused on the patients with carcinoma and pooled those patients that received chemotherapy with those that received chemotherapy and radiation to improve our statistical power." Radiation influences the gut microbiota by itself (see for instance Cui et al EMBO Mol Med 2017 and Gerassy-Vainberg et al, Gut 2017), so I don't think that the pooling is scientifically justified.

**Additional points**

* Fig. 1A shows the diversity pre-post for the 3 lesion types. It might help to have the same data from a control cohort (that provided samples at the same time intervals) to evaluate which portion of the change is actually do to the temporal variation in the gut microbiota and which portion relates to treatment. Here, all we can say is that carcinoma treatment and associated events induces greater changes than adenoma and advanced adenoma ones.


* Lines 226-229. "The most exciting future direction from the current study is the possibility that markers within the microbiota could be used to evaluate the effect of treatment and predict recurrence for those diagnosed with carcinoma. If such an approach is effective, it might be possible to target the microbiota as part of adjuvant therapy." A biomarker does not always play a key role in the disease.


* The study cohort (67 individuals, before and post treatment) is clearly presented. However, the presentation of the training cohort does not seem clear to me. From Table 2, we don't know if all samples are from post or pre-treatment. This cohort is not presented in the mat and meth section. From line 124, I guess they are all post-treatment? When were they collected?


* No information is provided on the type of treatment (which chemo for instance).



* Sampling described in Baxter et al Genome Med 2016 seems to suggest that samples were stored on ice for at least 24h before freezing. We cannot exclude that this sampling condition may have impact the gut microbiota composition. I would suggest mentioning at least the sampling condition in the main text of the manuscript.


* It might help to better describe Figure 2. What does the light grey line stand for?



* 	From the results section, it appears that no definitive statement can be made regarding the effect of a specific treatment on the gut microbiota due to the low number of samples, which I think is correct; however such statement are made at the end of the introduction (lines 84-85) and in the abstract (lines 29-33).



* Lines 69-71: "If the microbial community drives tumorigenesis then one would hypothesize that treatment to remove a lesion would affect the microbiota and risk of recurrence." I think the sentence should be rephrased, as it is known that treatment to remove lesion affects risk of recurrence.


* Line 83 : typo


* Line 142 : "showing are marked"








## Reviewer #2

* The statement in the abstract "There were large changes to the bacterial communities associated with treatment across the three groups" seems inconsistent with the statement "The change in the abundance of individual OTUs to treatment was not consistent within diagnosis groups (P-value > 0.05)."  If I understand the authors claims correctly, the changes within treatment group could be detected with machine learning approaches (random forest) but not inferential statistics.  Moreover, Figs. 1B-1C seem to suggest that the PCOA of pre-and post- treatment are superimposable.  Given this, it seems a better summary sentence would be that there were larger changes associated with colorectal cancer but much smaller, and more difficult to detect, changes associated with adenomas and advanced adenomas.



* "Only patients who had carcinomas experienced a significant decrease in positive probability of having a lesion after treatment (P-value < 0.05), indicating that the microbial milieu of the colon more closely resembled that of a normal colon."  This sentence in the abstract is confusing.  Do the authors mean that only patients with carcinomas had a significant decrease in the similarity of the microbial community to other patients in a previous cohort that also had colorectal cancer?




* I find it a little odd that for the adenomas and advanced adenomas groups, there was no difference between pre-treatment and treatment by change in the beta diversity metric (lines 94-96), PERMANOVA (line 102), number of observed OTUs, Shannon evenness, and Shannon diversity (line 104), and OTUs that were significantly different in the pre and post-treatment  groups (line 106).  It is of concern that differences in adenomas and advanced adenomas were only detectable by machine learning (Random Forest) techniques.  Could the authors use permutation tests (randomly permuting the pre- and post- treatment labels) and rerun their Random Forest pipeline over each permutation to ensure that the results for adenomas and advanced adenomas are not due to over-fitting?   In particular, I am concerned that the results could reflect the optimization of the mtry parameter described in lines 277-280.  Any permutation scheme would need to randomly assign treatment and pre-treatment
parameters and then run the entire optimization scheme over the 80/20 splits and through the rest of the pipeline.  It has been demonstrated in the literature that optimization over cross-training can still produce substantial overfitting ( e.g. https://doi.org/10.1093/bioinformatics/btr591; "Optimized application of penalized regression methods to diverse genomic data") and given the lack of any other observable differences between adenomas and post-adenomas, it seems likely that that kind of overfitting has occurred here.




* I would also think that language in the abstract and conclusion could better emphasize how small the effect was for adenomas and advanced adenomas.  Given that "treatment" for adenomas is presumably just removing the adenomas during colonoscopy, it is reasonable to expect a difference in the microbial community?  Can the authors more explicitly state what treatment was for the adenomas and pre-adenomas groups?



* Likewise, in lines 130-133, could the authors perform a permutation test to give p-values for the result of reproducibility with previous cohorts.  And in lines 149-152, could a permutation test be used to provide a p-value for statement "The positive probability for the pre and post-treatment samples from patients diagnosed with carcinoma significantly decreased with treatment, suggesting a shift toward a normal microbiota for most individuals".

