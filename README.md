## Differences in the Stool Microbiome Before and After Colorectal Cancer Treatment

**Background:** Colorectal cancer (CRC) continues to be a worldwide health problem with early detection being used as a key component in mitigating deaths due to the disease.  Previous research suggests a link between stool bacterial microbiome and CRC. The
overall objective was to investigate the changes in the bacterial microbiome after surgery in patients with lesion (i.e. adenoma or carcinoma). Specifically, we wanted to identify what within the community was different within those undergoing surgical removal of lesion. We also investigated the use of the bacterial microbiome and Fecal Immunoglobulin Test (FIT) to build models to either classify individuals as having a lesion or whether, based on the bacterial microbiome, the sample could be classified correctly as before or after surgery.   
**Results:** Adenoma individual’s bacterial microbiome were more similar to their pre-surgery sample then those with carcinoma (P-value = 0.00198) and this was also reflected in FIT as well (P-value = 2.15e-05). There was no significant difference in any indivdiual OTU between samples before and after surgery (P-value > 0.125). A model with a total of 37 variables was able to classify lesion with an AUC range of 0.847 to 0.791 while the model to classify samples as before and after had 33 with an AUC range of 0.79 to 0.651 for 100 20 repeated 10-fold cross-validated runs. Both models had a significant decrease in the positive probability of a lesion between individual’s before versus after surgery samples (P-value = 1.91e-11 and 6.72e-12). In total there were 14 OTUs that were common to both models and were mostly commensals with largest representation from OTUs belonging to Bacteroides, Blautia, Streptococcus, and Clostridiales.   
**Conclusions:** Our data suggests that treatment not only significantly reduces the probability of having a colonic lesion but also causes detectable changes in the bacterial microbiome. Further surveillance of these individuals will enable us to determine whether models such as the one we present here can also be used to predict recurrence of colorectal cancer.


### Overview
	project
	|- README # the top level description of content (this doc) - 
	|CONTRIBUTING # instructions for how to contribute to your 
	|project - LICENSE # the license for this project
	|
	|- submission/
	| |- study.Rmd # executable Rmarkdown for this study, if 
	| |applicable - study.md # Markdown (GitHub) version of the 
	| |*.Rmd file - study.tex # TeX version of *.Rmd file - 
	| |study.pdf # PDF version of *.Rmd file - header.tex # LaTeX 
	| |header file to format pdf version of manuscript - 
	| |references.bib # BibTeX formatted references - XXXX.csl # csl 
	| |file to format references for journal XXX
	|
	|- data # raw and primary data, are not changed once created
	| |- references/ # reference files to be used in analysis - raw/ 
	| |# raw data, will not be altered
	| |- mothur/ # mothur processed data
	| +- process/ # cleaned data, will not be altered once created;
	|                 # will be committed to repo
	|
	|- code/ # any programmatic code
	|
	|- results # all output from workflows and analyses
	| |- tables/ # text version of tables to be rendered with kable 
	| |in R - figures/ # graphs, likely designated for manuscript 
	| |figures
	| +- pictures/ # diagrams, images, and other non-graph graphics
	|
	|- exploratory/ # exploratory data analysis for study
	| |- notebook/ # preliminary analyses
	| +- scratch/ # temporary files that can be safely deleted or 
	| lost
	|
	+- Makefile # executable Makefile for this study, if applicable
### How to regenerate this repository
#### Dependencies and locations  
* Gnu Make (v3.81) should be located in the user's PATH  
* mothur (v1.37.5) should be located in the user's PATH
	* Note v1.37.0 causes all sorts of headaches  	
* R (v. 3.3.0) should be located in the user's PATH  
* LaTex should be in the user's PATH
* Pandoc should in the user's PATH
* PBS job scheduler for server (or something similar)
* ghostscript libtiff-tools (to convert tiff to pdf)
	* Can check on linux command line by typing `which tiff2pdf`

#### Running analysis  
```git clone https://github.com/SchlossLab/Baxter_followUps_2016```  
```make write.paper```
