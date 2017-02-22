## Differences in the Stool Microbiome Before and After Colorectal Cancer Treatment

**Background:** Colorectal cancer (CRC) continues to be a worldwide health problem with previous research suggesting that a link may exist between the fecal bacterial microbiome and CRC. The overall objective of our study was to test the hypothesis that changes in the bacterial microbiome occur after lesion (i.e. adenoma or carcinoma) removal. Specifically, we wanted to identify what within the community was different before and after removal of said lesion.  
**Results:** The bacterial microbiome changed more in response to lesion removal in carcinoma cases compared to adenoma cases (P-value < 0.05). There was no difference for either the adenoma or carcinoma group in the relative abundance of any OTU between pre and post lesion removal (P-value > 0.05). A model built to classify lesion had an AUC range of 0.811 - 0.866 while a model built to classify initial samples had an AUC range of 0.641 - 0.805. Follow up samples had a decrease in the positive probability of lesion (P-value < 0.05) suggesting a movement towards a more normal bacterial community. The initial sample model had a decrease in positive probability for the follow up samples to be an initial sample (P-value < 0.05). The lesion model used a total of 53 variables while the initial sample model used a total of 70 variables. A total of 23 OTUs were common to both models with the majority of these classifying to commensal bacteria (e.g. Bacteroides, Clostridiales, Blautia, and Ruminococcaceae).   
**Conclusions:** Our data supports the hypothesis that there are changes in the bacterial microbiome following colorectal cancer lesion removal. Individuals with carcinoma have more drastic differences to the overall community then those with adenoma. Changes to  commensal bacteria were some of the most important variables for model classification, suggesting that these bacteria may be central to initial polyp formation and transition to carcinoma.


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
