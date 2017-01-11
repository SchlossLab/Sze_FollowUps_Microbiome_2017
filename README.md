## Differences in the Stool Microbiome Before and After Colorectal Cancer Treatment

Colorectal cancer (CRC) continues to be a worldwide health problem with early detection being used as a key component in mitigating deaths due to the disease. Previous research suggests that the bacterial microbiome can be used as a biomarker for CRC. In this study,
we used a microbiome-based model that classified individuals as having an adenoma or carcinoma to assess whether the probability of having one of these lesions changed after surgical treatment. We also characterized the change in the gut microbiota before and
after the surgical treatment. This model was tested on a 66 person group that included samples before and after treatment to allow for the assessment of how the model adjusts risk after treatment. The model used for prediction had an AUC of 0.813. For the follow up
samples our Random Forest model significantly decreased the positive probability of lesion compared to the initial samples for both adenoma (P-value = 7.95e-08) and carcinoma (P-value = 7.95e-08). Our model predicted that 36.4% of the 67-person cohort had normal
colons and 63.6% had a lesion. Some OTUs that changed the most before and after treatment included OTUs that were affiliated with members of Blautia, Clostridium_XlVa and Escherichia/Shigella. Our model suggests that treatment does significantly reduce the
probability of having a colonic lesion. Further surveillance of these individuals will enable us to determine whether models such as the one we present here can also be used to predict recurrence of colorectal cancer.

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
