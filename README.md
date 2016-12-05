## Differences in the Stool Microbiome Before and After Colorectal Cancer Treatment

Colorectal cancer (CRC) continues to be a worldwide health problem with early detection being used as a key component in mitigating deaths due to the disease. Previous research suggests that the bacterial microbiome can be used as a biomarker for colorectal cancer (CRC). These reports have mostly focused on investigating how the bacterial microbiome is used at a single point in time to predict disease. In this study, we assessed whether a model built with bacterial microbiome data could accurately predict adenoma
or carcinoma and adjust positive probability for lesion (adenoma or carcinoma) after surgical treatment. This model was tested on a 66 person group that included samples before and after treatment to allow for the assessment of how the model adjusts risk after treatment. The model chosen had a 10-fold cross validated AUC of 0.819 and a test set AUC of 0.815. : For the follow up samples our Random Forest model significantly decreased the positive probability of lesion compared to the initial samples for both adenoma (P-value = 8.48e-07) and carcinoma (P-value = 5.96e-08). The top 5 most important operational taxonomic units (OTUs) for prediction of lesion classified to
Lachnospiraceae (Otu000013), Escherichia/Shigella (Otu000018), Ruminococcaceae (Otu000020), Ruminococcus (Otu000017), and Porphyromonas (Otu000153). In the test set, the prediction of lesion for the initial samples had a sensitivity of 98.5% (65/66) while follow up it was 100% (1 / 1). Overall, our model accurately predicted those with an adenoma or carcinoma and following treatment it also decreased the positive probability of having a lesion.


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

#### Running analysis  
```git clone https://github.com/SchlossLab/Baxter_followUps_2016```  
```make write.paper```
