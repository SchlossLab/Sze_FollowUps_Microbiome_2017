## The Fecal Microbiome Before and After Treatment for Colorectal Adenoma or Carcinoma

**Background:** Colorectal cancer (CRC) is a worldwide health problem and research suggests a correlation between the fecal bacterial microbiome and CRC.  Despite this, very little is known about what happens to the microbiome after treatment for an adenoma or carcinoma.  This study tested the hypothesis that treatment for adenoma or carcinoma results in changes towards a normal bacterial community.  Specifically, we tried to identify components within the community that were different before and after treatment of adenoma, advanced adenoma, and carcinoma.  
**Results:** There was a larger change in the bacterial community in response to treatment for carcinoma versus adenoma (P-value < 0.05) but not carcinoma versus advanced adenoma (P-value > 0.05). There was a trend for increasingly less community similarity, between 12 samples pre- and post-treatment from adenoma to advanced adenoma to carcinoma. Despite this, no difference was found in the relative abundance of any specific OTU before and after treatment for adenoma, advanced adenoma, or carcinoma groups (P-value > 0.05). Using Random Forest models to assess whether changes in post-treatment samples were towards a normal community, only those with carcinoma had a significant decrease in positive probability (P-value < 0.05); providing further evidence that treatment has the greatest effect in those with carcinoma. The adenoma model used a total of 62 OTUs, the SRN model used a total of 61 OTUs, and the carcinoma model used a total of 59 OTUs. A total of 26 OTUs were common to all three models with many classifying to commensal bacteria (e.g. Lachnospiraceae, Bacteroides, Anaerostipes, Blautia, and Dorea). Both chemotherapy and radiation did not provide any additional changes to the
bacterial community in those treated for carcinoma (Pvalue > 0.05).  
**Conclusions:** Our data partially supports the hypothesis that the microbiome changes after treatment towards a normal community. Individuals with carcinoma had more drastic differences to the overall community than those with adenoma. Common OTUs to all
models were overwhelmingly from commensal bacteria, suggesting that these bacteria may be important in initial polyp formation, development of advanced adenoma, and transition to carcinoma.

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
