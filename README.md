## Normalization of the microbiota in patients after treatment for colonic lesions

**Background.** Colorectal cancer is a worldwide health problem. Despite growing evidence that members of the gut microbiota can drive tumorigenesis, little is known about what happens to it after treatment for an adenoma or carcinoma. This study tested the hypothesis that treatment for adenoma or carcinoma alters the abundance of bacterial populations associated with disease to those associated with a normal colon. We tested this hypothesis by sequencing the 16S rRNA genes in the feces of 67 individuals before and after treatment
for adenoma (N = 22), advanced adenoma (N = 19), and carcinoma (N = 26).

**Results.** There were small changes to the bacterial community associated with adenoma or advanced adenoma and large changes associated with carcinoma. The communities from patients with carcinomas changed significantly more than those with adenoma following treatment (P-value < 0.001). Although treatment was associated with intrapersonal changes, the change in the abundance of individual OTUs to treatment was not consistent within diagnosis groups (P-value > 0.05). Because the distribution of OTUs across patients and diagnosis groups was irregular, we used the Random Forest machine learning algorithm to try and identify groups of OTUs that could distinguish between pre and post-treatment samples for each of the diagnosis groups. Although the three models could differentiate between the pre and post-treatment samples, there was little overlap between the OTUs that were indicative of treatment. Next, we used a larger cohort that contained individuals with normal colons and those with adenomas, advanced adenomas, and carcinomas to determine whether individuals who underwent treatment were more likely to have OTUs associated with normal colons. We again built Random Forest models and measured the change in the positive probability of having one of the three diagnoses to assess whether these changes were similar to their original diagnosis or not. Patients who had carcinomas changed towards a microbial milieu that resembles the normal colon (P-value < 0.05). Finally, we were unable to detect any significant differences in the microbiota between individuals treated with surgery alone and those treated with chemotherapy or chemotherapy and radiation (P-value > 0.05).  

**Conclusions.** Although it was difficult to identify significant differences between pre and post treatment samples for adenoma and advanced adenoma, this was not the case for carcinomas. Not only were there large changes in pre versus post treatment samples for those with carcinoma but also these changes were towards a more normal microbiota composition.



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
