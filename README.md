## The Fecal Microbiome Before and After Treatment for Colorectal Adenoma or Carcinoma

**Background.** Colorectal cancer (CRC) is a worldwide health problem. Despite growing evidence that members of the gut microbiota can drive tumorigenesis, little is known about what happens to the microbiota after treatment for an adenoma or carcinoma. This study tested the hypothesis that treatment for adenoma or carcinoma alters the abundance of bacterial populations associated with disease to those associated with a normal colon. We
tested this hypothesis by sequencing the 16S rRNA genes in the feces of 67 individuals
before and after treatment for adenoma (N=22), advanced adenoma (N=19), and carcinoma
(N=26).  
**Results.** There were large changes to the bacterial communities associated with treatment
across the three groups. The communities from patients with carcinomas changed
significantly more than those with adenoma following treatment (P-value=5.4e-05); however,
there was not a significant difference between those with advanced adenoma and those
with adenoma or carcinoma (P-value > 0.05). Although treatment brought about large
intrapersonal changes, the change in the abundance of individual OTUs to treatment
was not consistent within diagnosis groups (P-value > 0.05). Because the distribution
of OTUs across patients and diagnosis groups was patchy, we used the Random Forest
machine learning algorithm to identify groups of OTUs that allowed us to successfully
distinguish between pre- and post-treatment samples for each of the diagnosis groups.
However, across the three models, there was little overlap between the OTUs that were
indicative of treatment. Next, we used a larger cohort that contained individuals with normal colons and those with adenomas, advanced adenomas, and carcinomas to determine
whether individuals who underwent a treatment were more likely to have OTUs associated
with normal colons. We again built Random Forest models and measured the change
in the positive probability of having one of the three diagnoses. Although we could
clearly differentiate pre- and post-treatment communities from the three diagnosis groups,
only those patients that initially had carcinomas experienced a significant decrease in
positive probability of having a carcinoma (P-value < 0.05). Finally, we tested whether
the type of treatment impacted the microbiota of those diagnosed with carcinomas and
were unable to detect any significant differences in characteristics of these communities
between individuals treated with surgery alone and those treated with chemotherapy or
chemotherapy and radiation (P-value > 0.05).   
**Conclusions.** Further exploration of the relationship between diagnosis, treatment, and
the impact on the microbiota will yield improvements in disease management.

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
