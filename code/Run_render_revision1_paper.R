### List of commands to render Revision 1 pdf
### renders both main text and supplemental

# Load needed functions
source('code/functions.R')

# Load needed library
loadLibs(c("knitr", "rmarkdown"))

# Render the final pdfs
render('submission/manuscript_R1.Rmd', clean = FALSE)