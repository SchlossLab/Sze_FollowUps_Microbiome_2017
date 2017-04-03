### List of commands to render final pdf
### renders both main text and supplemental

# Load needed functions
source('code/functions.R')

# Load needed library
loadLibs(c("knitr", "rmarkdown"))

# Render the final pdfs
render('submission/manuscript_outline_20161024.Rmd', clean = FALSE)

render('submission/supplemental_outline_20161024.Rmd', clean=FALSE)

render('results/tables/Table1.Rmd', clean=FALSE)

render('results/tables/Table2.Rmd', clean=FALSE)

