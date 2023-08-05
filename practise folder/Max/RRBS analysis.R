############################################################################################
###                                                                                      ###
###  Analysis of RRBS with EpiStatProfiler                                               ###
###                                                                                      ###
############################################################################################


# https://github.com/BioinfoUninaScala/epistats/blob/main/README.md

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9636440/

library(devtools)
install_github("BioinfoUninaScala/epistats", 
               build_vignettes=FALSE, 
               repos=BiocManager::repositories(),
               dependencies=TRUE, type="source")

install.packages("vegan")
devtools::install_github("jeffkimbrel/jakR")

BiocManager::install("methylKit", force = TRUE)


library(methylKit)