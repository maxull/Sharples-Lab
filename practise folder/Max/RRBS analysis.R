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
BiocManager::install("genomation")

library(methylKit)
library(tidyverse)
library(ggplot2)
library(GEOquery)
library(genomation)
library(GenomicRanges)


setwd("/Users/maxul/Documents/Skole/Lab/CrosSys/")


#gunzip files

files <- list.files("./raw_data/", full.names = TRUE)

for (i in 1:length(files)) {
        gunzip(filename = files[i], remove = TRUE)
        print(i)
}


fq_files <- list.files("./raw_data/", pattern = "*.fq",full.names = TRUE) %>% head(10)

fq_names <- list.files("./raw_data/", pattern = "*.fq",full.names = FALSE) %>% head(10)






# test 
setwd("/Users/maxul/Documents/Skole/Lab/CrosSys/test/")

system("bismark --bowtie2 /Users/maxul/Documents/Skole/Lab/RRBS_analysis/GRCh376 -1 220048-II-01-011-PRE-ASAT_S1_L001_R1_001_val_1.fq_trimmed.fq -2 220048-II-01-011-PRE-ASAT_S1_L001_R2_001_val_2.fq_trimmed.fq")


# Extract methylation data
system("bismark_methylation_extractor --paired-end --bedGraph --comprehensive sample1_trimmed_R1_bismark_bt2_pe.bam")
system("bismark --version")
