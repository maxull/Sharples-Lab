################################################################
###                                                          ###
### Gene Set Enrichment analysis                             ###
###                                                          ###
################################################################

### Genome analysis and gene set enrichment analysis based on:

### Gene set enrichment analysis for genome-wide DNA methylation data
### https://doi.org/10.1101/2020.08.24.265702

### ChAMP: updated methylation analysis pipeline for Illumina BeadChips
### doi: 10.1093/bioinformatics/btx513


### set working directory to filepath that contains the data of interest
### this folder should contain the raw red and green signal .idat files and a CSV file of the targets

### save directory for later! 

mypath <- getwd()

#############################################################
### unzipping of zip files when downloading from GEO database (change file paths)

path <- "C:/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/memory_of_hypertrophy_data/GSE114763_RAW.tar"

untar(path, exdir = "C:/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/memory_of_hypertrophy_data/untared")



my_dir <- "/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/memory_of_hypertrophy_data/untared"

zip_file <- list.files(path = my_dir, pattern = "*.gz", full.names = TRUE)

ldply(.data = zip_file, .fun = gunzip)


### if CSV file is missing use code underneath (change file paths)

list <- as.data.frame(list.files(my_dir))

list <- list %>% 
        filter(row_number()<=80) %>% 
        mutate(Sentrix_ID = substr(list$`list.files(my_dir)`, start = 12, stop = 23)) %>% 
        strsplit(list$`list.files(my_dir)`, split = "_")

write.csv(list, "/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/memory_of_hypertrophy_data/names.csv", row.names = FALSE)


#############################################################
#############################################################

### load data from idat files and filter


myLoad <- champ.load(testDir, arraytype = "EPIC", method = "minfi")

myLoad <- champ.load(directory = my_dir,
           method="ChAMP",                #### method: "ChAMP" or "minfi", has to be minfi if you want to run functional normalization
           methValue="B",
           autoimpute=TRUE,
           filterDetP=TRUE,
           ProbeCutoff=0,
           SampleCutoff=0.1,
           detPcut=0.01,
           filterBeads=TRUE,
           beadCutoff=0.05,
           filterNoCG=TRUE,
           filterSNPs=TRUE,
           population=NULL,
           filterMultiHit=TRUE,
           filterXY=TRUE,
           force=FALSE,
           arraytype="EPIC")


################################################################################################################################

### create folder in which you want to store temporary files
### set working directory to that folder

setwd("/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/")

### save temporary files
saveRDS(Anno, file = "Anno.RDATA")

saveRDS(EPIC.manifest.hg19, file = "EPIC.manifest.hg19.RDATA")

saveRDS(multi.hit, file = "multi.hit.RDATA")

saveRDS(myLoad, file = "myLoad.RDATA")

saveRDS(probe.features, file = "probe.features.RDATA")


### reload temporary files (not all are neaded for downstream analysis!)

Anno <- readRDS("Anno.RDATA")

EPIC.manifest.hg19 <- readRDS("EPIC.manifest.hg19.RDATA")

multi.hit <- readRDS("multi.hit.RDATA")

myLoad <- readRDS("myLoad.RDATA")

probe.features <- readRDS("probe.features.RDATA")


### set working directory back to original


setwd(mypath)


################################################################################################################################






