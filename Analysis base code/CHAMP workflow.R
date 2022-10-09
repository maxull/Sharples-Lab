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


### extra functions https://github.com/YuanTian1991/ChAMP-Script

library(wateRmelon); library(methylumi);library(FDb.InfiniumMethylation.hg19);library(minfi); library(maxprobes); library(tidyverse);library(ggplot2)
library(tidyverse); library(pathview); library(gage); library(gageData); library(ChAMP); library(org.Hs.eg.db); library(AnnotationDbi);
library(reshape2)



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

myLoad2 <- champ.load(directory = my_dir,
           method="minfi",                #### method: "ChAMP" or "minfi", has to be minfi if you want to run functional normalization
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

CpG.GUI(arraytype = "EPIC")

################################################################################################################################

### create folder in which you want to store temporary files
### set working directory to that folder

setwd("/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/")

### save temporary files
saveRDS(Anno, file = "Anno.RDATA")

saveRDS(EPIC.manifest.hg19, file = "EPIC.manifest.hg19.RDATA")

saveRDS(multi.hit, file = "multi.hit.RDATA")

saveRDS(myLoad2, file = "myLoad.RDATA")

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

champ.QC()

QC.GUI()

### normalize data, chose nethod "SWAN" or functional normalization
### BMIQ is the standard normalization method
### to plot BMIQ: myNorm <- champ.norm(plotBMIQ=TRUE)
### this will save PDF of density curves

myNorm <- champ.norm(beta = myLoad$beta,
                     rgSet = myLoad$rgSet,
                     method ="FunctionalNormalization",
                     arraytype = "EPIC")                        ### does not work with this dataset

myNorm <- champ.norm(beta = myLoad$beta,
                     rgSet = myLoad$rgSet,
                     method ="BMIQ",
                     arraytype = "EPIC")                        ### bmiq worked, dont know why



champ.SVD()


################################################################################################################################

### visualize normalization





densityPlot((myNorm), sampGroups=myLoad$pd$Sample_Group,
            main="Normalized", legend=TRUE)
densityPlot(myLoad$beta, sampGroups=myLoad$pd$Sample_Group,
            main="Not-Normalized", legend=TRUE)



################################################################################################################################

### identify DMPs

champ.DMP(arraytype = "EPIC")                   ### no significant DMPs with BH adjusted p val of 0.05

myDMP<- champ.DMP(arraytype = "EPIC",
          adjust.method = "none",
          adjPVal = 0.05)               ### worked, no p val adjustement


### visualize DMPs

DMP.GUI()

setwd("/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/") ### where do you want to save DMPs

saveRDS(myDMP, myDMP.RDATA)




################################################################################################################################

### identify DMRs

myDMR <- champ.DMR(arraytype = "EPIC")          #no significant DMPs with BH adjusted p val of 0.05

myDMR <- champ.DMR(arraytype = "EPIC",
                   method = "Bumphunter",
                   adjPvalDmr = "none",
                   compare.group = c("Baseline", "7wk_Loading"),
                   cores = 4)

DMR.GUI()

setwd("/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/") ### where do you want to save DMRs

saveRDS(myDMR, myDMP.RDATA)

################################################################################################################################

### identify GSEA



### default test
myGSEA <- champ.GSEA(arraytype = "EPIC",
                     DMP = myDMP[["Baseline_to_7wk_Loading"]],
                     cores = 4)


GSEA <- as.data.frame(myGSEA$DMP)


### gOmeth method
myGSEA2 <- champ.GSEA(arraytype = "EPIC",
                      DMP = myDMP[["Baseline_to_7wk_Loading"]],
                      method = "gometh",            ### Note that gometh method would count numbers of CpGs in each genes and correct this bias.
                      adjPval = "none",
                      cores = 4)

GSEA2 <- as.data.frame(myGSEA2$DMP)

### ebayes method

############

### filtering bVals to the samples i want


bVals <- as.data.frame(myNorm)

pd <- as.data.frame(myLoad$pd)

names <- (paste(pd$Slide, pd$Array, sep = "_"))

colnames(bVals) <- names


pd %>% 
        mutate(colnames = paste(Slide, Array, sep = "_")) %>% 
        filter(Sample_Group %in% c("Baseline", "7wk_Loading")) -> pd1       ### add sample_groups you want to keep


bVals[pd1$colnames] -> f_bVals


#myGSEA3 <- champ.GSEA(beta = f_bVals,
#                      DMP = myDMP$Baseline_to_7wk_Loading,
#                      DMR = myDMR,
#                      arraytype = "EPIC",
#                      method = "ebayes",            ### fails with multiple comparison groups
#                      pheno = pd1$Sample_Group,
#                      cores = 8)

myGSEA3 <- champ.ebGSEA(beta = f_bVals,
                      arraytype = "EPIC",            ### fails with multiple comparison groups
                      pheno = pd1$Sample_Group,
                      cores = 8)

GSEA <- as.data.frame(myGSEA3[["GSEA"]][["Rank(AUC)"]])
GSEA_sig <- as.data.frame(myGSEA3[["GSEA"]][["Rank(P)"]])
################################################################################################################################
################################################################################################################################
################################################################################################################################

### map GSEA to GO and KEGG

data("go.sets.hs")
data("go.subs.hs")

data("kegg.sets.hs")
data("sigmet.idx.hs")

BiocManager::install("GIGSEA")
library(GIGSEA);library(tibble)
MSigDB.KEGG.Pathway

GSEA <- rownames_to_column(GSEA, "TERM")

df <- geneSet2Net(GSEA$TERM, geneset = myGSEA3$EnrichGene) 

# works, but still not mapped to pathway overview, or list of mosti significant pathways etc. 
# map "statistic from gtResults to kegg and GO pathways

################################################################################################################################

## cth correction based on blood methylation profiles


myRefbase <- champ.refbase(beta = f_bVals,                      ### returns b-vals adjusted for 5 main cellpopulations identified
                           arraytype = "EPIC")


## random data views


btoa <- myDMP[["Baseline_to_Acute"]]

GSEA2 <- myGSEA2[["DMP"]]

df <- myGSEA3[["GSEA"]][["Rank(AUC)"]]
df2 <-  myGSEA3[["gtResult"]]
df <- myGSEA3[["EnrichGene"]]
ebGSEA <- myGSEA3[["GSEA"]]

df1 <- as.data.frame(myGSEA3[["EnrichGene"]][["chr5q23"]])


################################################################################################################################
################################################################################################################################

BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")
library(org.Hs.eg.db)
library(AnnotationDbi)
library(gageData)
columns(org.Hs.eg.db)
data("kegg.sets.hs")

doGT(pheno.v = pd1, 
     data.m = f_bVals,
     array = "850k")
