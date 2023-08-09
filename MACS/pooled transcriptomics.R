############################################################################################
###                                                                                      ###
###  Pooled transcriptomics for MACS                                                     ###
###                                                                                      ###
############################################################################################

# purpose is to get gene expression changes from RT training intervent that performed transcriptomics


# method: 

#       1. identify chronic RT protocols (more than 4 weeks)

#       2. download raw transcriptome from arrays and RNAseq from GEO etc.

#       3. QC and normalize all the same

#       4. correct for study/batch effects

#       5. differential gene expression

#       step 2: overlap with MACS epigenommic changes

#       step 3: identify genes of interest, and confim differential gene expression by qPCR


library(GEOquery)
library(tidyverse)
library(tweeDEseqCountData)
library(limma)
library(Biobase)
library(biomaRt)
library(tweeDEseq)

# set working directory to a local, and not github connected folder

setwd("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/")

# get raw data and pheno data from GEO


# study 1 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28422
# Raue et al. (2012)

GSE28422 <- getGEO(GEO = "GSE28422", GSEMatrix = TRUE)

metadata_GSE28422 <- pData(phenoData(GSE28422[[1]]))

# have to manually download raw data

# untar the file and save in directory

untar(tarfile = "/Users/maxul/Downloads/GSE28422_RAW.tar",
        exdir = "/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/GSE28422")

# unpack .gz files

# get list of files

gz_files <- list.files("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/GSE28422", full.names = TRUE)


for (i in 1:length(gz_files)) {
        gunzip(filename = gz_files[i], remove = TRUE)
        print(i)
}

# list .CEL files

cel_files <- list.files("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/GSE28422", 
           full.names = TRUE, 
           pattern = "*.CEL")

library(affy)

GSE28422_raw <- read.affybatch(cel_files)

?"AffyBatch"

df <- rma(GSE28422_raw)

GSE28422_exprs <- exprs(df)

# pm = perfect match

pm <- as.data.frame(probes(GSE28422_raw, which="pm"))

names(pm)[1] <- "ProbesetID"



data <- expresso(GSE28422_raw, bg.correct = FALSE, 
                 normalize = FALSE, 
                 pmcorrect.method = "pmonly", 
                 summary.method = "avgdiff")


exprs(data)


# annotate probes

mart <- useMart("ENSEMBL_MART_ENSEMBL")

mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
        mart=mart,
        attributes=c(
                "affy_hg_u133_plus_2",
                "ensembl_gene_id",
                "gene_biotype",
                "external_gene_name"),
        filter = "affy_hg_u133_plus_2",
        values = rownames(GSE28422_exprs))


unique(annotLookup$external_gene_name) %>% 
        length()





# study 2 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47881
# Phillips et al. (2013)







