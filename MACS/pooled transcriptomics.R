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

BiocManager::install("GEOquery")



library(GEOquery)
library(tidyverse)
library(tweeDEseqCountData)
library(limma)
library(Biobase)
library(biomaRt)
library(tweeDEseq)
library(affy)
library(ggplot2)
library(ggrepel)

# set working directory to a local, and not github connected folder

setwd("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/")


# get raw data and pheno data from GEO




#


# annotate affymetrix probes

mart <- useMart("ENSEMBL_MART_ENSEMBL")

mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
        mart=mart,
        attributes=c(
                "affy_hg_u133_plus_2",
                "ensembl_gene_id",
                "gene_biotype",
                "external_gene_name",
                "hgnc_symbol",
                "entrezgene"),
        filter = "affy_hg_u133_plus_2",
        values = rownames(GSE47881_exprs))

listAttributes(mart)

saveRDS(annotLookup, file = "./annotLookup.RDATA")

annotLookup <- readRDS("./annotLookup.RDATA")

# Load the package
library(org.Hs.eg.db)
data(org.Hs.eg.db)

annotLookup <- annotLookup %>% 
        mutate(ENTREZID = mapIds(org.Hs.eg.db, keys = annotLookup$ensembl_gene_id, keytype="ENSEMBL", column = "ENTREZID"))
#################################################################################################


# study 1 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28422
# Raue et al. (2012)

#################################################################################################



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


# read data GSE28422

GSE28422_raw <- read.affybatch(cel_files)


GSE28422_eset <- expresso(GSE28422_raw,
                          bg.correct = TRUE, bgcorrect.method = "rma",               # RMA is background correct method correcting probe intensities by a global method
                          normalize.method = "quantiles",
                          pmcorrect.method = "pmonly",
                          summary.method = "avgdiff")                               # computed average

exprs(GSE28422_eset)

plotDensities(exprs(GSE28422_eset), legend = FALSE)


# plot pca


pca.out <- prcomp(t(log2(exprs(GSE28422_eset))), scale. = FALSE)

# check proportion of variability explained by pca 1-50

# Extract the proportion of variance explained by each principal component
var_explained <- pca.out$sdev^2 / sum(pca.out$sdev^2)

# Create a data frame for plotting
df <- data.frame(Component = 1:length(var_explained), Variance = var_explained)

# Keep only the first 50 principal components
df <- df[1:50,]

# Plot

ggplot(df, aes(x = Component, y = Variance)) +
        geom_bar(stat = "identity") +
        scale_x_continuous(breaks = 1:50) +
        labs(x = "Principal Component", y = "Proportion of Variance Explained",
             title = "Variance Explained by Principal Components")


plot(pca.out$x[,1:2])


# check for outliers

# Compute Mahalanobis distances for the first two principal components
distances <- mahalanobis(pca.out$x[,1:2], colMeans(pca.out$x[,1:2]), cov(pca.out$x[,1:2]))

# Identify outliers as samples with a Mahalanobis distance greater than a certain threshold
# Here, I'm using the 97.5 percentile of the Chi-square distribution with 3 degrees of freedom as the threshold
outliers <- which(distances > qchisq(0.975, df = 3))

# Print the row names of the outliers
print(colnames(exprs(GSE28422_eset))[outliers])

# Plot the first two principal components, highlighting the outliers

ggplot(data.frame(pca.out$x), aes(x = PC1, y = PC2)) +
        geom_point() +
        geom_text_repel(data = data.frame(pca.out$x)[outliers,], aes(label = rownames(data.frame(pca.out$x))[outliers])) +
        labs(x = "PC1", y = "PC2", title = "First Two Principal Components with Outliers Highlighted")


# remove outliers

GSE28422_exprs <- log2(exprs(GSE28422_eset))

GSE28422_exprs.f <- GSE28422_exprs[,!(colnames(GSE28422_exprs) %in% names(outliers))]

# remove .cel from name

colnames(GSE28422_exprs.f) <- gsub(pattern = "*.CEL", replacement = "", colnames(GSE28422_exprs.f))


# identify pre-post samples

metadata_GSE28422 %>% 
        filter(characteristics_ch1 == "age: Young", 
               characteristics_ch1.1 == "gender: Male",
               characteristics_ch1.3 == "time point: Basal") %>% 
        mutate(timepoint = ifelse(`training state:ch1`== "Untrained", "Pre", "Post")) %>% 
        dplyr::select(timepoint, description.1) %>% 
        mutate(FP = sapply(strsplit(description.1, split = "_"),"[", 6)) %>% 
        dplyr::select(timepoint, FP) %>% 
        mutate(FP = paste(FP, "GSE28422", sep = ""),
               Batch = "GSE28422") ->  df_GSE28422
        

FP <- factor(df$FP)
timepoint <- factor(df$timepoint, levels = c("Pre", "Post"))


design <- model.matrix(~FP+timepoint)


colnames(GSE28422_exprs) <- gsub(pattern = "*.CEL", replacement = "", colnames(GSE28422_exprs))

GSE28422_exprs %>% 
        as.data.frame() %>% 
        dplyr::select(rownames(df)) %>% 
        lmFit(., design) -> fit


fit <- eBayes(fit)


GSE28422_res <- topTable(fit, coef = "timepointPost",number = Inf) %>% 
        filter(P.Value < 0.05) %>% 
        rownames_to_column(var = "affy_hg_u133_plus_2") %>% 
        merge(.,annotLookup, by = "affy_hg_u133_plus_2") %>% 
        arrange(-abs(logFC)) 

GSE28422_res.f <- GSE28422_res[!duplicated(GSE28422_res$logFC),]


GSE28422_res %>% 
        distinct(affy_hg_u133_plus_2, logFC, .keep_all = TRUE) 





saveRDS(GSE28422_res.f, file = "./GSE28422_res.f.RDATA")


#################################################################################################

# study 2 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47881
# Phillips et al. (2013)

#################################################################################################


# get metadata

GSE47881 <- getGEO(GEO = "GSE47881", GSEMatrix = TRUE)

metadata_GSE47881 <- pData(phenoData(GSE47881[[1]]))


# have to manually download raw data

# untar the file and save in directory

untar(tarfile = "/Users/maxul/Downloads/GSE47881_RAW.tar",
      exdir = "/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/GSE47881")

# unpack .gz files

# get list of files

gz_files <- list.files("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/GSE47881", full.names = TRUE)


for (i in 1:length(gz_files)) {
        gunzip(filename = gz_files[i], remove = TRUE)
        print(i)
}

# list .CEL files

cel_files <- list.files("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/GSE47881", 
                        full.names = TRUE, 
                        pattern = "*.CEL")


# read data GSE28422

GSE47881_raw <- read.affybatch(cel_files)


GSE47881_eset <- expresso(GSE47881_raw,
                          bg.correct = TRUE, bgcorrect.method = "rma",               # RMA is background correct method correcting probe intensities by a global method
                          normalize.method = "quantiles",
                          pmcorrect.method = "pmonly",
                          summary.method = "avgdiff")                               # computed average


# create expression dataframe

exprs(GSE47881_eset)

plotDensities(log2(exprs(GSE47881_eset)), legend = FALSE)


# plot pca


pca.out <- prcomp(t(log2(exprs(GSE47881_eset))), scale. = FALSE)

# check proportion of variability explained by pca 1-50

# Extract the proportion of variance explained by each principal component
var_explained <- pca.out$sdev^2 / sum(pca.out$sdev^2)

# Create a data frame for plotting
df <- data.frame(Component = 1:length(var_explained), Variance = var_explained)

# Keep only the first 50 principal components
df <- df[1:50,]

# Plot

ggplot(df, aes(x = Component, y = Variance)) +
        geom_bar(stat = "identity") +
        scale_x_continuous(breaks = 1:50) +
        labs(x = "Principal Component", y = "Proportion of Variance Explained",
             title = "Variance Explained by Principal Components")


plot(pca.out$x[,1:2])


# check for outliers

# Compute Mahalanobis distances for the first two principal components
distances <- mahalanobis(pca.out$x[,1:2], colMeans(pca.out$x[,1:2]), cov(pca.out$x[,1:2]))

# Identify outliers as samples with a Mahalanobis distance greater than a certain threshold
# Here, I'm using the 97.5 percentile of the Chi-square distribution with 3 degrees of freedom as the threshold
outliers <- which(distances > qchisq(0.975, df = 3))

# Print the row names of the outliers
print(colnames(exprs(GSE47881_eset))[outliers])

# Plot the first two principal components, highlighting the outliers

ggplot(data.frame(pca.out$x), aes(x = PC1, y = PC2)) +
        geom_point() +
        geom_text_repel(data = data.frame(pca.out$x)[outliers,], aes(label = rownames(data.frame(pca.out$x))[outliers])) +
        labs(x = "PC1", y = "PC2", title = "First Two Principal Components with Outliers Highlighted")


# no outliers in GSE47881

GSE47881_exprs <- log2(exprs(GSE47881_eset))


# clean up col names/sample IDs så they are just the GSM kode

colnames(GSE47881_exprs) <- sapply(strsplit(colnames(GSE47881_exprs), split = "_"), "[", 1)



# identify pre-post samples

metadata_GSE47881 %>% 
        mutate(age = as.numeric(`age:ch1`),
               FP = `patientid:ch1`,
               timepoint = ifelse(`time:ch1` == "pre-training", "Pre", "Post")) %>% 
        dplyr::select(FP, timepoint, age) %>% 
        filter(age < 35, FP != "NB021") %>% 
        mutate(Batch = "GSE47881")->  df_GSE47881

# create designmatrix for ebayes model

FP <- factor(df$FP)
timepoint <- factor(df$timepoint, levels = c("Pre", "Post"))


design <- model.matrix(~FP+timepoint)



GSE47881_exprs %>% 
        as.data.frame() %>% 
        dplyr::select(rownames(df)) %>% 
        lmFit(., design) -> fit


fit <- eBayes(fit)



GSE47881_res <- topTable(fit, coef = "timepointPost",number = Inf) %>% 
        filter(P.Value < 0.05) %>% 
        rownames_to_column(var = "affy_hg_u133_plus_2") %>% 
        merge(.,annotLookup, by = "affy_hg_u133_plus_2") %>% 
        arrange(-abs(logFC)) 

GSE47881_res.f <- GSE47881_res[!duplicated(GSE47881_res$logFC),]





saveRDS(GSE47881_res.f, file = "./GSE47881_res.f.RDATA")


#################################################################################################


# study 3 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
# Bamman et al. 2007 / Thalacker-Mercer et al. 2013

#################################################################################################



GSE42507 <- getGEO(GEO = "GSE42507", GSEMatrix = TRUE)

metadata_GSE42507 <- pData(phenoData(GSE42507[[1]]))



write.csv(metadata_GSE42507, file = "/Users/maxul/Downloads/GSE42507.csv")


# there was only baseline data in GEO, waiting for answer from Dan/Adam incase there is pre-post data




#################################################################################################


# MetaMex studies

# GSE24235 # post is 24h after acute RT, but will include
# GSE28422 # already done
# GSE28422 # already done
# GSE28998 # baseline and post samples are 4h post acute exercise, so will skip this one
# EMEXP740
# GSE106865 # done
# GSE8479 # done (only older subjects have pre-post samples)

#################################################################################################




#################################################################################################

#metamex 1: GSE24235

#################################################################################################



GSE24235 <- getGEO(GEO = "GSE24235", GSEMatrix = TRUE)

metadata_GSE24235 <- pData(phenoData(GSE24235[[1]]))





# untar the file and save in directory

untar(tarfile = "/Users/maxul/Downloads/GSE24235_RAW.tar",
      exdir = "/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/GSE24235")

# unpack .gz files

# get list of files

gz_files <- list.files("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/GSE24235", full.names = TRUE)


for (i in 1:length(gz_files)) {
        gunzip(filename = gz_files[i], remove = TRUE)
        print(i)
}

# list .CEL files

cel_files <- list.files("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/GSE24235", 
                        full.names = TRUE, 
                        pattern = "*.CEL")


# read data GSE28422

GSE24235_raw <- read.affybatch(cel_files)


GSE24235_eset <- expresso(GSE24235_raw,
                          bg.correct = TRUE, bgcorrect.method = "rma",               # RMA is background correct method correcting probe intensities by a global method
                          normalize.method = "quantiles",
                          pmcorrect.method = "pmonly",
                          summary.method = "avgdiff")                               # computed average


# create expression dataframe

exprs(GSE24235_eset)

plotDensities(log2(exprs(GSE24235_eset)), legend = FALSE)


# plot pca


pca.out <- prcomp(t(log2(exprs(GSE24235_eset))), scale. = FALSE)

# check proportion of variability explained by pca 1-50

# Extract the proportion of variance explained by each principal component
var_explained <- pca.out$sdev^2 / sum(pca.out$sdev^2)

# Create a data frame for plotting
df <- data.frame(Component = 1:length(var_explained), Variance = var_explained)

# Keep only the first 50 principal components
df <- df[1:50,]

# Plot

ggplot(df, aes(x = Component, y = Variance)) +
        geom_bar(stat = "identity") +
        scale_x_continuous(breaks = 1:50) +
        labs(x = "Principal Component", y = "Proportion of Variance Explained",
             title = "Variance Explained by Principal Components")


plot(pca.out$x[,1:2])


# check for outliers

# Compute Mahalanobis distances for the first two principal components
distances <- mahalanobis(pca.out$x[,1:2], colMeans(pca.out$x[,1:2]), cov(pca.out$x[,1:2]))

# Identify outliers as samples with a Mahalanobis distance greater than a certain threshold
# Here, I'm using the 97.5 percentile of the Chi-square distribution with 3 degrees of freedom as the threshold
outliers <- which(distances > qchisq(0.975, df = 3))

# Print the row names of the outliers
print(colnames(exprs(GSE24235_eset))[outliers])

# Plot the first two principal components, highlighting the outliers

ggplot(data.frame(pca.out$x), aes(x = PC1, y = PC2)) +
        geom_point() +
        geom_text_repel(data = data.frame(pca.out$x)[outliers,], aes(label = rownames(data.frame(pca.out$x))[outliers])) +
        labs(x = "PC1", y = "PC2", title = "First Two Principal Components with Outliers Highlighted")


# no outliers in GSE24235

GSE24235_exprs <- log2(exprs(GSE24235_eset))


# clean up col names/sample IDs så they are just the GSM kode

colnames(GSE24235_exprs) <- gsub(".CEL", "", colnames(GSE24235_exprs))


sapply(strsplit(colnames(GSE24235_exprs), split = "."), "[", 1)



# identify pre-post samples

metadata_GSE24235 %>% 
        mutate(FP = sapply(strsplit(title, split = "_"),"[", 3),
               FP = paste("FP", substr(FP, start = 1, stop = 4), sep = "")) %>% 
        mutate(timepoint = ifelse(`condition:ch1` == "24h post acute resistance exercise following 12-week resistance training",
                                  "Post", ifelse(`condition:ch1` == "resting", "Pre", "4h"))) %>% 
        filter(timepoint != "4h") %>% 
        group_by(FP) %>% 
        filter(any(timepoint == "Pre") & any(timepoint == "Post")) %>% 
        ungroup() %>% 
        dplyr::select(geo_accession, FP, timepoint) %>% 
        mutate(Batch = "GSE24235") %>% 
        as.data.frame()-> df_GSE24235
        
rownames(df_GSE24235) <- df_GSE24235$geo_accession


# create designmatrix for ebayes model

FP <- factor(df$FP)
timepoint <- factor(df$timepoint, levels = c("Pre", "Post"))


design <- model.matrix(~FP+timepoint)



GSE24235_exprs %>% 
        as.data.frame() %>% 
        dplyr::select(df$geo_accession) %>% 
        lmFit(., design) -> fit


fit <- eBayes(fit)



GSE24235_res <- topTable(fit, coef = "timepointPost",number = Inf) %>% 
        filter(P.Value < 0.05) %>% 
        rownames_to_column(var = "affy_hg_u133_plus_2") %>% 
        merge(.,annotLookup, by = "affy_hg_u133_plus_2") %>% 
        arrange(-abs(logFC)) 

GSE24235_res.f <- GSE24235_res[!duplicated(GSE24235_res$logFC),]






saveRDS(GSE24235_res.f, file = "./GSE24235_res.f.RDATA")










#################################################################################################

#metamex 4: GSE28998

#################################################################################################



GSE28998 <- getGEO(GEO = "GSE28998", GSEMatrix = TRUE)

metadata_GSE28998 <- pData(phenoData(GSE28998[[1]]))

# baseline and post samples are 4h post acute exercise, so will skip this one



#################################################################################################

#metamex 5: GSE106865

#################################################################################################



GSE106865 <- getGEO(GEO = "GSE106865", GSEMatrix = TRUE)

metadata_GSE106865 <- pData(phenoData(GSE106865[[1]]))


# read raw .idat files


untar("/Users/maxul/Downloads/GSE106865_RAW.tar", exdir = "./GSE106865/")

files <- list.files("./GSE106865/", full.names = TRUE)

for (i in 1:length(files)) {
        gunzip(filename = files[i], remove = TRUE)
        print(i)
        
}

idat_files <- list.files("./GSE106865/", pattern = "*.idat", full.names = TRUE)

GSE106865_idat <- read.idat(idat_files, bgxfile = "HumanHT-12_V4_0_R2_15002873_B.bgx", bgxpath = "/Users/maxul/Downloads/")


GSE106865_idat$E %>% colnames() -> x
        sapply(strsplit(x, split = "/"),"[", 3) -> x
                sapply(strsplit(x, split = "_"),"[",1)  -> names


colnames(GSE106865_idat$E) <- names

GSE106865_exprs <-  neqc(GSE106865_idat)  # normalize and background correct expression values, and log2 transformed


GSE106865_exprs$E







plotDensities(GSE106865_exprs$E, legend = FALSE)


# plot pca


pca.out <- prcomp(t(GSE106865_exprs$E), scale. = FALSE)

# check proportion of variability explained by pca 1-50

# Extract the proportion of variance explained by each principal component
var_explained <- pca.out$sdev^2 / sum(pca.out$sdev^2)

# Create a data frame for plotting
df <- data.frame(Component = 1:length(var_explained), Variance = var_explained)

# Keep only the first 50 principal components
df <- df[1:50,]

# Plot

ggplot(df, aes(x = Component, y = Variance)) +
        geom_bar(stat = "identity") +
        scale_x_continuous(breaks = 1:50) +
        labs(x = "Principal Component", y = "Proportion of Variance Explained",
             title = "Variance Explained by Principal Components")


plot(pca.out$x[,1:2])


# check for outliers

# Compute Mahalanobis distances for the first two principal components
distances <- mahalanobis(pca.out$x[,1:2], colMeans(pca.out$x[,1:2]), cov(pca.out$x[,1:2]))

# Identify outliers as samples with a Mahalanobis distance greater than a certain threshold
# Here, I'm using the 97.5 percentile of the Chi-square distribution with 3 degrees of freedom as the threshold
outliers <- which(distances > qchisq(0.975, df = 3))

# Print the row names of the outliers
print(colnames(GSE106865_exprs$E)[outliers])

# Plot the first two principal components, highlighting the outliers

ggplot(data.frame(pca.out$x), aes(x = PC1, y = PC2)) +
        geom_point() +
        geom_text_repel(data = data.frame(pca.out$x)[outliers,], aes(label = rownames(data.frame(pca.out$x))[outliers])) +
        labs(x = "PC1", y = "PC2", title = "First Two Principal Components with Outliers Highlighted")


# remove outliers

GSE28422_exprs <- log2(exprs(GSE28422_eset))

GSE106865_exprs.f <- GSE106865_exprs$E[,!(colnames(GSE106865_exprs$E) %in% names(outliers))]



pca.out <- prcomp(t(GSE106865_exprs.f), scale. = FALSE)

plot(pca.out$x[,1:2])



# get the GSM's for pre post samples

metadata_GSE106865 %>% 
        dplyr::select(1) %>% 
        mutate(FP = paste("FP", sapply(strsplit(title, split = "_"), "[", 5), sep = ""),
               timepoint = sapply(strsplit(title, split = "_"), "[", 7),
               timepoint = ifelse(timepoint == "pre", "Pre", "Post"),
               trainingstatus = sapply(strsplit(title, split = "_"), "[", 6)) %>% 
        filter(timepoint == "Pre") %>% 
        mutate(timepoint = factor(ifelse(trainingstatus == "untrained", "Pre", "Post"), levels = c("Pre", "Post"))) %>%
        dplyr::select(2,3) %>% 
        mutate(Batch = "GSE106865") -> df_GSE106865







FP <- factor(df$FP)
timepoint <- factor(df$timepoint, levels = c("Pre", "Post"))


design <- model.matrix(~FP+timepoint)


GSE106865_exprs$E %>% 
        as.data.frame() %>% 
        dplyr::select(rownames(df)) %>% 
        lmFit(., design) -> fit


fit <- eBayes(fit)



GSE106865_res <- topTable(fit, coef = "timepointPost",number = Inf, adjust.method = "BH") %>% 
        filter(P.Value < 0.05) %>%  
        arrange(-abs(logFC)) 


# annotate illumina identifyers to gene name and ensembl

BiocManager::install("illuminaHumanv4.db")

library(illuminaHumanv4.db)

illumina_keys <- as.data.frame(illuminaHumanv4ENSEMBLREANNOTATED)

GSE106865_res %>% 
        mutate(IlluminaID = paste("ILMN_",rownames(GSE106865_res), sep = "")) %>% 
        merge(.,illumina_keys, by = "IlluminaID")





# annotate with older package based on https://doi.org/10.1093/nar/gkp942

BiocManager::install("illuminaHumanv4BeadID.db")
library(illuminaHumanv2BeadID.db)


x <- as.data.frame(illuminaHumanv2BeadIDALIAS2PROBE)


GSE106865_res %>% 
        mutate(probe_id = rownames(GSE106865_res)) %>% 
        merge(.,x, by = "probe_id") -> GSE106865_res.f


saveRDS(GSE106865_res.f, file = "./GSE106865_res.f.RDATA")


#################################################################################################

#metamex 6: GSE8479 (only older subjects have pre-post samples)

#################################################################################################


GSE8479 <- getGEO(GEO = "GSE8479", GSEMatrix = TRUE)

metadata_GSE8479 <- pData(phenoData(GSE8479[[1]]))




GSE8479_csv <- read.csv("/Users/maxul/Downloads/GSE8479_NonNormalizedMelovetalCoded.csv")

# format csv data for normalization and background correction with neqc function

exprs <- as.matrix(GSE8479_csv[,grep("AVG_Signal", colnames(GSE8479_csv))])

detP <- as.matrix(GSE8479_csv[,grep("Detection.Pval", colnames(GSE8479_csv))])


GSE8479_exprs <- neqc(exprs, detection.p = detP)

rownames(GSE8479_exprs) <- GSE8479_csv[,3] 


plotDensities(GSE8479_exprs, legend = FALSE)







# plot pca


pca.out <- prcomp(t(GSE8479_exprs), scale. = FALSE)

# check proportion of variability explained by pca 1-50

# Extract the proportion of variance explained by each principal component
var_explained <- pca.out$sdev^2 / sum(pca.out$sdev^2)

# Create a data frame for plotting
df <- data.frame(Component = 1:length(var_explained), Variance = var_explained)

# Keep only the first 50 principal components
df <- df[1:50,]

# Plot

ggplot(df, aes(x = Component, y = Variance)) +
        geom_bar(stat = "identity") +
        scale_x_continuous(breaks = 1:50) +
        labs(x = "Principal Component", y = "Proportion of Variance Explained",
             title = "Variance Explained by Principal Components")


plot(pca.out$x[,1:2])


# check for outliers

# Compute Mahalanobis distances for the first two principal components
distances <- mahalanobis(pca.out$x[,1:2], colMeans(pca.out$x[,1:2]), cov(pca.out$x[,1:2]))

# Identify outliers as samples with a Mahalanobis distance greater than a certain threshold
# Here, I'm using the 97.5 percentile of the Chi-square distribution with 3 degrees of freedom as the threshold
outliers <- which(distances > qchisq(0.975, df = 3))

# Print the row names of the outliers
print(colnames(GSE8479_exprs)[outliers])

# Plot the first two principal components, highlighting the outliers

ggplot(data.frame(pca.out$x), aes(x = PC1, y = PC2)) +
        geom_point() +
        geom_text_repel(data = data.frame(pca.out$x)[outliers,], aes(label = rownames(data.frame(pca.out$x))[outliers])) +
        labs(x = "PC1", y = "PC2", title = "First Two Principal Components with Outliers Highlighted")



# weird pca plot results (clustering) 

# label samples to check for patterns






ggplot(data.frame(pca.out$x), aes(x = PC1, y = PC2)) +
        geom_point() +
        geom_text_repel(data = data.frame(pca.out$x), aes(label = rownames(data.frame(pca.out$x)))) +
        labs(x = "PC1", y = "PC2", title = "First Two Principal Components with Outliers Highlighted")


# post exercise is clusterd together... 


metadata_GSE8479 <- metadata_GSE8479 %>% 
        mutate(ID = sapply(strsplit(title, split =  " "), "[", 3)) %>% 
        mutate(timepoint = ifelse((grepl("EB", ID)), "Post", "Pre")) %>%  
        mutate(ID = gsub("EB", "", ID))

# get ID of post samples

metadata_GSE8479 %>% 
        filter(timepoint == "Post") %>% 
        pull(ID) -> x

metadata_GSE8479 %>% 
        filter(ID %in% x) %>% rownames()-> z


# relabel expression set to GSM code


metadata_GSE8479 %>% 
        mutate(ID = sapply(strsplit(title, split =  " "), "[", 3)) -> df


colnames(GSE8479_exprs) %>% 
        as.data.frame() %>% 
        mutate(ID = gsub("AVG_Signal", "", .),
               ID = gsub("\\.", "", ID)) %>% 
        merge(df, ., by = "ID") %>% 
        dplyr::select(geo_accession,"colnames(GSE8479_exprs)"  = ".") %>% 
        merge(as.data.frame(colnames(GSE8479_exprs)), ., by = "colnames(GSE8479_exprs)") -> df2
        


# Assuming your first dataframe is named df2
new_colnames <- df2$geo_accession[match(colnames(GSE8479_exprs), df2$'colnames(GSE8479_exprs)')]

# Rename the columns of GSE8479_exprs
colnames(GSE8479_exprs) <- new_colnames



# keep only pre-post samples

GSE8479_exprs.f <- GSE8479_exprs[,z]


# rerun pca on only pre-post samples


pca.out <- prcomp(t(GSE8479_exprs.f), scale. = FALSE)

plot(pca.out$x[,1:2])




metadata_GSE8479 %>% 
        filter(geo_accession %in% z) %>% 
        dplyr::select("FP" = ID, timepoint) %>% 
        mutate(Batch = "GSE8479") -> df_GSE8479


        


FP <- factor(metadata_GSE8479[z,]$ID)
timepoint <- factor(metadata_GSE8479[z,]$timepoint, levels = c("Pre", "Post"))


design <- model.matrix(~FP+timepoint)


GSE8479_exprs.f %>% 
        as.data.frame() %>% 
        lmFit(., design) -> fit


fit <- eBayes(fit)



GSE8479_res <- topTable(fit, coef = "timepointPost",number = Inf, adjust.method = "BH") %>% 
        filter(P.Value < 0.05) %>%  
        arrange(-abs(logFC)) 




saveRDS(GSE8479_res, file = "./GSE8479_res.RDATA")


##########################################################################################

# GSE154846: Timmons and Atherton 2020

##########################################################################################

# https://doi.org/10.1371/journal.pgen.1003389

BiocManager::install("oligo")

BiocManager::install("pd.hta.2.0")

library(oligo)


GSE154846 <- getGEO(GEO = "GSE154846", GSEMatrix = TRUE)

metadata_GSE154846 <- pData(phenoData(GSE154846[[1]]))

metadata_GSE154846_2 <- pData(phenoData(GSE154846[[2]]))




metadata_GSE154846_2 %>% 
        dplyr::select(1, characteristics_ch1.1) %>% 
        filter(nchar(title) > 4 & nchar(title) <9) %>% 
        mutate(FP = as.numeric(str_extract(title, "\\d+")))  %>% 
        mutate(timepoint = sapply(strsplit(characteristics_ch1.1, split =  " "), "[", 3))



        # untar the file and save in directory
        
        untar(tarfile = "/Users/maxul/Downloads/GSE154846_RAW.tar",
              exdir = "/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/GSE154846")
        
        # unpack .gz files
        
        # get list of files
        
        gz_files <- list.files("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/GSE154846", full.names = TRUE)
        
        
        for (i in 1:length(gz_files)) {
                gunzip(filename = gz_files[i], remove = TRUE)
                print(i)
        }
        
        # list .CEL files
        
        cel_files <- list.files("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Pooled transcriptomics/GSE154846", 
                                full.names = TRUE, 
                                pattern = "*.CEL")
        
        
        # read data GSE28422

        GSE154846_raw <- oligo::read.celfiles(cel_files)

        
        
        GSE154846_eset <- rma(GSE154846_raw, target = "core")
        
        # create expression dataframe
        
        exprs(GSE154846_eset)
        
        plotDensities(log2(exprs(GSE154846_eset)), legend = FALSE)
        
        
        # plot pca
        
        
        pca.out <- prcomp(t(log2(exprs(GSE154846_eset))), scale. = FALSE)
        
        # check proportion of variability explained by pca 1-50
        
        # Extract the proportion of variance explained by each principal component
        var_explained <- pca.out$sdev^2 / sum(pca.out$sdev^2)
        
        # Create a data frame for plotting
        df <- data.frame(Component = 1:length(var_explained), Variance = var_explained)
        
        # Keep only the first 50 principal components
        df <- df[1:50,]
        
        # Plot
        
        ggplot(df, aes(x = Component, y = Variance)) +
                geom_bar(stat = "identity") +
                scale_x_continuous(breaks = 1:50) +
                labs(x = "Principal Component", y = "Proportion of Variance Explained",
                     title = "Variance Explained by Principal Components")
        
        
        plot(pca.out$x[,1:2])
        
        
        # check for outliers
        
        # Compute Mahalanobis distances for the first two principal components
        distances <- mahalanobis(pca.out$x[,1:2], colMeans(pca.out$x[,1:2]), cov(pca.out$x[,1:2]))
        
        # Identify outliers as samples with a Mahalanobis distance greater than a certain threshold
        # Here, I'm using the 97.5 percentile of the Chi-square distribution with 3 degrees of freedom as the threshold
        outliers <- which(distances > qchisq(0.975, df = 3))
        
        # Print the row names of the outliers
        print(colnames(exprs(GSE154846_eset))[outliers])
        
        # Plot the first two principal components, highlighting the outliers
        
        ggplot(data.frame(pca.out$x), aes(x = PC1, y = PC2)) +
                geom_point() +
                geom_text_repel(data = data.frame(pca.out$x)[outliers,], aes(label = rownames(data.frame(pca.out$x))[outliers])) +
                labs(x = "PC1", y = "PC2", title = "First Two Principal Components with Outliers Highlighted")
        
        
        # no outliers in GSE154846
        
        GSE154846_exprs <- log2(exprs(GSE154846_eset))
        
        
        # clean up col names/sample IDs så they are just the GSM kode
        
        colnames(GSE154846_exprs) <- gsub(".CEL", "", colnames(GSE154846_exprs))
        
        
        colnames(GSE154846_exprs) <- sapply(strsplit(colnames(GSE154846_exprs), split = "_"), "[", 2)
        
        

        # identify pre-post samples
        
        metadata_GSE154846_2 %>% 
                dplyr::select(1, characteristics_ch1.1) %>% 
                filter(nchar(title) > 4 & nchar(title) <9) %>% 
                mutate(FP = as.numeric(str_extract(title, "\\d+")))  %>% 
                mutate(timepoint = sapply(strsplit(characteristics_ch1.1, split =  " "), "[", 3)) %>% 
                filter(timepoint != "POST_Unloading") %>% 
                mutate(Batch = "GSE154846")-> df_GSE154846
        

        
        
        # create designmatrix for ebayes model
        
        FP <- factor(df$FP)
        timepoint <- factor(df$timepoint, levels = c("PRE", "POST"))
        
        
        design <- model.matrix(~FP+timepoint)
        
        
        
        GSE154846_exprs %>% 
                as.data.frame() %>% 
                dplyr::select(df$title) %>% 
                lmFit(., design) -> fit
        
        
        fit <- eBayes(fit)
        
        
        
        GSE154846_res <- topTable(fit, coef = "timepointPOST",number = Inf) %>% 
                filter(P.Value < 0.05) %>% 
                rownames_to_column(var = "PROBEID") %>% 
                merge(., HTA_anno, by = "PROBEID")%>% 
                arrange(-abs(logFC)) 
        

        
        
        
        
        
        saveRDS(GSE154846_res, file = "./GSE154846_res.RDATA")
        
        
        

BiocManager::install("affycoretools")

library(hta20transcriptcluster.db)
library(affycoretools)

GSE154846_eset2 <- annotateEset(GSE154846_eset, hta20transcriptcluster.db)


fData(GSE154846_eset2) %>% 
        na.omit() -> HTA_anno


## change probe annotation (row name to SYMBOL)

GSE154846_exprs %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "PROBEID") %>% 
        merge(., HTA_anno, by = "PROBEID") %>% 
        na.omit() %>% 
        distinct(SYMBOL, .keep_all = TRUE) %>% 
        dplyr::select(df_GSE154846$title, "SYMBOL") -> GSE154846_df


rownames(GSE154846_df) <- GSE154846_df$SYMBOL

GSE154846_df <- GSE154846_df[,(-ncol(GSE154846_df))]

##########################################################################################

# summarize pre-post results in venn diagram

##########################################################################################


library(ggvenn)


# make list with GSE names and Gene symbol names


venn_df <- list(GSE106865 = GSE106865_res.f %>% filter(abs(logFC) > 0) %>% pull(alias_symbol),
                GSE24235 = GSE24235_res.f %>% filter(abs(logFC) > 0)  %>% pull(external_gene_name),
                GSE28422 = GSE28422_res.f%>% filter(abs(logFC) > 0)  %>% pull(external_gene_name),
                GSE47881 = GSE47881_res.f%>% filter(abs(logFC) > 0)  %>% pull(external_gene_name))


ggvenn(venn_df, text_size = 8)


# remake venn with negative fold change


venn_df_neg <- list(GSE106865 = GSE106865_res.f %>% filter(logFC <= 0) %>% pull(alias_symbol)%>% unique(),
                GSE24235 = GSE24235_res.f %>% filter(logFC <= 0) %>% pull(external_gene_name)%>% unique(),
                GSE28422 = GSE28422_res.f%>% filter(logFC <= 0) %>% pull(external_gene_name)%>% unique(),
                GSE47881 = GSE47881_res.f%>% filter(logFC <= 0) %>% pull(external_gene_name)%>% unique())

ggvenn(venn_df_neg)

# remake venn with positive fold change


venn_df_pos <- list(GSE106865 = GSE106865_res.f %>% filter(logFC >= 0) %>% pull(alias_symbol)%>% unique(),
                    GSE24235 = GSE24235_res.f %>% filter(logFC >= 0) %>% pull(external_gene_name)%>% unique(),
                    GSE28422 = GSE28422_res.f%>% filter(logFC >= 0) %>% pull(external_gene_name)%>% unique(),
                    GSE47881 = GSE47881_res.f%>% filter(logFC >= 0) %>% pull(external_gene_name)%>% unique())

ggvenn(venn_df_pos)



# Compute overlaps
overlaps <- list(
        A = setdiff(venn_df$GSE106865, unlist(venn_df[-1])),
        B = setdiff(venn_df$GSE24235, unlist(venn_df[c(-2)])),
        C = setdiff(venn_df$GSE28422, unlist(venn_df[c(-3)])),
        D = setdiff(venn_df$GSE47881, unlist(venn_df[c(-4)])),
        AB = intersect(venn_df$GSE106865, venn_df$GSE24235),
        AC = intersect(venn_df$GSE106865, venn_df$GSE28422),
        AD = intersect(venn_df$GSE106865, venn_df$GSE47881),
        BC = intersect(venn_df$GSE24235, venn_df$GSE28422),
        BD = intersect(venn_df$GSE24235, venn_df$GSE47881),
        CD = intersect(venn_df$GSE28422, venn_df$GSE47881),
        ABC = Reduce(intersect, venn_df[c("GSE106865", "GSE24235", "GSE28422")]),
        ABD = Reduce(intersect, venn_df[c("GSE106865", "GSE24235", "GSE47881")]),
        ACD = Reduce(intersect, venn_df[c("GSE106865", "GSE28422", "GSE47881")]),
        BCD = Reduce(intersect, venn_df[c("GSE24235", "GSE28422", "GSE47881")]),
        ABCD = Reduce(intersect, venn_df)
)



# Compute overlaps
overlaps_neg <- list(
        A = setdiff(venn_df_neg$GSE106865, unlist(venn_df_neg[-1])),
        B = setdiff(venn_df_neg$GSE24235, unlist(venn_df_neg[c(-2)])),
        C = setdiff(venn_df_neg$GSE28422, unlist(venn_df_neg[c(-3)])),
        D = setdiff(venn_df_neg$GSE47881, unlist(venn_df_neg[c(-4)])),
        AB = intersect(venn_df_neg$GSE106865, venn_df_neg$GSE24235),
        AC = intersect(venn_df_neg$GSE106865, venn_df_neg$GSE28422),
        AD = intersect(venn_df_neg$GSE106865, venn_df_neg$GSE47881),
        BC = intersect(venn_df_neg$GSE24235, venn_df_neg$GSE28422),
        BD = intersect(venn_df_neg$GSE24235, venn_df_neg$GSE47881),
        CD = intersect(venn_df_neg$GSE28422, venn_df_neg$GSE47881),
        ABC = Reduce(intersect, venn_df_neg[c("GSE106865", "GSE24235", "GSE28422")]),
        ABD = Reduce(intersect, venn_df_neg[c("GSE106865", "GSE24235", "GSE47881")]),
        ACD = Reduce(intersect, venn_df_neg[c("GSE106865", "GSE28422", "GSE47881")]),
        BCD = Reduce(intersect, venn_df_neg[c("GSE24235", "GSE28422", "GSE47881")]),
        ABCD = Reduce(intersect, venn_df_neg)
)


# Compute overlaps
overlaps_pos <- list(
        A = setdiff(venn_df_pos$GSE106865, unlist(venn_df_pos[-1])),
        B = setdiff(venn_df_pos$GSE24235, unlist(venn_df_pos[c(-2)])),
        C = setdiff(venn_df_pos$GSE28422, unlist(venn_df_pos[c(-3)])),
        D = setdiff(venn_df_pos$GSE47881, unlist(venn_df_pos[c(-4)])),
        AB = intersect(venn_df_pos$GSE106865, venn_df_pos$GSE24235),
        AC = intersect(venn_df_pos$GSE106865, venn_df_pos$GSE28422),
        AD = intersect(venn_df_pos$GSE106865, venn_df_pos$GSE47881),
        BC = intersect(venn_df_pos$GSE24235, venn_df_pos$GSE28422),
        BD = intersect(venn_df_pos$GSE24235, venn_df_pos$GSE47881),
        CD = intersect(venn_df_pos$GSE28422, venn_df_pos$GSE47881),
        ABC = Reduce(intersect, venn_df_pos[c("GSE106865", "GSE24235", "GSE28422")]),
        ABD = Reduce(intersect, venn_df_pos[c("GSE106865", "GSE24235", "GSE47881")]),
        ACD = Reduce(intersect, venn_df_pos[c("GSE106865", "GSE28422", "GSE47881")]),
        BCD = Reduce(intersect, venn_df_pos[c("GSE24235", "GSE28422", "GSE47881")]),
        ABCD = Reduce(intersect, venn_df_pos)
)


overlaps$ABCD
overlaps_neg$ABCD
overlaps_pos$ABCD


GSE47881_res.f %>% 
        filter(external_gene_name %in% overlaps$ABCD)

GSE24235_res.f %>% 
        filter(external_gene_name %in% overlaps$ABCD)


sGSE28422_res.f %>% 
        filter(external_gene_name %in% overlaps$ABCD)

GSE106865_res.f %>% 
        filter(alias_symbol %in% overlaps$ABCD)

GSE24235_res.f %>% 
        filter(external_gene_name == "LDLRAD4")



venn_df <- list(GSE24235 = GSE24235_res.f %>% filter(abs(logFC) > 0)  %>% pull(external_gene_name),
                GSE28422 = GSE28422_res.f%>% filter(abs(logFC) > 0)  %>% pull(external_gene_name),
                GSE47881 = GSE47881_res.f%>% filter(abs(logFC) > 0)  %>% pull(external_gene_name))


ggvenn(venn_df, text_size = 8)

overlaps_neg$BCD
overlaps_pos$BCD




# extract overlaps that overlap between at least two studies

overlaps_2 <- paste(overlaps$AB, overlaps$AC, overlaps$AD, overlaps$BC, overlaps$BD, overlaps$CD) %>% 
        strsplit(split = " ") %>% 
        unlist() %>% 
        unique() 

library(gage)
KEGG_new <- kegg.gsets(species = "hsa", id.type = "kegg", check.new=TRUE)

# merge all datasets into expression set: genes as rownames and columns as samples

annotLookup %>% 
        dplyr::select(ENTREZID,"alias_symbol" = external_gene_name) -> x

GSE106865_res.f %>% 
        filter(alias_symbol %in% overlaps_2) %>% 
        merge(., x, by = "alias_symbol") %>% 
        dplyr::select(ENTREZID, logFC) %>% 
        distinct(ENTREZID, .keep_all = TRUE) %>% 
        na.omit() -> z
        
GSE24235_res.f %>% 
        mutate(alias_symbol = external_gene_name) %>% 
        filter(alias_symbol %in% overlaps_2) %>% 
        merge(., annotLookup, by = "ensembl_gene_id") %>% 
        dplyr::select(ENTREZID, logFC) %>% 
        distinct(ENTREZID, .keep_all = TRUE) %>% 
        na.omit() -> z


GSE28422_res.f %>% 
        mutate(alias_symbol = external_gene_name) %>% 
        filter(alias_symbol %in% overlaps_2) %>% 
        merge(., annotLookup, by = "ensembl_gene_id") %>% 
        dplyr::select(ENTREZID, logFC) %>% 
        distinct(ENTREZID, .keep_all = TRUE) %>% 
        na.omit() -> z


GSE47881_res.f %>% 
        mutate(alias_symbol = external_gene_name) %>% 
        filter(alias_symbol %in% overlaps_2) %>% 
        merge(., annotLookup, by = "ensembl_gene_id") %>% 
        dplyr::select(ENTREZID, logFC) %>% 
        distinct(ENTREZID, .keep_all = TRUE) %>% 
        na.omit() -> z

rownames(z) <- z$ENTREZID


z <- z[-1]


library(org.Hs.eg.db)
data(kegg.sets.hs)


### homogenate

exprsMat <- as.matrix(z) 

subset <- KEGG_new$kg.sets[KEGG_new[["sigmet.idx"]]]



kegg_res_3 <- gage(exprs = exprsMat, gsets = subset, same.dir = TRUE, ref = NULL, samp = NULL)

view(kegg_res_4$less)                              ### view less methylated KEGG pathways post vs. Baseline
view(kegg_res_1$greater)     


# did not work, will do two and two

kegg_res_1$greater %>% 
        head(6) %>% 
        rbind(.,kegg_res_2$greater %>% 
                      head(10)) %>% 
        rbind(.,kegg_res_3$greater %>% 
                      head(10)) %>% 
        rbind(.,kegg_res_4$greater %>% 
                      head(10)) %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "pathway") %>% 
        distinct(pathway, .keep_all = TRUE) -> top_pathways


kegg_res_2$greater %>% 
        head(10)

kegg_res_3$greater %>% 
        head(10)

kegg_res_4$greater %>% 
        head(10)







##########################################################################################

# combine and run batch correct and pre-post diffrential expression

##########################################################################################


# start by adding SYMBOL gene name to all expression sets



GSE106865_exprs.f %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "probe_id") %>%
        merge(., x, by = "probe_id") %>% 
        distinct(alias_symbol, .keep_all = TRUE) %>% 
        dplyr::select(rownames(df_GSE106865), "SYMBOL" = alias_symbol) -> df1

GSE154846_exprs %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "PROBEID") %>% 
        merge(., HTA_anno, by = "PROBEID") %>% 
        na.omit() %>% 
        distinct(SYMBOL, .keep_all = TRUE) %>% 
        dplyr::select(df_GSE154846$title, "SYMBOL") -> df2

colnames(df2) <- c(rownames(df_GSE154846), "SYMBOL")

GSE24235_exprs %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "affy_hg_u133_plus_2") %>% 
        merge(.,annotLookup, by = "affy_hg_u133_plus_2") %>% 
        dplyr::select(rownames(df_GSE24235), "SYMBOL" = external_gene_name) %>% 
        distinct(SYMBOL, .keep_all = TRUE) %>% 
        filter(SYMBOL != "") -> df3

GSE28422_exprs %>% 
        as.data.frame() %>% 
        dplyr::select(rownames(df_GSE28422)) %>% 
        rownames_to_column(var = "affy_hg_u133_plus_2") %>% 
        merge(.,annotLookup, by = "affy_hg_u133_plus_2") %>% 
        dplyr::select(rownames(df_GSE28422), "SYMBOL" = external_gene_name) %>% 
        distinct(SYMBOL, .keep_all = TRUE) %>% 
        filter(SYMBOL != "") -> df4

GSE47881_exprs %>% 
        as.data.frame() %>% 
        dplyr::select(rownames(df_GSE47881)) %>% 
        rownames_to_column(var = "affy_hg_u133_plus_2") %>% 
        merge(.,annotLookup, by = "affy_hg_u133_plus_2") %>% 
        dplyr::select(rownames(df_GSE47881), "SYMBOL" = external_gene_name) %>% 
        distinct(SYMBOL, .keep_all = TRUE) %>% 
        filter(SYMBOL != "") -> df5

GSE8479_exprs %>% 
        as.data.frame() %>% 
        dplyr::select(rownames(df_GSE8479)) %>% 
        rownames_to_column(var = "SYMBOL") -> df6


# merge all df 1-6 and pheno data files

df1 %>% 
        merge(., df2, by = "SYMBOL") %>% 
        merge(., df3, by = "SYMBOL") %>% 
        merge(., df4, by = "SYMBOL") %>% 
        merge(., df5, by = "SYMBOL") %>% 
        merge(., df6, by = "SYMBOL") -> pooled_matrix_7610_genes_120_samples 

rownames(pooled_matrix_7610_genes_120_samples) <- pooled_matrix_7610_genes_120_samples$SYMBOL
pooled_matrix_7610_genes_120_samples <- pooled_matrix_7610_genes_120_samples[,-1]

# pool matrix without df6

df1 %>% 
        merge(., df2, by = "SYMBOL") %>% 
        merge(., df3, by = "SYMBOL") %>% 
        merge(., df4, by = "SYMBOL") %>% 
        merge(., df5, by = "SYMBOL") -> pooled_matrix_13014_genes_92_samples


df2 %>% 
        merge(., df3, by = "SYMBOL") %>% 
        merge(., df4, by = "SYMBOL") %>% 
        merge(., df5, by = "SYMBOL") -> pooled_matrix_19851_genes_72_samples


rownames(pooled_matrix_19851_genes_72_samples) <- pooled_matrix_19851_genes_72_samples$SYMBOL
pooled_matrix_19851_genes_72_samples <- pooled_matrix_19851_genes_72_samples[,-1]

# clean up annotaton files

df_GSE106865 %>% 
        dplyr::select(FP, timepoint, Batch) -> y

df_GSE154846 %>% 
        dplyr::select(FP, timepoint, Batch) %>% 
        mutate(timepoint = ifelse(timepoint == "PRE", "Pre", "Post")) %>% 
        rbind(y, .) -> y

df_GSE24235 %>% 
        dplyr::select(FP, timepoint, Batch)%>% 
        rbind(y, .) -> y

df_GSE28422 %>% 
        dplyr::select(FP, timepoint, Batch)%>% 
        rbind(y, .) -> y

df_GSE47881 %>% 
        dplyr::select(FP, timepoint, Batch)%>% 
        rbind(y, .) -> y

df_GSE8479 %>% 
        dplyr::select(FP, timepoint, Batch)%>% 
        rbind(y, .) -> y


# save files

saveRDS(pooled_matrix_7610_genes_120_samples, "./pooled_matrix_7610_genes_120_samples.RDATA")
saveRDS(y, "./pooled_pheno_data.RDATA")



pca.out <- prcomp(t(pooled_matrix_7610_genes_120_samples), scale. = FALSE)

plot(pca.out$x[,1:2])

# very pronounced batch effects

# remove batch effects with combat
library(sva)

batch = y$Batch

FP <- factor(y$FP)
timepoint <- factor(y$timepoint, levels = c("Pre", "Post"))


design <- model.matrix(~timepoint)

combat_eDATA <- ComBat(dat = pooled_matrix_7610_genes_120_samples,
                       batch = batch,
                       mod = design)


# rerun pca



pca.out <- prcomp(t(combat_eDATA), scale. = FALSE)

plot(pca.out$x[,1:2])


# batch correct worked perfectly

# run ebayes model


# create designmatrix for ebayes model



design <- model.matrix(~FP+timepoint)



combat_eDATA %>% 
        as.data.frame() %>% 
        dplyr::select(rownames(y)) %>% 
        lmFit(., design) -> fit


fit <- eBayes(fit)



pooled_res <- topTable(fit, coef = "timepointPost",number = Inf) %>% 
        filter(P.Value < 0.05) %>% 
        arrange(-abs(logFC)) 

write.csv(pooled_res, "./pooled_res.csv")


# re_run with max genes


pooled_matrix_19851_genes_72_samples

y[rownames(y) %in% colnames(pooled_matrix_19851_genes_72_samples),] -> y2




batch = factor(y2$Batch)

FP <- factor(y2$FP)
timepoint <- factor(y2$timepoint, levels = c("Pre", "Post"))


design <- model.matrix(~timepoint)

combat_eDATA2 <- ComBat(dat = pooled_matrix_19851_genes_72_samples,
                       batch = batch,
                       mod = design)


# rerun pca



pca.out <- prcomp(t(combat_eDATA2), scale. = FALSE)

plot(pca.out$x[,1:2])


# batch correct worked perfectly

# run ebayes model


# create designmatrix for ebayes model



design <- model.matrix(~FP+timepoint)



combat_eDATA2 %>% 
        as.data.frame() %>% 
        dplyr::select(rownames(y2)) %>% 
        lmFit(., design) -> fit


fit <- eBayes(fit)



pooled_res2 <- topTable(fit, coef = "timepointPost",number = Inf) %>% 
        filter(P.Value < 0.05) %>% 
        arrange(-abs(logFC))


write.csv(pooled_res2, "./pooled_res2.csv")


# run without batch correct, but as variable in design



design <- model.matrix(~FP+timepoint+batch)



pooled_matrix_19851_genes_72_samples %>% 
        as.data.frame() %>% 
        dplyr::select(rownames(y2)) %>% 
        lmFit(., design) -> fit


fit <- eBayes(fit)



pooled_res3 <- topTable(fit, coef = "timepointPost",number = Inf) %>% 
        filter(P.Value < 0.05) %>% 
        arrange(-abs(logFC))


## for some reason "design matrix doesnt work anymore


FP = y$FP
timepoint <- factor(y$timepoint, levels = c("Pre", "Post"))
batch = factor(y$Batch)


design <- model.matrix(~FP+timepoint+batch)

z <-  c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

design <- cbind(design, "batchGSE106865" = z)


pooled_matrix_7610_genes_120_samples %>% 
        as.data.frame() %>% 
        dplyr::select(rownames(y)) %>% 
        lmFit(., design) -> fit


fit <- eBayes(fit)



pooled_res3 <- topTable(fit, coef = "timepointPost",number = Inf) %>% 
        filter(P.Value < 0.05) %>% 
        arrange(-abs(logFC))



saveRDS(pooled_res3, "./pooled_transcriptomics.RDATA")
write.csv(pooled_res3, "./pooled_transcriptomics.csv")




##########################################################################################

# GSEA og pooled genes

##########################################################################################


# annotate with entrezID

pooled_res3 %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "gene_name") -> df_res
        



entrezIDs <- mapIds(org.Hs.eg.db, keys = df_res$gene_name, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

df_res$entrezid <- entrezIDs

library(gage)
KEGG_new <- kegg.gsets(species = "hsa", id.type = "kegg", check.new=TRUE)



rownames(df_res) <- entrezIDs

df_res %>% 
        dplyr::select(logFC) %>% 
        as.matrix() -> exprsMat


kegg_res <- gage(exprs = exprsMat, gsets = KEGG_new$kg.sets, same.dir = TRUE, ref = NULL, samp = NULL)


kegg_res$greater
kegg_res$less

# save kegg sets

write.csv(kegg_res$greater, "./pooled_kegg_greater.csv")
write.csv(kegg_res$less, "./pooled_kegg_less.csv")






##########################################################################################

# overlap with MACS

##########################################################################################

        
pooled_res3 %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "gene_name") -> df_res

DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        mutate(gene_name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        merge(., df_res, by = "gene_name") %>% 
        filter(Relation_to_Island == "Island") %>% 
        ggplot(aes(x = delta_M, y = logFC))+
        geom_point()+
        geom_smooth(method = "lm")+
        geom_text_repel(aes(label = gene_name))



DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        mutate(gene_name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        merge(., df_res, by = "gene_name") %>% 
        filter(Relation_to_Island == "Island") %>% 
        ggplot(aes(x = delta_M, y = logFC))+
        geom_point()+
        geom_smooth(method = "lm")+
        geom_text_repel(aes(label = gene_name))

# visualize hypomethylated and increased gene expression


DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        mutate(gene_name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        merge(., df_res, by = "gene_name") %>% 
        filter(Relation_to_Island == "Island") %>% 
        filter(delta_M <0 & logFC >0) %>% 
        ggplot(aes(x = delta_M, y = logFC))+
        geom_point()+
        geom_smooth(method = "lm")+
        geom_text_repel(aes(label = gene_name))



DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        mutate(gene_name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        merge(., df_res, by = "gene_name") %>% 
        filter(Relation_to_Island == "Island") %>% 
        filter(delta_M <0 & logFC >0) %>% 
        ggplot(aes(x = delta_M, y = logFC))+
        geom_point()+
        geom_smooth(method = "lm")+
        geom_text_repel(aes(label = gene_name))



# visualize hypermethylated and decreased gene expression


DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        mutate(gene_name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        merge(., df_res, by = "gene_name") %>% 
        filter(Relation_to_Island == "Island") %>% 
        filter(delta_M >0 & logFC <0) %>% 
        ggplot(aes(x = delta_M, y = logFC))+
        geom_point()+
        geom_smooth(method = "lm")+
        geom_text_repel(aes(label = gene_name))



DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        mutate(gene_name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        merge(., df_res, by = "gene_name") %>% 
        filter(Relation_to_Island == "Island") %>% 
        filter(delta_M >0 & logFC <0) %>% 
        ggplot(aes(x = delta_M, y = logFC))+
        geom_point()+
        geom_smooth(method = "lm")+
        geom_text_repel(aes(label = gene_name))























