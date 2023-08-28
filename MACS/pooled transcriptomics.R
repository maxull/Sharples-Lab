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
                "external_gene_name"),
        filter = "affy_hg_u133_plus_2",
        values = rownames(GSE47881_exprs))


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
        dplyr::select(timepoint, FP)  ->  df
        

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

GSE28422_res[!duplicated(GSE28422_res$logFC),]


GSE28422_res %>% 
        distinct(affy_hg_u133_plus_2, logFC, .keep_all = TRUE) 








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
        filter(age < 35, FP != "NB021")->  df

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

GSE47881_res[!duplicated(GSE47881_res$logFC),]


GSE28422_res %>% 
        distinct(affy_hg_u133_plus_2, logFC, .keep_all = TRUE) 





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

# GSE24235 # post is 24h after acute RT, so will skip this one
# GSE28422 # already done
# GSE28422 # already done
# GSE28998 # baseline and post samples are 4h post acute exercise, so will skip this one
# EMEXP740
# GSE106865 # done
# GSE8479 

#################################################################################################




#################################################################################################

#metamex 1: GSE24235

#################################################################################################



GSE24235 <- getGEO(GEO = "GSE24235", GSEMatrix = TRUE)

metadata_GSE24235 <- pData(phenoData(GSE24235[[1]]))

# post is 24h after acute RT, so will skip this one



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
        dplyr::select(2,3) -> df







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


x <- as.data.frame(illuminaHumanv2BeadIDENSEMBL)


GSE106865_res %>% 
        mutate(probe_id = rownames(GSE106865_res)) %>% 
        merge(.,x, by = "probe_id")


#################################################################################################

#metamex 6: GSE8479

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






































