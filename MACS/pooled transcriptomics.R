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
library(ggplot2)
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
library(ggplot2)
library(ggrepel)
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


# there are less annotated "gene name", so will use ensmbl id for now

# remove .CEL from colnames

names(pm) <- gsub(pattern = ".CEL", replacement = "", names(pm))




# filter GSE28422 for low read counts

GSE28422.f <- filterCounts(pm,
                           cpm.cutoff = 0.5,
                           n.samples.cutoff = 2,
                           mean.cpm.cutoff = 0)


# normalize 

scaled <- scale(GSE28422.f, scale = FALSE)

plotDensities(object = scaled, legend = FALSE)




# run PCA

pca.out <- prcomp(t(GSE28422.f), scale. = TRUE)

# check proportion of variability explained by pca 1-50

# Extract the proportion of variance explained by each principal component
var_explained <- pca.out$sdev^2 / sum(pca.out$sdev^2)

# Create a data frame for plotting
df <- data.frame(Component = 1:length(var_explained), Variance = var_explained)

# Keep only the first 50 principal components
df <- df[1:50,]

# Plot
library(ggplot2)
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
print(colnames(GSE28422.f)[outliers])

# Plot the first two principal components, highlighting the outliers
library(ggplot2)
library(ggrepel)
ggplot(data.frame(pca.out$x), aes(x = PC1, y = PC2)) +
        geom_point() +
        geom_text_repel(data = data.frame(pca.out$x)[outliers,], aes(label = rownames(data.frame(pca.out$x))[outliers])) +
        labs(x = "PC1", y = "PC2", title = "First Two Principal Components with Outliers Highlighted")


GSE28422.ff <- GSE28422.f[,!(colnames(NormByCQN.f) %in% names(outliers))]



# study 2 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47881
# Phillips et al. (2013)







