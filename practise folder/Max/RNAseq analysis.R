############################################################################################
###                                                                                      ###
###  Analysis of RNAseq                                                                  ###
###                                                                                      ###
############################################################################################

# "Omic Association Studies with R and Bioconductor" by Gonzálex, Juan R. & Cáceres, Alejandro



BiocManager::install("tweeDEseqCountData")
BiocManager::install("tweeDEseq")

library(tweeDEseqCountData)
library(limma)
library(Biobase)
library(biomaRt)
library(tweeDEseq)

data("pickrell")

counts <- exprs(pickrell.eset)


# normalize by total counts per sample/column

lib.size <- colSums(counts)

NormByTotalNrReads <- sweep(counts, 2, FUN = "/", lib.size)


# load read length annotation dataframe from biomaRt package

data("annotEnsembl63")
head(annotEnsembl63)


# add annotation data to count data

genes.ok <- intersect(as.character(rownames(counts)), as.character(rownames(annotEnsembl63)))

geneAnnot <- annotEnsembl63[genes.ok,]
counts.ok <- counts[genes.ok,]

# check if count data and annotation data match

identical(rownames(geneAnnot), rownames(counts.ok))

# they match


###  normalize by RPKM

# get gene length

width <- geneAnnot$Length

NormByRPKM <- t(t(counts.ok / width * 1000) / colSums(counts.ok) * 1^6)


### normalize by TMM (trimmed mean of M-values) from the tweeDEseq package

NormByTMM <- normalizeCounts(counts.ok, method = "TMM")


## normalize by "CQN" in tweeDEseq package and cqn package

library(cqn)

NormByCQN <- normalizeCounts(counts.ok , method = "cqn", annot = geneAnnot[, c("Length", "GCcontent")])



# compare normalizations

# assumption: most genes are not differentially expressed between samples, thus correct normalization of samples should center data arround 0

MbyT <- log2(NormByTotalNrReads[,1] / NormByTotalNrReads[,2])

MbyRPKM <- log2(NormByRPKM[,1] / NormByRPKM[,2])


# draw histogram of read distributions 

par(mfrow = c(1,2))

hist(MbyT, xlab = "log2-ratio", main = "Total reads")
abline(v=0, col = "red")
hist(MbyRPKM, xlab = "log2-ratio", main = "RPKM")
abline(v=0, col = "red")

# histogram does not show us if the normalization worked at any expression level, therefore use MA-plot


library(edgeR)

par(mfrow = c(1,1))
maPlot(counts[,1], counts[,2], pch = 19, cex = .5, ylim = c(-8,8), 
       allCol = "darkgray", lowess = TRUE,
       xlab = expression(A == log[2] (sqrt(Sample1 %.% Sample2))),
       ylab = expression(M == log[2] (Sample1/Sample2)))
grid(col = "black")

# compare all 4 methods with non-nprmalized

par(mfrow = c(3,2))

maPlot(counts[,1], counts[,2], pch = 19, cex = .5, ylim = c(-8,8), 
       allCol = "darkgray", lowess = TRUE,
       xlab = expression(A == log[2] (sqrt(Sample1 %.% Sample2))),
       ylab = expression(M == log[2] (Sample1/Sample2)))
grid(col = "black")

maPlot(counts.ok[,1], counts.ok[,2], pch = 19, cex = .5, ylim = c(-8,8), 
       allCol = "darkgray", lowess = TRUE,
       xlab = expression(A == log[2] (sqrt(Sample1 %.% Sample2))),
       ylab = expression(M == log[2] (Sample1/Sample2)))
grid(col = "black")

maPlot(NormByCQN[,1], NormByCQN[,2], pch = 19, cex = .5, ylim = c(-8,8), 
       allCol = "darkgray", lowess = TRUE,
       xlab = expression(A == log[2] (sqrt(Sample1 %.% Sample2))),
       ylab = expression(M == log[2] (Sample1/Sample2)))
grid(col = "black")

NormByRPKM <- NormByRPKM[!is.na(NormByRPKM[,1]) & !is.na(NormByRPKM[,2]),]

maPlot(NormByRPKM[,1], NormByRPKM[,2], pch = 19, cex = .5, ylim = c(-8,8), 
       allCol = "darkgray", lowess = TRUE,
       xlab = expression(A == log[2] (sqrt(Sample1 %.% Sample2))),
       ylab = expression(M == log[2] (Sample1/Sample2)))
grid(col = "black")

maPlot(NormByTMM[,1], NormByTMM[,2], pch = 19, cex = .5, ylim = c(-8,8), 
       allCol = "darkgray", lowess = TRUE,
       xlab = expression(A == log[2] (sqrt(Sample1 %.% Sample2))),
       ylab = expression(M == log[2] (Sample1/Sample2)))
grid(col = "black")

maPlot(NormByTotalNrReads[,1], NormByTotalNrReads[,2], pch = 19, cex = .5, ylim = c(-8,8), 
       allCol = "darkgray", lowess = TRUE,
       xlab = expression(A == log[2] (sqrt(Sample1 %.% Sample2))),
       ylab = expression(M == log[2] (Sample1/Sample2)))
grid(col = "black")


# maPlot is interesting for comparison of normalizations 


### Filtering

# filter out genes with low read count (less than 5 per million) is common practice

NormByCQN.f <- filterCounts(NormByCQN, 
                            cpm.cutoff = 0.5,
                            n.samples.cutoff = 2,
                            mean.cpm.cutoff = 0)



# check filtering effect

dim(NormByCQN)
dim(NormByCQN.f)


### Run PCA on samples


pca.out <- prcomp(t(NormByCQN.f), scale. = TRUE)


plot(pca.out$x[,2:3])

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


# check for outliers

# Compute Mahalanobis distances for the first two principal components
distances <- mahalanobis(pca.out$x[,1:2], colMeans(pca.out$x[,1:2]), cov(pca.out$x[,1:2]))

# Identify outliers as samples with a Mahalanobis distance greater than a certain threshold
# Here, I'm using the 97.5 percentile of the Chi-square distribution with 3 degrees of freedom as the threshold
outliers <- which(distances > qchisq(0.975, df = 3))

# Print the row names of the outliers
print(row.names(NormByCQN.f)[outliers])

# Plot the first two principal components, highlighting the outliers
library(ggplot2)
library(ggrepel)
ggplot(data.frame(pca.out$x), aes(x = PC1, y = PC2)) +
        geom_point() +
        geom_text_repel(data = data.frame(pca.out$x)[outliers,], aes(label = rownames(data.frame(pca.out$x))[outliers])) +
        labs(x = "PC1", y = "PC2", title = "First Two Principal Components with Outliers Highlighted")



# filter out samples more than 3 degrees of freedom from the mahalanobis point

NormByCQN.ff <- NormByCQN.f[,!(colnames(NormByCQN.f) %in% names(outliers))]
dim(NormByCQN.ff)
dim(NormByCQN.f)




### differential gene expression

# RNAseq data is basically counts of reads, and as such is discrete and not continuous. Therefore we have to employ "negative binomial dispersion" instead of linear models.
# could we use logistic regression?


# get sex data

pheno <- pData(pickrell.eset)

# estimate common dispursion between genders

# remove filtered columns/outliers from sex df

pheno <- pheno[!rownames(pheno) %in% names(outliers),]

sex <- pheno$gender

d <- DGEList(counts = NormByCQN.ff, group = as.factor(sex))

d <- estimateCommonDisp(d)

# common dispersion

d$common.dispersion

sqrt(d$common.dispersion)


# estimate diffreence between each gene

d <- estimateTagwiseDisp(d)

names(d)


# use exact test to estimate difference between common and tag wise dispersion between genes ~ sex

resEdgeR.common <- exactTest(d, pair = c("female" , "male"), dispersion = "common")


resEdgeR.tagwise <- exactTest(d, pair = c("female" , "male"), dispersion = "tagwise")


# get results

topTags(resEdgeR.common)

topTags(resEdgeR.tagwise)

library(tidyverse)
as.data.frame(resEdgeR.common$table) %>% 
        mutate(FDR = p.adjust(PValue, method = "fdr")) %>% 
        filter(FDR <=0.05) %>% 
        arrange(FDR) %>%  rownames() -> common

as.data.frame(resEdgeR.tagwise$table) %>% 
        mutate(FDR = p.adjust(PValue, method = "fdr")) %>% 
        filter(FDR <=0.05) %>% 
        arrange(FDR) %>%  rownames() -> tagwise      

# check overlap between common and tagwise genes

df <- list(common = common, 
           tagwise = tagwise)

library(ggvenn)

ggvenn(df)


# add gene symbol name to differential gene expressions 


# get ensembl and gene names

ensembl <- useEnsembl(biomart = "genes")


ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')

gene_ids <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.con) 

colnames(gene_ids) <- c("ENTREZID", "GENE_NAME")
        


# merge dataframes

as.data.frame(resEdgeR.common$table) %>% 
        mutate(FDR = p.adjust(PValue, method = "fdr")) %>% 
        filter(FDR <=0.05) %>%
        rownames_to_column(var = "ENTREZID") %>% 
        merge(., gene_ids, by = "ENTREZID") %>% 
        arrange(-abs(logFC)) 
        

as.data.frame(resEdgeR.tagwise$table) %>% 
        mutate(FDR = p.adjust(PValue, method = "fdr")) %>% 
        filter(FDR <=0.05) %>%
        rownames_to_column(var = "ENTREZID") %>% 
        merge(., gene_ids, by = "ENTREZID") %>% 
        arrange(-abs(logFC)) 
