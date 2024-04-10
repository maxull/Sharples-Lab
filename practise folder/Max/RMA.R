# RMA


library(cowplot)
library(ggrepel)
library(scales)
library(tidyverse)
library(ggplot2)
library(pathview)
library(viridis)
library(ENmix)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(sva)

counts <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/RMA_RNAseq_counts.RDATA")


# annotate Ensebl to gene name
counts <- counts %>% 
        mutate(Ensembl = gsub("\\..*.", "", Gene.ID))

counts$Gene <-  mapIds(org.Hs.eg.db, keys = counts$Ensembl, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")


########################################################################

### PCA

################################################################

colnames(counts)

# PCA 1 correlates with cell population

pca.out <- prcomp(t(counts[,2:37]), scale. = FALSE)

plot(pca.out$x[,1:2])

loadings <- pca.out$rotation[,1]



###################################################################
### DEseq2 analysis
###############################################################

data <- counts[,2:37]
# Convert data frame to DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = data,
                              design = ~ Participant + Timepoint)

# Estimate size factors
dds <- estimateSizeFactors(dds)

# Subset genes with average count > 1
dat <- counts(dds, normalized = TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx, ]

# Specify full and null model
mod <- model.matrix(~ Participant + Timepoint, colData(dds))
mod0 <- model.matrix(~ Participant, colData(dds))

# Determine number of SVs needed
n.sv <- num.sv(dat, mod, method = "be")  # Test with 'be' and 'leek' method, choose the lowest n
n.sv  # Output the number of SVs

# Perform SVAseq
svseq <- svaseq(dat, mod, mod0, n.sv = n.sv)

# Add new SVs to design
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
design(ddssva) <- ~ Participant + SV1 + Timepoint

# Perform differential expression analysis
ddssva <- DESeq(ddssva)

# DESeq QC: Dispersion plot
par(mfrow=c(1,1))
plotDispEsts(ddssva)



















