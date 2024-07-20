###
###
###             DNA methyaltion deconvolution with EpiSCORE
###
###             DID NOT WORK! 






# Install packages

library(devtools)
devtools::install_github("immunogenomics/presto")

options(timeout = 5*60)  # extend timeout to 5 min
devtools::install_github("aet21/EpiSCORE")

BiocManager::install("EpiDISH")


library(tidyverse)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MASS)
library(EpiDISH)
library(presto)
library(EpiSCORE)
library(Seurat)


##########################################################################################
##########      Load data                       ##########################################
##########################################################################################

# DNA methylation profile
beta <- readRDS("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/Skole/M.Sc/Master 21-22/Master/DATA/Epigenetics/beta.RDATA")

# single cell RNAseq
# singlecell_counts <- readRDS("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/Ph.D/deconvolution/data/singlecell_counts.RDATA")
# singlecell_metadata <- readRDS("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/Ph.D/deconvolution/data/singlecell_metadata.RDATA")
# 
# singlefiber_counts <- readRDS("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/Ph.D/deconvolution/data/singlefiber_counts.RDATA")
# singlefiber_metadata <- readRDS("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/Ph.D/deconvolution/data/singlefiber_metadata.RDATA")

# use these counts instead which contain myonuclei
load("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/Ph.D/deconvolution/data/integrated_seurat_only_turiel_nonPAD.RData")


##########################################################################################
##########      extract log_morm count data       ########################################
##########################################################################################


# # add 1 and log transform count matrix
# log_sc_counts <- log2(singlecell_counts+1)
# 
# log_sf_counts <- log2(singlefiber_counts+1)

sc_counts <- seurat_all@assays$RNA$data #counts is raw data and "data" is normalized data

# extract metadata

singlecell_metadata <- seurat_all@meta.data

tail(singlecell_metadata)

unique(singlecell_metadata$cluster)

##########################################################################################
##########      load marker genes from all cells in single cell data       ###############
##########################################################################################

# Used find all markers from Seurat package on single cell data
# markers <- read_csv("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/Ph.D/Singe_fiber_exercise/markers_incl_myonuclei.csv")


##########################################################################################
##########      build expression reference matrix.             ###########################
##########################################################################################

# example dataset

example <- load("/Users/maxullrich/Downloads/Skin_scRNAseq_construction.rda")



###
# constract index vecor and name vector for single cell counts dataset
###

# filter single cell data to only include "control"
sc_labes <- singlecell_metadata %>% 
        rownames_to_column(var = "Sample_ID") %>% 
        dplyr::select(Sample_ID, cluster) %>% summary

unique(sc_labes$cluster)

sc_labes$idx[grep("Muscle fibers", sc_labes$cluster)] <- 8
sc_labes$idx[grep("FAPs", sc_labes$cluster)] <- 1
sc_labes$idx[grep("Endothelial cells", sc_labes$cluster)] <- 2
sc_labes$idx[grep("Satellite cells", sc_labes$cluster)] <- 3
sc_labes$idx[grep("Myeloid cells", sc_labes$cluster)] <- 4
sc_labes$idx[grep("Smooth muscle/Pericytes", sc_labes$cluster)] <- 5
sc_labes$idx[grep("NK/T/B cells", sc_labes$cluster)] <- 6
sc_labes$idx[grep("Mast cells", sc_labes$cluster)] <- 7

summary(factor(sc_labes$idx))


###                                                             ###
# drop this part                                        ################################### here
###                                                             ###

sc_labes <- sc_labes %>% 
        filter(Sample_ID %in% colnames(sc_counts)) 



# filter and reorder sc counts so labels and counts match

sc_counts_flt <- sc_counts[,colnames(sc_counts) %in% sc_labes$Sample_ID]

colnames(sc_counts_flt[,1:10]) == sc_labes$Sample_ID %>% head(10)

# log normalize counts


log_sc_counts <- log2(sc_counts_flt+1)

dim(log_sc_counts)

log_sc_counts@x <- round(log_sc_counts@x, 3) 

###############################################################################################
expref.o <- ConstExpRef(exp.m = log_sc_counts, 
                        celltype.idx = sc_labes$idx, 
                        namesCellT.v = c("FAPs",
                             "Endothelial cells",
                             "Satellite cells",
                             "Myoeloid cells",
                             "Smooth muscle/Pericytes",
                             "NK/T/B cells",
                             "Mast cells",
                             "Muscle fibers"))                

# check reference matrix

dim(expref.o$ref$med)

head(expref.o$ref$med)

view((expref.o$ref$med))


##########################################################################################
##########      build DNAm reference matrix                    ###########################
##########################################################################################



refMscm2.m <- ImputeDNAmRef(expref.o$ref$med,db="SCM2",geneID="SYMBOL")

refMrmap.m <- ImputeDNAmRef(expref.o$ref$med,db="RMAP",geneID="SYMBOL")

# merge together

refMmg.m <- ConstMergedDNAmRef(refMscm2.m,refMrmap.m)

dim(refMmg.m)

# rownames are change to entrez gene IDs
# a final row of weights has been added which is basically a score of how sure are we that the promoter is hypermethylated


##########################################################################################
##########      validate DNAm reference matrix                 ###########################
##########################################################################################

data(dataExampleLung)

dim(avSIM.m)

# how many of our cell marker genes are abova threshold 0.4
plot(density(refMmg.m[,9]),lwd=2,xlab="Weight",main="")
abline(v=0.3,lwd=2,col="red")

paste("Number of selected genes=",length(which(refMmg.m[,5]>0.6)),sep="")

# use slightly lower threshold

estF.o <- wRPC(data=avSIM.m,ref=refMmg.m,useW=TRUE,wth=0.6,maxit=200)

cor(estF.o$estF,trueW.m)

# Estimation of proportions of endothelial cells and immune cells is well correlated, which are the only cell types present in both datasets



##########################################################################################
##########       Estimate cell population proportions in my data.     ####################
##########################################################################################

# average the promoter methylation in my samples
average_promoter.m <- constAvBetaTSS(beta.m = beta, type = "850k")

# estimate cell proportions
estf.m <- wRPC(data=average_promoter.m,ref=refMmg.m,useW=TRUE,wth=0.6,maxit=200)$estF


