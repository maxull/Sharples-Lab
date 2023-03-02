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

library(wateRmelon); library(methylumi);library(FDb.InfiniumMethylation.hg19);library(minfi); library(tidyverse);library(ggplot2)
library(pathview); library(gage); library(gageData); library(ChAMP); library(org.Hs.eg.db); library(AnnotationDbi);
library(reshape2)


BiocManager::install("ChAMP")

library(R.utils); library(plyr);library(minfi); library(tidyverse);library(ggplot2); library(pathview); library(gage); library(gageData); library(ChAMP); library(org.Hs.eg.db); library(AnnotationDbi)

library(ENmix); library(FedData); library(cowplot)
library(limma)
library(RColorBrewer)

library(pheatmap)
library(viridis)

library(devtools)
install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest") 
install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")

library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)  # newest illumina annotation data


### set working directory to filepath that contains the data of interest
### this folder should contain the raw red and green signal .idat files and a CSV file of the targets

### save directory for later! 

setwd("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/")

mypath <- getwd()

#############################################################
### unzipping of zip files when downloading from GEO database (change file paths)


gunzip("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/GC-AS-10179.tar.gz", remove=FALSE)

untar("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/GC-AS-10179.tar", exdir = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/raw idat/")



#############################################################
#############################################################

### load data from idat files and filter

# load raw data 

myData <- champ.import(directory = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/raw idat/", 
                       arraytype = "EPIC")

myLoad <- champ.load(directory = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/raw idat/", 
                     arraytype = "EPIC", method = "minfi")

# add sample group to pd sheets

myLoad$pd$Sample_Group <- substr_right(myLoad$pd$Sample_Name,2)
myData$pd$Sample_Group <- substr_right(myLoad$pd$Sample_Name,2)

# identify principal components

### PCA analysis

# https://bioinfo4all.wordpress.com/2021/01/31/tutorial-6-how-to-do-principal-component-analysis-pca-in-r/

install.packages(c("factoextra", "FactoMineR"))

library(factoextra)
library(FactoMineR)


pca.data <- PCA(t(myData$beta), scale.unit = FALSE, graph = FALSE)


fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 100))

fviz_pca_ind(pca.data, addEllipses = FALSE,
             col.ind = "cos2", 
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
             repel = TRUE)


### first three principal components explain 57.6 % of variance in Dataset, therefore normalize with 3 principal components
### Normalize data

# loag rgSet to do functional normalization

myLoad <- champ.load(directory = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/raw idat/", 
                     arraytype = "EPIC", method = "minfi")

# change annotation to new annotation file


myLoad$rgSet@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")

# normalize


myNorm_fun <- preprocessFunnorm(myLoad$rgSet, 
                                nPCs=3, 
                                sex = NULL, 
                                bgCorr = TRUE,
                                dyeCorr = TRUE, 
                                keepCN = FALSE, 
                                ratioConvert = TRUE,
                                verbose = TRUE)



# notes from running processFunnorm
####################################################################

# > myNorm_fun <- preprocessFunnorm(myLoad$rgSet, 
#                                   +                                 nPCs=3, 
#                                   +                                 sex = NULL, 
#                                   +                                 bgCorr = TRUE,
#                                   +                                 dyeCorr = TRUE, 
#                                   +                                 keepCN = FALSE, 
#                                   +                                 ratioConvert = TRUE,
#                                   +                                 verbose = TRUE)
# [preprocessFunnorm] Background and dye bias correction with noob
# Loading required package: IlluminaHumanMethylationEPICanno.ilm10b4.hg19
# [preprocessFunnorm] Mapping to genome
# [preprocessFunnorm] Quantile extraction
# [preprocessFunnorm] Normalization
# Warning message:
#         In .getSex(CN = CN, xIndex = xIndex, yIndex = yIndex, cutoff = cutoff) :
#         An inconsistency was encountered while determining sex. One possibility is that only one sex is present. We recommend further checks, for example with the plotSex function.
#  
# 






####################################################################



####################################################################

# filter normalized values

####################################################################

# add  the normalized beta values to myData

myImport = myData

myImport$beta <- getBeta(myNorm_fun)

length(rownames(myImport$beta)) # all cpgs are still there, move on with filtering

flt_beta <- champ.filter(beta = myImport$beta,
                         pd = myImport$pd,
                         detP = myImport$detP[rownames(myImport$beta),],
                         beadcount = myImport$beadcount[rownames(myImport$beta),],
                         Meth = myImport$Meth[rownames(myImport$beta),],
                         UnMeth = myImport$UnMeth[rownames(myImport$beta),],
                         arraytype = "EPIC")

# notes from running processFunnorm
####################################################################

# > flt_beta <- champ.filter(beta = myImport$beta,
#                            +                          pd = myImport$pd,
#                            +                          detP = myImport$detP[rownames(myImport$beta),],
#                            +                          beadcount = myImport$beadcount[rownames(myImport$beta),],
#                            +                          Meth = myImport$Meth[rownames(myImport$beta),],
#                            +                          UnMeth = myImport$UnMeth[rownames(myImport$beta),],
#                            +                          arraytype = "EPIC")
# [===========================]
# [<<<< ChAMP.FILTER START >>>>>]
#         
#         In New version ChAMP, champ.filter() function has been set to do filtering on the result of champ.import(). You can use champ.import() + champ.filter() to do Data Loading, or set "method" parameter in champ.load() as "ChAMP" to get the same effect.
# 
# This function is provided for user need to do filtering on some beta (or M) matrix, which contained most filtering system in champ.load except beadcount. User need to input beta matrix, pd file themselves. If you want to do filterintg on detP matrix and Bead Count, you also need to input a detected P matrix and Bead Count information.
# 
# Note that if you want to filter more data matrix, say beta, M, intensity... please make sure they have exactly the same rownames and colnames.
# 
# 
# [ Section 1:  Check Input Start ]
# You have inputed beta,Meth,UnMeth for Analysis.
# 
# pd file provided, checking if it's in accord with Data Matrix...
#     pd file check success.
# 
#   Parameter filterDetP is TRUE, checking if detP in accord with Data Matrix...
#     detP check success.
# 
#   Parameter filterBeads is TRUE, checking if beadcount in accord with Data Matrix...
#     beadcount check success.
# 
#   parameter autoimpute is TRUE. Checking if the conditions are fulfilled...
#     !!! ProbeCutoff is 0, which means you have no needs to do imputation. autoimpute has been reset FALSE.
# 
#   Checking Finished :filterDetP,filterBeads,filterMultiHit,filterSNPs,filterNoCG,filterXY would be done on beta,Meth,UnMeth.
#   You also provided :detP,beadcount .
# [ Section 1: Check Input Done ]
# 
# 
# [ Section 2: Filtering Start >>
# 
#   Filtering Detect P value Start
#     The fraction of failed positions per sample
#     You may need to delete samples with high proportion of failed probes:
# 
#      Failed CpG Fraction.
# 4BH          0.0008846706
# 2BH          0.0006617706
# 5PM          0.0005624472
# 5BM          0.0006479115
# 12BM         0.0005636022
# 1PM          0.0006687001
# 5BH          0.0008927551
# 8PH          0.0018998474
# 2BM          0.0005820809
# 12BH         0.0005231799
# 12PM         0.0006028695
# 7PM          0.0004504198
# 1PH          0.0004561944
# 12PH         0.0005497431
# 7PH          0.0005185602
# 6PH          0.0007714882
# 2PM          0.0006144187
# 8BM          0.0007587841
# 7BH          0.0004873773
# 4BM          0.0004712084
# 4PM          0.0004446451
# 6BH          0.0004157721
# 8PM          0.0012496261
# 1BM          0.0011745561
# 5PH          0.0007888120
# 7BM          0.0005497431
# 1BH          0.0005301094
# 6BM          0.0004596591
# 4PH          0.0005000814
# 6PM          0.0006536861
# 8BH          0.0011757111
# 2PH          0.0011861054
# 
#     Filtering probes with a detection p-value above 0.01.
#     Removing 5278 probes.
#     If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples
# 
#   Filtering BeadCount Start
#     Filtering probes with a beadcount <3 in at least 5% of samples.
#     Removing 15736 probes
# 
#   Filtering NoCG Start
#     Only Keep CpGs, removing 2896 probes from the analysis.
# 
#   Filtering SNPs Start
#     Using general EPIC SNP list for filtering.
#     Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
# Removing 95783 probes from the analysis.
# 
# Filtering MultiHit Start
# Filtering probes that align to multiple locations as identified in Nordlund et al
# Removing 11 probes from the analysis.
# 
# Filtering XY Start
# Filtering probes located on X,Y chromosome, removing 16369 probes from the analysis.
# 
# Updating PD file
# 
# Fixing Outliers Start
# Replacing all value smaller/equal to 0 with smallest positive value.
# Replacing all value greater/equal to 1 with largest value below 1..
# [ Section 2: Filtering Done ]
# 
# All filterings are Done, now you have 729786 probes and 32 samples.
# 
# [<<<<< ChAMP.FILTER END >>>>>>]
# [===========================]
# [You may want to process champ.QC() next.]




####################################################################

# visualize normalization and distribution + QC

####################################################################


# mean and sd raw data

long_beta <- myData$beta %>% 
        as.data.frame() %>% 
        pivot_longer(names_to = "FP", values_to = "beta", cols = 1:32)

beta_mean <- long_beta %>% 
        mutate(FP = ifelse(nchar(FP) < 4, paste(0, FP, sep = ""), FP)) %>%
        group_by(FP) %>% 
        summarise(m = mean(beta),
                  s = sd(beta))



reso <- 600
length <- 3.25*reso/72

png(file="C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/Figures/raw beta.png",
    units="in",res=reso,height=length,width=length*2)


long_beta %>% 
        mutate(FP = ifelse(nchar(FP) < 4, paste(0, FP, sep = ""), FP)) %>% 
        mutate(NR = substr(FP,1,2)) %>% 
        ggplot(aes(x = FP, y = beta, color = NR))+
        geom_violin()+
        geom_boxplot(color = "black", width = 0.2)+
        theme(axis.text.x = element_text(angle = -90))+
        labs(title = "Raw beta values")

dev.off()

long_beta %>% 
        group_by(FP) %>% 
        summarise(m = mean(beta),
                  s = sd(beta))%>% 
        ggplot(aes(x = FP, y = m))+
        geom_point()+
        geom_errorbar(aes(ymin = m-s, ymax = m+s), width = 0.1)


# density plots of normalized and pre-normalized values, and after filtering

par(mfrow = c(3, 1))

densityPlot(myData$beta, sampGroups=myData$pd$Sample_Group,
            main="Pre-Normalized", legend=TRUE)

densityPlot(myImport$beta, sampGroups=myImport$pd$Sample_Group,
            main="Normalized", legend=TRUE)

densityPlot(flt_beta$beta, sampGroups=flt_beta$pd$Sample_Group,
            main="Normalized+filtered", legend=TRUE)

dev.off()


champ.QC(beta = flt_beta$beta)


### visualize normalization and save density plots


reso <- 600
length <- 3.25*reso/72

png(file="C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/Figures/Normalized density plot.png",
    units="in",res=reso,height=length,width=length*2)
densityPlot(myNorm@assays@data@listData[["Beta"]], sampGroups=myLoad$pd$Sample_Group,
            main="Normalized", legend=TRUE)
dev.off()

png(file="C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/Figures/Non-Normalized density plot.png",
    units="in",res=reso,height=length,width=length*2)
densityPlot(myLoad$beta, sampGroups=myLoad$pd$Sample_Group,
            main="Not-Normalized", legend=TRUE)
dev.off()


### run quality control



pal <- brewer.pal(8,"Dark2")
plotMDS(B2M(flt_beta$beta), top=1000, gene.selection="common",
        col=pal[factor(flt_beta$pd$Sample_Group)])



### PCA analysis

# https://bioinfo4all.wordpress.com/2021/01/31/tutorial-6-how-to-do-principal-component-analysis-pca-in-r/

install.packages(c("factoextra", "FactoMineR"))

library(factoextra)
library(FactoMineR)

flt_beta$beta

pca.data <- PCA(t(myData$beta), scale.unit = FALSE, graph = FALSE)


fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 100))

fviz_pca_ind(pca.data, addEllipses = FALSE,
             col.ind = "cos2", 
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
             repel = TRUE)


champ.QC()      # default QC on myLoad

QC.GUI()



getwd()

setwd("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/")


### save temporary files

saveRDS(EPIC.manifest.hg19, file = "EPIC.manifest.hg19.RDATA")

saveRDS(multi.hit, file = "multi.hit.RDATA")

saveRDS(myLoad, file = "myLoad.RDATA")

saveRDS(probe.features, file = "probe.features.RDATA")

saveRDS(myNorm_fun, file = "myNorm_fun.RDATA")

saveRDS(myData, file = "myData.RDATA")

saveRDS(flt_beta, file = "flt_beta.RDATA")


### reload temporary files (not all are neaded for downstream analysis!)

EPIC.manifest.hg19 <- readRDS("EPIC.manifest.hg19.RDATA")

multi.hit <- readRDS("multi.hit.RDATA")

myLoad <- readRDS("myLoad.RDATA")

probe.features <- readRDS("probe.features.RDATA")

myNorm <- readRDS("myNorm.RDATA")

myData <- readRDS("myData.RDATA")

flt_beta <- readRDS("flt_beta.RDATA")

### set working directory back to original


setwd(mypath)


################################################################################################################################
################################################################################################################################
################################################################################################################################

### identify DMPs from M vals

####################################################################################




myDMP_BH_BM <- champ.DMP(beta = B2M(flt_beta$beta),
                   pheno = flt_beta$pd$Sample_Group,
                   compare.group = c("BH", "BM"),
                   arraytype = "EPIC")                   ### significant DMPs with BH adjusted p val of 0.05

myDMP_PH_PM <- champ.DMP(beta = B2M(flt_beta$beta),
                         pheno = flt_beta$pd$Sample_Group,
                         compare.group = c("PH", "PM"),
                         arraytype = "EPIC")                   ### significant DMPs with BH adjusted p val of 0.05

myDMP_BM_PM <- champ.DMP(beta = B2M(flt_beta$beta),
                         pheno = flt_beta$pd$Sample_Group,
                         compare.group = c("BM", "PM"),
                         adjust.method = "none",
                         arraytype = "EPIC")                   ### significant DMPs with un-adjusted p val of 0.05

myDMP_BH_PH <- champ.DMP(beta = B2M(flt_beta$beta),
                         pheno = flt_beta$pd$Sample_Group,
                         compare.group = c("BH", "PH"),
                         adjust.method = "none",
                         arraytype = "EPIC")                   ### significant DMPs with un-adjusted p val of 0.05


# code to check DMPs with differenct p values
myDMP$BH_to_BM %>% 
        as.data.frame() %>% 
        filter(adj.P.Val < 0.001) %>% 
        rownames() %>% 
        length()


###

# add newer anntoaton data
library(dplyr)
data(Other)
data("Islands.UCSC")

anno <- as.data.frame(Other[,c("UCSC_RefGene_Name", "Regulatory_Feature_Group")]) %>% 
        rownames_to_column(var = "cpg")

Islands <- Islands.UCSC[, 2] %>% 
        rownames()


Islands <-data.frame(cpg = rownames(Islands.UCSC),
                     Relation_to_Island = Islands.UCSC[,2])


anno <- merge(anno, Islands, by = "cpg")



anno <- anno %>% 
        mutate(UCSC_RefGene_Name = ifelse(UCSC_RefGene_Name == "", paste("NA"), paste(UCSC_RefGene_Name)))



library(stringr)


anno <- anno %>% 
        mutate(split = str_split(UCSC_RefGene_Name, ";")) %>% # split
        mutate(split = map(.$split, ~ unique(.x))) %>% # drop duplicates
        mutate(UCSC_RefGene_Name = map_chr(.$split, ~paste(.x, collapse = ";"))) # recombine

anno <- anno[,1:4]

saveRDS(anno, "anno.RDATA")

anno <- readRDS("anno.RDATA")


# BH adjusted p val of 0.05 gave no significant DMPs between BH to PH, and BM to PM



### save DMPs


saveRDS(myDMP_BH_BM, "myDMP_BH_BM.RDATA")

saveRDS(myDMP_PH_PM, "myDMP_PH_PM.RDATA")

saveRDS(myDMP_BM_PM, "myDMP_BM_PM.RDATA")

saveRDS(myDMP_BH_PH, "myDMP_BH_PH.RDATA")

# reload DMPs

myDMP_BH_BM <- readRDS("myDMP_BH_BM.RDATA")

myDMP_PH_PM <- readRDS("myDMP_PH_PM.RDATA")

myDMP_BM_PM <- readRDS("myDMP_BM_PM.RDATA")

myDMP_BH_PH <- readRDS("myDMP_BH_PH.RDATA")

# done to here

####

# to do
# - DMR
# - GSEA
# - global test etc.



################################################################################################################################


# create heatmap of dmps

###################################################################################

# add annotation of condition

dfh <- data.frame(sample = as.character(colnames(flt_beta$beta)), condition = "Timepoint") %>% 
        column_to_rownames("sample") %>% 
        mutate(condition = substr_right(rownames(.),2))


# filter bvals for only conditions of interest

dfh_1 <- dfh %>% 
        filter(condition %in% c("BH", "BM"))

dmp_list <- rownames(myDMP_BH_BM$BH_to_BM)

b_vals <- flt_beta$beta

b_vals_1 <- b_vals[dmp_list, rownames(dfh_1)]
        

# isolate DMPs within islands

annEPIC <- read_excel("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/Annotation file Illumina/infinium-methylationepic-v-1-0-b5-manifest-file.xlsx")

ann_df <- annEPIC@listData %>% 
        as.data.frame()


Island <- ann_df %>% 
        select(Name, Relation_to_Island) %>% 
        filter(Relation_to_Island == "Island") %>% 
        column_to_rownames(var = "Name")

# check that island cpg exists in b_vals list


b_vals_isl <- b_vals[rownames(b_vals) %in% c(rownames(Island)),] # keep cpg rows that map to cpg islands

# filter DMPs BH to BM

b_vals_3 <- b_vals_isl[, rownames(dfh_1)]

# 160000 + rows is still to much, lets keep highest fold change DMPs

dmp_isl <- myDMP$BH_to_BM %>%
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        filter(cgi == "island") %>% 
        arrange(logFC)

b_vals_isl <- b_vals[rownames(b_vals) %in% c(dmp_isl$cpg),rownames(dfh_1)]

pheatmap(t(b_vals_isl), annotation_row = dfh_1, cutree_rows = 2, 
         show_colnames = FALSE, annotation_names_row = FALSE, 
         color = viridis(n = 100), scale = "none")





# heatmap all timepoints and conditions

dfh_2 <- dfh %>% 
        filter(condition %in% c("PH", "PM"))

# change direction of post dmps, so it matches with baseline

dmp_isl_2 <- myDMP$PM_to_PH %>%
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        filter(cgi == "island") %>% 
        mutate(logFC = logFC*-1) %>% 
        arrange(logFC)

# plot post DMPs


b_vals_isl_2 <- b_vals[rownames(b_vals) %in% c(dmp_isl_2$cpg),rownames(dfh_2)]

pheatmap(t(b_vals_isl_2), annotation_row = dfh_2, cutree_rows = 2, 
         show_colnames = FALSE, annotation_names_row = FALSE, 
         color = viridis(n = 100), scale = "none", )


# plot baseline to post change

# homogenate

dmp_isl_change <- myDMP_BH_PH$BH_to_PH %>%
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>%  
        arrange(logFC)

dfh_3 <- dfh %>% 
        filter(condition %in% c("BH", "PH"))

b_vals_change <- b_vals[rownames(b_vals) %in% c(dmp_isl_change$cpg),rownames(dfh_3)]

pheatmap(t(b_vals_change), annotation_row = dfh_3, cutree_rows = 2, 
         show_colnames = FALSE, annotation_names_row = FALSE, 
         color = viridis(n = 100), scale = "none")

# myonuclei

dmp_isl_change2 <- myDMP2$PM_to_BM %>%
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        mutate(logFC = logFC*-1) %>%
        arrange(logFC)

dfh_4 <- dfh %>% 
        filter(condition %in% c("BM", "PM"))

b_vals_change2 <- b_vals[rownames(b_vals) %in% c(dmp_isl_change2$cpg),rownames(dfh_4)]

pheatmap(t(b_vals_change2), annotation_row = dfh_4, cutree_rows = 2, 
         show_colnames = FALSE, annotation_names_row = FALSE, 
         color = viridis(n = 100), scale = "none")


#################################################################################

### export csv files of DMPs

#################################################################################



myDMP_BH_PH$BH_to_PH %>%
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        arrange(logFC) %>% 
        merge(., anno, by = "cpg") %>% 
        write.csv("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/DMPs/BH_to_PH_unadj. P 0.05.csv", row.names = FALSE)


myDMP_BH_BM$BH_to_BM %>% 
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        arrange(logFC) %>%
        merge(., anno, by = "cpg") %>% 
        write.csv("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/DMPs/BH_to_BM_FDR P 0.05.csv", row.names = FALSE)


myDMP_PH_PM$PH_to_PM %>% 
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        arrange(logFC) %>%
        merge(., anno, by = "cpg") %>%
        write.csv("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/DMPs/PH_to_PM_FDR P 0.05.csv", row.names = FALSE)

myDMP_BM_PM$BM_to_PM %>% 
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        arrange(logFC) %>%
        merge(., anno, by = "cpg") %>%
        write.csv("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/DMPs/BM_to_PM_unadj. P 0.05.csv", row.names = FALSE)


#################################################################################

# cell pop correction

#################################################################################

dfh_2<- dfh %>% 
        filter(condition %in% c("BH", "PH"))

b_vals_2 <- b_vals[, rownames(dfh_2)]

myRefbase <- champ.refbase(beta = b_vals_2,                      ### returns b-vals adjusted for 5 main cellpopulations identified
                           arraytype = "EPIC")

myRefbase$CorrectedBeta

densityPlot(myNorm@assays@data@listData[["Beta"]], sampGroups=myLoad$pd$Sample_Group,
            main="Normalized", legend=TRUE)

densityPlot(myRefbase$CorrectedBeta, sampGroups=dfh_2$condition,
            main="Normalized", legend=TRUE)





getwd()

setwd(mypath)










#################################################################################

### identify DMRs

#################################################################################

myDMR <- champ.DMR(arraytype = "EPIC",
                   compare.group = c("BH", "BM"))          #no significant DMPs with BH adjusted p val of 0.05

myDMR_BH_PH <- champ.DMR(beta = B2M(flt_beta$beta),
                   pheno = flt_beta$pd$Sample_Group,
                   compare.group = c("BH", "PH"),
                   arraytype = "EPIC",
                   method = "Bumphunter",
                   adjPvalDmr = "none",
                   cores = 4)

dmr <- as.data.frame(myDMR_BH_PH$BumphunterDMR)


### DMR did not work


# global test

################################################################################


BiocManager::install("globaltest")

library(globaltest)

# prepare subsets

unique_cpg <- list()

unique_gene <- unique(anno$UCSC_RefGene_Name)



for (i in 1:length(unique_gene)) {
        
        index <- which(anno$UCSC_RefGene_Name == unique_gene[i])
        
        z <- anno[index,1]
        
        unique_cpg[[(as.character(unique_gene[i]))]] <- z
        
        
        print(i)
}


# remove "NA" from unique_cpg list

unique_cpg2 <- unique_cpg[-2]

cpgs <- rownames(flt_beta$beta)
x <- intersect(unique_cpg[[1]], cpgs)



# remove cpg that do not exist in filtered b vals

for (i in 1:length(unique_cpg)) {
        unique_cpg2[[i]] <- intersect(unique_cpg[[i]], cpgs)
        print(i)
}

subsets <- unique_cpg2
saveRDS(subsets, "subsets.RDATA")

unique_cpg2

subsets <- list()


for (i in 1:lenth(unique_cpg)) {
        unique_cpg
        
        print(i)
        
}



for (i in 1:length(uniquecg$`unique(CGtoGENE2$ENTREZ)`)) {
        index <- which(array2$ENTREZ == uniquecg[i,])
        z <- array2[index,]
        mygenes[[(as.character(uniquecg[i,]))]] <- z$rowname
        print(i)
        
}













saveRDS(unique_cpg, "unique_cpg.RDATA")


unique_cpg <- readRDS("unique_cpg.RDATA")


### isolate M values of homogenate and myonuclei in separate meatrix

m_vals <- B2M(flt_beta$beta)

# list of homogenate samples

homogenate_samples <- flt_beta$pd %>% 
        filter(Sample_Group %in% c("BH", "PH")) %>% 
        as.data.frame() %>% 
        dplyr::select(1,4) 

myonuclei_samples <- flt_beta$pd %>% 
        filter(Sample_Group %in% c("BM", "PM")) %>% 
        as.data.frame() %>% 
        dplyr::select(1,4)

m_vals_homo <- m_vals[,homogenate_samples$Sample_Name]

m_vals_myo <- m_vals[,myonuclei_samples$Sample_Name]


# create df with baseline = 0 and post = 1

condition_homo <- flt_beta$pd %>% 
        filter(Sample_Group %in% c("BH", "PH")) %>%
        mutate(condition = substr_right(Sample_Group, 2),
               condition = substr(condition,1,1),
               condition = ifelse(condition == "B", 0, 1)) %>% 
        dplyr::select(1,8)

condition_myo <- flt_beta$pd %>% 
        filter(Sample_Group %in% c("BM", "PM")) %>%
        mutate(condition = substr_right(Sample_Group, 2),
               condition = substr(condition,1,1),
               condition = ifelse(condition == "B", 0, 1)) %>% 
        dplyr::select(1,8)

rownames(condition_homo) <- condition_homo[,1]

condition_homo <- condition_homo %>% 
        dplyr::select(2)

rownames(condition_myo) <- condition_myo[,1]

condition_myo <- condition_myo %>% 
        dplyr::select(2)

# ready for global test

# m_vals = m_vals_homo and m_vals_myo
# subsets = unique_cpg
# pheno = condition_homo and condition_myo


library(globaltest)
gt_homo <- gt(response = condition_homo$condition, alternative = t(m_vals_homo), subsets = unique_cpg, model = "logistic", direction = FALSE, permutations = 0)

################################################################################################################################

### identify GSEA



### default test
myGSEA <- champ.GSEA(arraytype = "EPIC",
                     DMP = myDMP2[["BH_to_PH"]],
                     cores = 4)


GSEA <- as.data.frame(myGSEA$DMP)


### gOmeth method
myGSEA2 <- champ.GSEA(arraytype = "EPIC",
                      DMP = myDMP[["BH_to_BM"]],
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


myRefbase <- champ.refbase(beta = b_Vals,                      ### returns b-vals adjusted for 5 main cellpopulations identified
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



################################################################################################################################





