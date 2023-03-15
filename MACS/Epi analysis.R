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

# Load the champ package
library(champ)

# Load the idat files from the directory
myLoad <- champ.load(directory = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/raw idat/", arraytype = "EPIC", method = "minfi")




# change annotation to new annotation file

# annotation headers explained: https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001568

myLoad$rgSet@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")

Illumina_anno <- read.csv("Annotation file Illumina/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7, header = TRUE)



count(Illumina_anno$UCSC_RefGene_Name)



# normalize


# preprocessFunnorm is a function that normalizes the data
# myLoad$rgSet is the data that is being normalized
# nPCs is the number of principal components to be used for normalization
# sex is the sex of the samples
# bgCorr is the background correction
# dyeCorr is the dye bias correction
# keepCN is the copy number correction
# ratioConvert is the ratio conversion
# verbose is the verbose output

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
myDMP_BH_BM$BH_to_BM %>% 
        as.data.frame() %>% 
        filter(adj.P.Val < 0.001) %>% 
        rownames() %>% 
        length()


###

# add newer anntoaton data
library(dplyr)
data(Other)
data("Islands.UCSC")

annotation <- as.data.frame(Other)

data(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)

as.data.frame(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
rm(manifest)
data("Locations")
manifest <- as.data.frame(Locations)


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



# run BMIQ normalization and identify DMPs across time

bmiq_norm <- champ.norm(beta = myLoad$beta, arraytype = "EPIC")

rownames(bmiq_norm) %>% 
        length()

myDMP_BM_PM_bmiq <- champ.DMP(beta = B2M(bmiq_norm),
                         pheno = flt_beta$pd$Sample_Group,
                         compare.group = c("BM", "PM"),
                        adjust.method = "none",
                         arraytype = "EPIC")                   ### significant DMPs with un-adjusted p val of 0.05

myDMP_BH_PH_bmiq <- champ.DMP(beta = B2M(bmiq_norm),
                         pheno = flt_beta$pd$Sample_Group,
                         compare.group = c("BH", "PH"),
                         adjust.method = "none",
                         arraytype = "EPIC")                   ### significant DMPs with un-adjusted p val of 0.05


########################################################################

# create horizontal barchart of DMPs in different positions

########################################################################


# start with baseline samples

baseline_dmps <- myDMP_BH_BM$BH_to_BM %>% 
        as.data.frame() %>%
        rownames_to_column(var = "cpg") %>% 
         dplyr::select(cpg, P.Value, deltaBeta)

# merge with newer annotation

baseline_dmps <- merge(baseline_dmps, anno, by = "cpg") %>% 
        as.data.frame() %>% 
        mutate("dir" = as.factor(ifelse(deltaBeta <0, "Negative_skew", "Positive_skew"))) 
        
dfh_1 <- dfh %>% 
        filter(condition %in% c("BH", "BM"))        

change_b_vals_baseline <- b_vals[,rownames(dfh_1)] %>% 
        as.data.frame() %>% 
        mutate(average_homo = (`1BH`+ `2BH`+ `4BH`+ `5BH`+ `6BH`+ `7BH`+ `8BH`+`12BH`)/8,
               average_myo = (`1BM`+ `2BM`+ `4BM`+ `5BM`+ `6BM`+ `7BM`+ `8BM`+`12BM`)/8,
               change = average_myo-average_homo) %>% 
        
        dplyr::select(change) %>% 
        rownames_to_column(var = "cpg")


n_count <- merge(baseline_dmps, change_b_vals_baseline, by.x = "cpg") %>% 
        filter(dir == "Negative_skew") %>% 
        mutate(Relation_to_Island = factor(Relation_to_Island, levels = c("S_Shelf","N_Shelf","Island", "S_Shore", "N_Shore", "OpenSea"))) %>% 
        summarise(count(Relation_to_Island)) %>% 
        as.data.frame()

p_count <- merge(baseline_dmps, change_b_vals_baseline, by.x = "cpg") %>% 
        filter(dir == "Positive_skew") %>%
        mutate(Relation_to_Island = factor(Relation_to_Island, levels = c("S_Shelf","N_Shelf","Island", "S_Shore", "N_Shore", "OpenSea"))) %>% 
        summarise(count(Relation_to_Island)) %>% 
        as.data.frame() 

df <- data.frame(
        category = n_count[6:1,]$x,
        neg_skew = (n_count[6:1,]$freq)*-1,
        pos_skew = p_count[6:1,]$freq)

p1 <- ggplot(df, aes(x = category)) +
        geom_bar(aes(y = neg_skew, fill = "neg_skew"), stat = "identity") +
        geom_bar(aes(y = pos_skew, fill = "pos_skew"), stat = "identity") +
        scale_fill_manual(values = c("#453781FF", "#DCE319FF")) + # set colors
        geom_label(aes(label = as.integer(neg_skew)*-1, x = category, y = as.integer(neg_skew)), size = 4, alpha = 0.8, nudge_y = -3000)+
        geom_label(aes(label = as.integer(pos_skew), x = category, y = as.integer(pos_skew)), size = 4, alpha = 0.8, nudge_y = 3000)+
        coord_flip()+
        theme_classic()+
        theme(legend.title = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_text(hjust = 0.58))+
        labs(y = "Number of DMPs",
             title = "Baseline DMPs in Homogenate vs. Myonuclei")


# find number of significant cpgs that were altered more than <0.01 B < 0.1 B < 0.2 B <


change_b_vals_baseline %>% 
        mutate(change = abs(change)) %>% # change all numbers to positive
        merge(baseline_dmps,.,by.x = "cpg") %>% 
        # find values in ranges
        mutate(skew = ifelse(change <0.01, "<0.01", 
                             ifelse(change >0.01 & change < 0.1, "0.01<0.1", 
                                    ifelse(change > 0.1 & change < 0.2, "0.1<0.2", 
                                           ifelse(change > 0.2 & change < 0.3, "0.2<0.3",">0.3"))))) %>% 
        summarise(count(skew)) -> df1
        
df2 <- data.frame(category = factor(df1[,1]$x, levels = c("<0.01", "0.01<0.1", "0.1<0.2", "0.2<0.3", ">0.3")),
                  count = df1[,1]$freq)

library(scales)
p3 <- ggplot(data = df2, aes(x = category))+
        geom_bar(aes(y = count), stat = "identity")+
        geom_text(aes(label = count, y = count), vjust = -0.5)+
        theme_classic()+
        labs(title = "number of DMPs by skewed B-values at baseline",
             x = "B-value skew")+
        theme(axis.title.y = element_blank())+
        scale_y_continuous(trans = "log2", expand  = c(0, 1), n.breaks = 9)


# repeat for post samples

post_dmps <- myDMP_PH_PM$PH_to_PM %>% 
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        dplyr::select(cpg, P.Value, deltaBeta)

# merge with newer annotation

post_dmps <- merge(post_dmps, anno, by = "cpg") %>% 
        as.data.frame() %>% 
        mutate("dir" = as.factor(ifelse(deltaBeta <0, "Negative_skew", "Positive_skew"))) 

dfh_2 <- dfh %>% 
        filter(condition %in% c("PH", "PM"))        

change_b_vals_post <- b_vals[,rownames(dfh_2)] %>% 
        as.data.frame() %>% 
        mutate(average_homo = (`1PH`+ `2PH`+ `4PH`+ `5PH`+ `6PH`+ `7PH`+ `8PH`+`12PH`)/8,
               average_myo = (`1PM`+ `2PM`+ `4PM`+ `5PM`+ `6PM`+ `7PM`+ `8PM`+`12PM`)/8,
               change = average_myo-average_homo) %>% 
        dplyr::select(change) %>% 
        rownames_to_column(var = "cpg")


n_count2 <- merge(post_dmps, change_b_vals_post, by.x = "cpg") %>% 
        filter(dir == "Negative_skew") %>%    nrow()
        mutate(Relation_to_Island = factor(Relation_to_Island, levels = c("S_Shelf","N_Shelf","Island", "S_Shore", "N_Shore", "OpenSea"))) %>% 
        summarise(count(Relation_to_Island)) %>% 
        as.data.frame()

p_count2 <- merge(post_dmps, change_b_vals_post, by.x = "cpg") %>% 
        filter(dir == "Positive_skew") %>% 
        mutate(Relation_to_Island = factor(Relation_to_Island, levels = c("S_Shelf","N_Shelf","Island", "S_Shore", "N_Shore", "OpenSea"))) %>% 
        summarise(count(Relation_to_Island)) %>% 
        as.data.frame() 

df2 <- data.frame(
        category = n_count2[6:1,]$x,
        neg_skew = (n_count2[6:1,]$freq)*-1,
        pos_skew = p_count2[6:1,]$freq)

p2 <- ggplot(df2, aes(x = category)) +
        geom_bar(aes(y = neg_skew, fill = "neg_skew"), stat = "identity") +
        geom_bar(aes(y = pos_skew, fill = "pos_skew"), stat = "identity") +
        scale_fill_manual(values = c("#453781FF", "#DCE319FF")) + # set colors
        geom_label(aes(label = as.integer(neg_skew)*-1, x = category, y = as.integer(neg_skew)), size = 4, alpha = 0.8, nudge_y = -3000)+
        geom_label(aes(label = as.integer(pos_skew), x = category, y = as.integer(pos_skew)), size = 4, alpha = 0.8, nudge_y = 3000)+
        coord_flip()+
        theme_classic()+
        theme(legend.title = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_text(hjust = 0.58))+
        labs(y = "Number of DMPs",
             title = "Post DMPs in Homogenate vs. Myonuclei")


# find number of significant cpgs that were altered more than <0.01 B < 0.1 B < 0.2 B <


change_b_vals_post %>% 
        mutate(change = abs(change)) %>% # change all numbers to positive
        merge(post_dmps,.,by.x = "cpg") %>% 
        # find values in ranges
        mutate(skew = ifelse(change <0.01, "<0.01", 
                             ifelse(change >0.01 & change < 0.1, "0.01<0.1", 
                                    ifelse(change > 0.1 & change < 0.2, "0.1<0.2", 
                                           ifelse(change > 0.2 & change < 0.3, "0.2<0.3",">0.3"))))) %>% 
        summarise(count(skew)) -> df3

df4 <- data.frame(category = factor(df3[,1]$x, levels = c("<0.01", "0.01<0.1", "0.1<0.2", "0.2<0.3", ">0.3")),
                  count = df3[,1]$freq)

p4 <- ggplot(data = df4, aes(x = category))+
        geom_bar(aes(y = count), stat = "identity")+
        geom_text(aes(label = count, y = count), vjust = -0.5)+
        theme_classic()+
        labs(title = "number of DMPs by skewed B-values at post",
             x = "B-value skew")+
        theme(axis.title.y = element_blank())+
        scale_y_continuous(trans = "log2", expand  = c(0, 1), n.breaks = 9)


plot_grid(p1,p3,p2,p4, rel_widths = c(1.5,1), labels = c("A", "B", "C", "D"))



# count baseline DMPs for islands within promoters

myDMP_BH_BM$BH_to_BM %>%
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        merge(., anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island") %>% 
        nrow()

# count post dmps for islnds within promoters

myDMP_PH_PM$PH_to_PM %>%
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        merge(., anno, by = "cpg") %>% 
        filter(Relation_to_Island == "Island") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island") %>% 
        nrow()

unique(anno$Regulatory_Feature_Group)

################################################################################################################################

# create heatmap of dmps

###################################################################################

# add annotation of condition

dfh <- data.frame(sample = as.character(colnames(baseline_dmps_in_islands)), condition = "Timepoint") %>% 
        column_to_rownames("sample") %>% 
        mutate(condition = substr_right(rownames(.),2))


# filter bvals for only conditions of interest

dfh_1 <- dfh %>% 
        filter(condition %in% c("BH", "BM"))

dmp_list <- rownames(myDMP_BH_BM$BH_to_BM)

b_vals <- flt_beta$beta

# create list of cpgs that map to islands

island_cpg <- anno %>% 
        filter(Relation_to_Island == "Island",
               cpg %in% rownames(b_vals))


# filter b vals for significant island cpgs

b_vals_1 <- b_vals[island_cpg$cpg,rownames(dfh_1)] %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% dmp_list) 


# repeat for post samples

dmp_post <- rownames(myDMP_PH_PM$PH_to_PM)

b_vals_2 <- b_vals[island_cpg$cpg,rownames(dfh_2)] %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% dmp_post)

dfh$condition <- factor(dfh$condition, levels = c("BH","BM","PH","PM"))




baseline_dmps_in_islands <- merge(b_vals_1, b_vals_2, by.x = "cpg") %>% 
        na.omit() %>% 
        dplyr::select("1BH","2BH","4BH","5BH","6BH","7BH","8BH","12BH",
                      "1BM","2BM","4BM","5BM","6BM","7BM","8BM","12BM",
                      "1PH","2PH","4PH","5PH","6PH","7PH","8PH","12PH",
                      "1PM","2PM","4PM","5PM","6PM","7PM","8PM","12PM")


# plot baseline

baseline_dmps <- pheatmap(t(baseline_dmps_in_islands[,rownames(dfh_1)]), annotation_row = dfh_1, cluster_rows = FALSE, gaps_row = 8, 
         show_colnames = FALSE, annotation_names_row = FALSE, 
         color = viridis(n = 100), scale = "none")

# plot post

post_dmps <- pheatmap(t(baseline_dmps_in_islands[,rownames(dfh_2)]), annotation_row = dfh_2, cluster_rows = FALSE, gaps_row = 8, 
                          show_colnames = FALSE, annotation_names_row = FALSE, 
                          color = viridis(n = 100), scale = "none")

# plot together



plot_heat <- plot_grid(baseline_dmps[[4]], post_dmps[[4]], ncol=1, labels = c("A", "B"))   

ggsave2(plot_heat, filename = "heatmap_homo_vs_myo.pdf",units = "cm", width = 19, height = 21, bg = "white")




# plot promoter assosiated islands in a heatmap


heat_1 <- b_vals_1 %>% 
        dplyr::select(`cpg`,`1BH`,`2BH`,`4BH`,`5BH`,`6BH`,`7BH`,`8BH`,`12BH`,
                      `1BM`,`2BM`,`4BM`,`5BM`,`6BM`,`7BM`,`8BM`,`12BM`) %>% 
        merge(., anno, by = "cpg") %>% 
        filter(Relation_to_Island == "Island" & Regulatory_Feature_Group == "Promoter_Associated") %>% 
        dplyr::select(`1BH`,`2BH`,`4BH`,`5BH`,`6BH`,`7BH`,`8BH`,`12BH`,
                      `1BM`,`2BM`,`4BM`,`5BM`,`6BM`,`7BM`,`8BM`,`12BM`) %>% 
        t() %>% 
        pheatmap(., annotation_row = dfh_1, gaps_row = 8, 
                 show_colnames = FALSE, annotation_names_row = FALSE, 
                 color = viridis(n = 100), scale = "none", cluster_rows = FALSE)


heat_2 <- b_vals_2 %>% 
        dplyr::select(`cpg`, `1PH`,`2PH`,`4PH`,`5PH`,`6PH`,`7PH`,`8PH`,`12PH`,
                      `1PM`,`2PM`,`4PM`,`5PM`,`6PM`,`7PM`,`8PM`,`12PM`) %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Relation_to_Island == "Island" & Regulatory_Feature_Group == "Promoter_Associated") %>% 
        dplyr::select(`1PH`,`2PH`,`4PH`,`5PH`,`6PH`,`7PH`,`8PH`,`12PH`,
                      `1PM`,`2PM`,`4PM`,`5PM`,`6PM`,`7PM`,`8PM`,`12PM`) %>% 
        t() %>% 
        pheatmap(., annotation_row = dfh_2, gaps_row = 8, 
                 show_colnames = FALSE, annotation_names_row = FALSE, 
                 color = viridis(n = 100), scale = "none", cluster_rows = FALSE)

promoter_heat <- plot_grid(heat_1[[4]], heat_2[[4]], ncol=1, labels = c("A","B"))






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

# DMP venn diagrams

#################################################################################

# create 4 cpg lists of DMPs at baseline homo, baseline myo, post homo and post myo


baseline_dmps <- myDMP_BH_BM$BH_to_BM %>% 
        as.data.frame() %>%
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, deltaBeta) 


# merge with newer annotation

baseline_dmps <- merge(baseline_dmps, anno, by = "cpg") %>% 
        as.data.frame() %>% 
        mutate("dir" = as.factor(ifelse(deltaBeta <0, "Negative_skew", "Positive_skew")))

post_dmps <- myDMP_PH_PM$PH_to_PM %>% 
        as.data.frame() %>%
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, deltaBeta)

# merge with newer annotation

post_dmps <- merge(post_dmps, anno, by = "cpg") %>% 
        as.data.frame() %>% 
        mutate("dir" = as.factor(ifelse(deltaBeta <0, "Negative_skew", "Positive_skew")))




base_neg <- merge(baseline_dmps, change_b_vals_baseline, by.x = "cpg") %>% 
        filter(dir == "Negative_skew")

base_pos <- merge(baseline_dmps, change_b_vals_baseline, by.x = "cpg") %>% 
        filter(dir == "Positive_skew")

post_neg <- merge(post_dmps, change_b_vals_post, by.x = "cpg") %>% 
        filter(dir == "Negative_skew")

post_pos <- merge(post_dmps, change_b_vals_post, by.x = "cpg") %>% 
        filter(dir == "Positive_skew")


# add the four different cpg lists to a shared list


DMPs <- list(
        " " = base_neg$cpg,
        Baseline = base_pos$cpg,
        Post = post_neg$cpg,
        "  " = post_pos$cpg
)

devtools::install_github("yanlinlin82/ggvenn")

library(ggvenn)

# plot venn-diagram with only baseline and post sides

venn <- ggvenn(DMPs, set_name_size = 10, stroke_size = 1, 
       fill_color = c("#453781FF", "#DCE319FF","#453781FF", "#DCE319FF"),text_size = 10,stroke_alpha = 0.8)





# find the two cpgs that were positively skewed at baseline and negatively at post

merge(base_pos, post_neg, by = "cpg") 

# find the significant position in homo DMPs across time

myDMP_BH_PH$BH_to_PH %>% 
        rownames_to_column(var = "cpg") %>% 
        rownames_to_column(var = "number") %>% 
        filter(cpg == "cg00408117" | cpg == "cg13897145" | cpg == "cg06394887")



Illumina_anno %>% 
        filter(IlmnID == "cg00408117" | IlmnID == "cg13897145"| IlmnID == "cg06394887")


# find the one probe that is neg baseline and pos post

merge(base_neg,post_pos, by = "cpg")

# find the same pobes in myonuclear DMPs list

myDMP_BM_PM$BM_to_PM %>% 
        rownames_to_column(var = "cpg") %>% 
        rownames() %>% 
        length()
        rownames_to_column(var = "number") %>% 
        filter(cpg == "cg00408117" | cpg == "cg13897145"| cpg == "cg06394887")

        

##############################################################################
        
# volcano plot of DMPs after RT
        
#############################################################################
        
# plot DMPs based on M-value delta change and p value
        
myDMP_BH_PH_bmiq$BH_to_PH %>% 
                as.data.frame() %>% 
                dplyr::select(P.Value, deltaBeta) %>% 
                mutate(color_comb = log2(abs(deltaBeta))) %>% 
                ggplot(aes(x = log2(deltaBeta), y = log10(P.Value)))+
                geom_point(aes(color = (color_comb)))+
                scale_color_viridis()+
                scale_y_reverse()

###############################################################################

# self organizing map (som)

################################################################################
library(kohonen)

# create dataframe with group averages

som_data <-  B2M(b_vals)%>% 
        as.data.frame() %>% 
        mutate(avg_BH = round((`1BH`+`2BH`+`4BH`+`5BH`+`6BH`+`7BH`+`8BH`+`12BH`)/8,2),
               avg_BM = round((`1BM`+`2BM`+`4BM`+`5BM`+`6BM`+`7BM`+`8BM`+`12BM`)/8,2),
               avg_PH = round((`1PH`+`2PH`+`4PH`+`5PH`+`6PH`+`7PH`+`8PH`+`12PH`)/8,2),
               avg_PM = round((`1PM`+`2PM`+`4PM`+`5PM`+`6PM`+`7PM`+`8PM`+`12PM`)/8,2)) %>% 
        dplyr::select(avg_BH, avg_BM, avg_PH, avg_PM) %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::filter(cpg %in% c(rownames(myDMP_BH_BM$BH_to_BM),rownames(myDMP_PH_PM$PH_to_PM)))


rownames(som_data) <- som_data[,1]
som_data <- som_data[,-1]

# transform into matrix   
        
som_data1 <- as.matrix(t(som_data))


# create grid of dimentions which will give us the total amount of groups to cluster into. 2 x 2 = 4 groups

som_grid <- somgrid(xdim = 2, ydim = 2, topo = "rectangular")

# set seed for reproducebility
set.seed(1)

# run model
som_model <- som(som_data1, grid = som_grid)

som_data[,"cluster"] <- som_model$unit.classif


count(som_model$unit.classif)

som_data %>% 
        pivot_longer(names_to = "condition", values_to = "mean_M", cols = 1:4) %>% 
        mutate(condition = factor(condition, levels = c("avg_BH", "avg_BM", "avg_PH", "avg_PM")),
               cluster = factor(cluster)) %>% 
        dplyr::group_by(cluster,condition) %>% 
        dplyr::summarize(m = mean(mean_M),
                  s = sd(mean_M)) %>% 
        ggplot(aes(x = condition, y = m))+
        geom_line(aes(group = cluster))+
        geom_errorbar(aes(ymin = m-s,
                          ymax = m+s), width = 0.2)+
        facet_grid(~cluster)

# do kmeans clustering instead of som

### run kmeans clustering on som_data

set.seed(1)


xyss <- vector()

for (i in 1) {
        xyss[i] <- sum(kmeans(som_data,i)$withinss)
}
plot(2:10, xyss, type = "b", main = "clusters of M-value profiles", xlab = "number of clusters", ylab = "XYSS")

# elbow at ~4 clusters

# make cluster with identified elbow
set.seed(1)
kmeans <- kmeans(som_data, 4, iter.max = 300, nstart = 10)

# add cluster number to som data and replot

k_means <-as.data.frame(kmeans$cluster)

som_data[,"kmeans"] <- k_means[,1]

count(k_means[,1])
#kmeans_plot <- 
som_data %>% 
        pivot_longer(names_to = "condition", values_to = "mean_M", cols = 1:4) %>% 
        mutate(condition = factor(condition, levels = c("avg_BH", "avg_BM", "avg_PH", "avg_PM")),
               kmeans = factor(kmeans))%>% 
        dplyr::group_by(kmeans,condition) %>% 
        dplyr::summarize(m = mean(mean_M),
                         s = sd(mean_M)) %>% 
        mutate(sample = substr_right(as.character(condition), 1),
               timepoint = substr_right(as.character(condition), 2)) %>% 
        mutate(sample = ifelse(sample == "M", "Myonuclei", "Homogenate"),
               timepoint = substr(timepoint, 1,1)) %>% 
        mutate(timepoint = ifelse(timepoint == "B", "Baseline", "Post")) %>% 
        mutate(timepoint = factor(timepoint, levels = c("Baseline", "Post")),
               sample = factor(sample, levels = c("Myonuclei","Homogenate"))) %>% 
        ggplot(aes(x = timepoint, y = m, color = sample))+
        geom_line(aes(group = sample), size = 1.5, position = position_dodge(width = 0.2))+
        geom_errorbar(aes(ymin = m-s,
                          ymax = m+s, x = timepoint, group = sample),inherit.aes = FALSE, width = 0.2, size = 1, position = position_dodge(width = 0.2), alpha = 0.5)+
        facet_grid(~kmeans, labeller = as_labeller(c(`1`="1: N = 56239", `2`="2: N = 49432",`3`="3: N = 52972",`4`="4: N = 29792")))+
        labs(y = "M-value")+
        theme_bw(base_size = 20)+
        theme(axis.text.x = element_text(size = 15),
              axis.title.x = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank())+
        geom_text(aes(label = round(m,2), y = m),position = position_dodge(width = 1), vjust = 1.5,fontface = "bold", size = 4.5)+
        labs(tag = "D")+
        theme(plot.tag.position = c(0.01,0.97),
              plot.tag = element_text(face = "bold"))

                
        
# remake som plot with individual y axis scaling


k1 <- som_data %>% 
        pivot_longer(names_to = "condition", values_to = "mean_M", cols = 1:4) %>% 
        mutate(condition = factor(condition, levels = c("avg_BH", "avg_BM", "avg_PH", "avg_PM")),
               kmeans = factor(kmeans))%>% 
        dplyr::group_by(kmeans,condition) %>% 
        dplyr::summarize(m = mean(mean_M),
                         s = sd(mean_M)) %>% 
        mutate(sample = substr_right(as.character(condition), 1),
               timepoint = substr_right(as.character(condition), 2)) %>% 
        mutate(sample = ifelse(sample == "M", "Myonuclei", "Homogenate"),
               timepoint = substr(timepoint, 1,1)) %>% 
        mutate(timepoint = ifelse(timepoint == "B", "Baseline", "Post")) %>% 
        mutate(timepoint = factor(timepoint, levels = c("Baseline", "Post")),
               sample = factor(sample, levels = c("Myonuclei","Homogenate"))) %>% 
        filter(kmeans == "1") %>% 
        ggplot(aes(x = timepoint, y = m, color = sample))+
        geom_line(aes(group = sample), size = 1.5, position = position_dodge(width = 0.2))+
        geom_errorbar(aes(ymin = m-s,
                          ymax = m+s, x = timepoint, group = sample),inherit.aes = FALSE, width = 0.2, size = 1, position = position_dodge(width = 0.2), alpha = 0.5)+
        labs(y = "M-value",
             title = "Cluster 1: N = 56239")+
        theme_bw(base_size = 20)+
        theme(axis.text.x = element_text(size = 15),
              legend.position = "none",
              axis.title.x = element_blank())+
        geom_text(aes(label = round(m,2), y = m),position = position_dodge(width = 1), vjust = 1.5,fontface = "bold", size = 8)


k2 <- som_data %>% 
        pivot_longer(names_to = "condition", values_to = "mean_M", cols = 1:4) %>% 
        mutate(condition = factor(condition, levels = c("avg_BH", "avg_BM", "avg_PH", "avg_PM")),
               kmeans = factor(kmeans))%>% 
        dplyr::group_by(kmeans,condition) %>% 
        dplyr::summarize(m = mean(mean_M),
                         s = sd(mean_M/sqrt(length(mean_M)))) %>% 
        mutate(sample = substr_right(as.character(condition), 1),
               timepoint = substr_right(as.character(condition), 2)) %>% 
        mutate(sample = ifelse(sample == "M", "Myonuclei", "Homogenate"),
               timepoint = substr(timepoint, 1,1)) %>% 
        mutate(timepoint = ifelse(timepoint == "B", "Baseline", "Post")) %>% 
        mutate(timepoint = factor(timepoint, levels = c("Baseline", "Post")),
               sample = factor(sample, levels = c("Myonuclei","Homogenate"))) %>% 
        filter(kmeans == "2") %>% 
        ggplot(aes(x = timepoint, y = m, color = sample))+
        geom_line(aes(group = sample), size = 1.5, position = position_dodge(width = 0.2))+
        geom_errorbar(aes(ymin = m-s,
                          ymax = m+s, x = timepoint, group = sample),inherit.aes = FALSE, width = 0.2, size = 1, position = position_dodge(width = 0.2), alpha = 0.5)+
        labs(y = "M-value",
             title = "Cluster 2: N = 49432")+
        theme_bw(base_size = 20)+
        theme(axis.text.x = element_text(size = 15),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.minor = element_blank())+
        scale_y_continuous(limits = c(-2.25,-0.90))
        geom_text(aes(label = round(m,2), y = m),position = position_dodge(width = 1), vjust = 1.5,fontface = "bold", size = 8)


k3 <- som_data %>% 
        pivot_longer(names_to = "condition", values_to = "mean_M", cols = 1:4) %>% 
        mutate(condition = factor(condition, levels = c("avg_BH", "avg_BM", "avg_PH", "avg_PM")),
               kmeans = factor(kmeans))%>% 
        dplyr::group_by(kmeans,condition) %>% 
        dplyr::summarize(m = mean(mean_M),
                         s = sd(mean_M/sqrt(length(mean_M)))) %>% 
        mutate(sample = substr_right(as.character(condition), 1),
               timepoint = substr_right(as.character(condition), 2)) %>% 
        mutate(sample = ifelse(sample == "M", "Myonuclei", "Homogenate"),
               timepoint = substr(timepoint, 1,1)) %>% 
        mutate(timepoint = ifelse(timepoint == "B", "Baseline", "Post")) %>% 
        mutate(timepoint = factor(timepoint, levels = c("Baseline", "Post")),
               sample = factor(sample, levels = c("Myonuclei","Homogenate"))) %>% 
        filter(kmeans == "3") %>% 
        ggplot(aes(x = timepoint, y = m, color = sample))+
        geom_line(aes(group = sample), size = 1.5, position = position_dodge(width = 0.2))+
        geom_errorbar(aes(ymin = m-s,
                          ymax = m+s, x = timepoint, group = sample),inherit.aes = FALSE, width = 0.2, size = 1, position = position_dodge(width = 0.2), alpha = 0.5)+
        labs(y = "M-value",
             title = "Cluster 3: N = 52972")+
        theme_bw(base_size = 20)+
        theme(axis.text.x = element_text(size = 15),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.minor = element_blank())+
        geom_text(aes(label = round(m,2), y = m),position = position_dodge(width = 1), vjust = 1.5,fontface = "bold", size = 8)


k4 <- som_data %>% 
        pivot_longer(names_to = "condition", values_to = "mean_M", cols = 1:4) %>% 
        mutate(condition = factor(condition, levels = c("avg_BH", "avg_BM", "avg_PH", "avg_PM")),
               kmeans = factor(kmeans))%>% 
        dplyr::group_by(kmeans,condition) %>% 
        dplyr::summarize(m = mean(mean_M),
                         s = sd(mean_M/sqrt(length(mean_M)))) %>% 
        mutate(sample = substr_right(as.character(condition), 1),
               timepoint = substr_right(as.character(condition), 2)) %>% 
        mutate(sample = ifelse(sample == "M", "Myonuclei", "Homogenate"),
               timepoint = substr(timepoint, 1,1)) %>% 
        mutate(timepoint = ifelse(timepoint == "B", "Baseline", "Post")) %>% 
        mutate(timepoint = factor(timepoint, levels = c("Baseline", "Post")),
               sample = factor(sample, levels = c("Myonuclei","Homogenate"))) %>% 
        filter(kmeans == "4") %>% 
        ggplot(aes(x = timepoint, y = m, color = sample))+
        geom_line(aes(group = sample), size = 1.5, position = position_dodge(width = 0.2))+
        geom_errorbar(aes(ymin = m-s,
                          ymax = m+s, x = timepoint, group = sample),inherit.aes = FALSE, width = 0.2, size = 1, position = position_dodge(width = 0.2), alpha = 0.5)+
        labs(y = "M-value",
             title = "Cluster 4: N = 29792")+
        theme_bw(base_size = 20)+
        theme(axis.text.x = element_text(size = 15),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.minor = element_blank())+
        scale_y_continuous(limits = c(-4.15,-3))+
        geom_text(aes(label = round(m,2), y = m),position = position_dodge(width = 1), vjust = 1.5,fontface = "bold", size = 8)

grobs <- ggplotGrob(k4)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


library(ggpubr)

ggarrange(k1,k2,k3,k4,common.legend = TRUE, legend = "right", nrow = 1,widths = c(1,0.9,0.9,0.9))


plot_grid(k1,k2,k3,k4, nrow = 1)

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

# DMR on BMIQ norm data

myDMR_BH_PH_bmiq <- champ.DMR(beta = B2M(bmiq_norm),
                         pheno = flt_beta$pd$Sample_Group,
                         compare.group = c("BH", "PH"),
                         arraytype = "EPIC",
                         method = "Bumphunter",
                         adjPvalDmr = "none",
                         cores = 4)

myDMR_BM_PM_bmiq <- champ.DMR(beta = B2M(bmiq_norm),
                              pheno = flt_beta$pd$Sample_Group,
                              compare.group = c("BM", "PM"),
                              arraytype = "EPIC",
                              method = "Bumphunter",
                              adjPvalDmr = "none",
                              cores = 4)


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

saveRDS(unique_cpg, "unique_cpg.RDATA")


unique_cpg <- readRDS("unique_cpg.RDATA")


# remove "NA" from unique_cpg list

unique_cpg2 <- unique_cpg[-2]




# remove cpg that do not exist in filtered b vals

cpgs <- rownames(flt_beta$beta)

for (i in 1:length(unique_cpg)) {
        unique_cpg2[[i]] <- intersect(unique_cpg[[i]], cpgs)
        print(i)
}

subsets <- unique_cpg2
saveRDS(subsets, "subsets.RDATA")





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

gt.options(trim = TRUE) # quietly removes/skips cpgs in subset file that are not present in the filtered m_vals


gt_homo <- gt(response = condition_homo$condition, 
              alternative = t(m_vals_homo), 
              subsets = subsets, 
              model = "linear", 
              direction = FALSE, 
              permutations = 0,
              standardize = TRUE) # equal weighting of cpgs in gene 

saveRDS(gt_homo, "gt_homo.RDATA")


gt_myo <- gt(response = condition_myo$condition, 
              alternative = t(m_vals_myo), 
              subsets = subsets, 
              model = "linear", 
              direction = FALSE, 
              permutations = 0,
              standardize = TRUE) # equal weighting of cpgs in gene 

saveRDS(gt_myo, "gt_myo.RDATA")


gt_homo <- readRDS("gt_homo.RDATA")

gt_myo <- readRDS("gt_myo.RDATA")


gt_homo_res <- gt_homo@result %>% 
        as.data.frame() %>% 
        na.omit() %>% 
        dplyr::select(1:4, Cov = '#Cov') %>% 
        filter(Cov > 1 ) %>% 
        filter(`p-value` < 0.05) 

gt_myo_res <- gt_myo@result %>% 
        as.data.frame() %>% 
        na.omit() %>% 
        dplyr::select(1:4, Cov = '#Cov') %>% 
        filter(Cov > 1 ) %>% 
        filter(`p-value` < 0.05)



# add direction by finding the ???





# make ggplot of the global test

# homogenate

for (i in 1:length(rownames(gt_homo_res))) {
        
        x <- rownames(gt_homo_res[i,])
        
        z <- subsets[[x]]
        
        y <- b_vals_homo[z,] %>%        # change this to myo b_vals for myonuclei
                as.data.frame() %>% 
                rownames_to_column(var = "cpg")
        
        y <- y %>% 
                pivot_longer(names_to = "FP", values_to = "b_vals", cols = !1)
         
        y <- merge(y, anno, by = "cpg")
        
        y <- merge(y,rownames_to_column(condition_homo, var = "FP"), by = "FP")
        
        y <- y %>% 
                mutate(condition = factor(ifelse(condition == 0, "Baseline", "Post"), levels = c("Baseline", "Post")))
        
        p1 <- ggplot(data = y, aes(x = cpg, y = b_vals, group = condition, color = condition))+
                geom_point(aes(shape = Relation_to_Island), size = 1)+
                geom_smooth(method = "loess", se = FALSE, size = 1)+
                theme_classic(base_size = 8)+
                theme(axis.text.x = element_text(angle = -90),
                      axis.title.x = element_blank())+ 
                scale_y_continuous(limits = c(0,1))+
                labs(title = paste(x, "homogenate", "(", "p-value  = ", round(gt_homo_res[i,1], 4),")" ))
        
        ggsave(p1, filename = paste("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/Figures/global test homo genes/", "p=",round(gt_homo_res[i,1], 4),"_",x,"_", "homogenate",".png", sep = ""),)
        
        
        
        print(i)
}

# myonuclei

b_vals_myo <- flt_beta$beta[,myonuclei_samples$Sample_Name]


for (i in 1:length(rownames(gt_myo_res))) {
        
        x <- rownames(gt_myo_res[i,])
        
        z <- subsets[[x]]
        
        y <- b_vals_myo[z,] %>%        
                as.data.frame() %>% 
                rownames_to_column(var = "cpg")
        
        y <- y %>% 
                pivot_longer(names_to = "FP", values_to = "b_vals", cols = !1)
        
        y <- merge(y, anno, by = "cpg")
        
        y <- merge(y,rownames_to_column(condition_myo, var = "FP"), by = "FP")
        
        y <- y %>% 
                mutate(condition = factor(ifelse(condition == 0, "Baseline", "Post"), levels = c("Baseline", "Post")))
        
        p1 <- ggplot(data = y, aes(x = cpg, y = b_vals, group = condition, color = condition))+
                geom_point(aes(shape = Relation_to_Island), size = 1)+
                geom_smooth(method = "loess", se = FALSE, size = 1)+
                theme_classic(base_size = 8)+
                theme(axis.text.x = element_text(angle = -90),
                      axis.title.x = element_blank())+ 
                scale_y_continuous(limits = c(0,1))+
                labs(title = paste(x, "myonuclei", "(", "p-value  = ", round(gt_myo_res[i,1], 4),")" ))
        
        ggsave(p1, filename = paste("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/Figures/global test myo genes/", "p=",round(gt_myo_res[i,1], 4),"_",x,"_", "myonuclei",".png", sep = ""),)
        
        
        
        print(i)
}



# add direction on global delta-Beta value change

# mean beta post - mean beta baseline = negative (hypomethylation) positive (hyper)

gt_homo_res

b_vals_homo

unique_cpg2

x <- as.data.frame(rowMeans(t(b_vals_homo[unique_cpg2[[1]],]))) %>% 
        rownames_to_column()

y <- condition_homo %>% 
        rownames_to_column() %>% 
        as.data.frame() %>% 
        mutate(condition = ifelse(condition == 0, "Baseline", "Post"))

z <- merge(x,y, by = "rowname") %>% 
        group_by(condition) %>% 
        summarize(mean_b = mean(rowMeans(t(b_vals_homo[unique_cpg2[[1]], ])))) %>% 
        pivot_wider(names_from = condition, values_from = mean_b, data = .) %>% 
        mutate(change = Post - Baseline,
               percent_change = (change/Baseline)*100) %>% 
        mutate(gene = names(unique_cpg2[1]))

change_df <- z

# calculate change for all genes

for (i in 1:length(unique_cpg2)) {
        tryCatch({
        
        x <- as.data.frame(rowMeans(t(b_vals_homo[unique_cpg2[[i]],]))) %>% 
                rownames_to_column()
        
        colnames(x) <- c("rowname", "means")
        
        z <- merge(x,y, by = "rowname") %>% 
                group_by(condition) %>% 
                summarize(mean_b = mean(means, na.rm = TRUE)) %>% 
                pivot_wider(names_from = condition, values_from = mean_b, data = .) %>% 
                mutate(change = Post - Baseline,
                       percent_change = (change/Baseline)*100) %>% 
                mutate(gene = names(unique_cpg2[i]))
        
        change_df <- change_df %>% 
                rows_insert(z)
        
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        print(i)
}


# filter for change

change_df[change_df$change > 0.01 | change_df$change < -0.01,]

change_df_homo <- change_df

# merge and save unfiltered global test, and mean change

gt_homo_results2 <- gt_homo_results %>% 
        rownames_to_column(var = "gene") %>%
        merge(., change_df_homo, by = "gene", all.x = TRUE)

saveRDS(gt_homo_results2, "gt_homo_results2.RDATA")


# repeat for myonuclei

x <- as.data.frame(rowMeans(t(b_vals_myo[unique_cpg2[[1]],]))) %>% 
        rownames_to_column()

y <- condition_myo %>% 
        rownames_to_column() %>% 
        as.data.frame() %>% 
        mutate(condition = ifelse(condition == 0, "Baseline", "Post"))

z <- merge(x,y, by = "rowname") %>% 
        group_by(condition) %>% 
        summarize(mean_b = mean(rowMeans(t(b_vals_myo[unique_cpg2[[1]], ])))) %>% 
        pivot_wider(names_from = condition, values_from = mean_b, data = .) %>% 
        mutate(change = Post - Baseline,
               percent_change = (change/Baseline)*100) %>% 
        mutate(gene = names(unique_cpg2[1]))

change_df <- z

# calculate change for all genes

for (i in 1:length(unique_cpg2)) {
        tryCatch({
                
                x <- as.data.frame(rowMeans(t(b_vals_myo[unique_cpg2[[i]],]))) %>% 
                        rownames_to_column()
                
                colnames(x) <- c("rowname", "means")
                
                z <- merge(x,y, by = "rowname") %>% 
                        group_by(condition) %>% 
                        summarize(mean_b = mean(means, na.rm = TRUE)) %>% 
                        pivot_wider(names_from = condition, values_from = mean_b, data = .) %>% 
                        mutate(change = Post - Baseline,
                               percent_change = (change/Baseline)*100) %>% 
                        mutate(gene = names(unique_cpg2[i]))
                
                change_df <- change_df %>% 
                        rows_insert(z)
                
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        print(i)
}

change_df[change_df$change > 0.01 | change_df$change < -0.01,]

change_df_myo <- change_df

gt_myo_results <- gt_myo@result %>%
        as.data.frame() %>% 
        rownames_to_column(var = "gene") %>%
        merge(., change_df_homo, by = "gene", all.x = TRUE)




################################################################################################################################

# identify direction of DMP affect of interstitial cells

################################################################################

b_vals <- as.data.frame(flt_beta$beta)

head_b <- head(b_vals) %>% 
        as.data.frame()

b_vals <- b_vals %>% 
        mutate("1BI" = b_vals$"1BM" - b_vals$"1BH",
               "1PI" = b_vals$"1PM" - b_vals$"1PH",
               "2BI" = b_vals$"2BM" - b_vals$"2BH",
               "2PI" = b_vals$"2PM" - b_vals$"2PH",
               "4BI" = b_vals$"4BM" - b_vals$"4BH",
               "4PI" = b_vals$"4PM" - b_vals$"4PH",
               "5BI" = b_vals$"5BM" - b_vals$"5BH",
               "5PI" = b_vals$"5PM" - b_vals$"5PH",
               "6BI" = b_vals$"6BM" - b_vals$"6BH",
               "6PI" = b_vals$"6PM" - b_vals$"6PH",
               "7BI" = b_vals$"7BM" - b_vals$"7BH",
               "7PI" = b_vals$"7PM" - b_vals$"7PH",
               "8BI" = b_vals$"8BM" - b_vals$"8BH",
               "8PI" = b_vals$"8PM" - b_vals$"8PH",
               "12BI" = b_vals$"12BM" - b_vals$"12BH",
               "12PI" = b_vals$"12PM" - b_vals$"12PH")          # subtract homogenate profile from myonuclear profile to identify in which direction interstitial cells are dragging the myonuclear profile

summary(b_vals)

b_vals_int <- b_vals %>% 
        mutate("1BI" = b_vals$"1BM" - b_vals$"1BI",
               "1PI" = b_vals$"1PM" - b_vals$"1PI",
               "2BI" = b_vals$"2BM" - b_vals$"2BI",
               "2PI" = b_vals$"2PM" - b_vals$"2PI",
               "4BI" = b_vals$"4BM" - b_vals$"4BI",
               "4PI" = b_vals$"4PM" - b_vals$"4PI",
               "5BI" = b_vals$"5BM" - b_vals$"5BI",
               "5PI" = b_vals$"5PM" - b_vals$"5PI",
               "6BI" = b_vals$"6BM" - b_vals$"6BI",
               "6PI" = b_vals$"6PM" - b_vals$"6PI",
               "7BI" = b_vals$"7BM" - b_vals$"7BI",
               "7PI" = b_vals$"7PM" - b_vals$"7PI",
               "8BI" = b_vals$"8BM" - b_vals$"8BI",
               "8PI" = b_vals$"8PM" - b_vals$"8PI",
               "12BI" = b_vals$"12BM" - b_vals$"12BI",
               "12PI" = b_vals$"12PM" - b_vals$"12PI") %>% 
        dplyr::select(33:48)

summary(b_vals_int)

b <- as.matrix(b_vals_int)
densityPlot(b)

unique(anno$Regulatory_Feature_Group)
# use champ.refbase to identify and correct for immune cell profiles from blood


myRefbase_int <- champ.refbase(beta = b_vals_int,                      ### returns b-vals adjusted for 5 main cellpopulations identified
                           arraytype = "EPIC")

myRefbase_myo <- champ.refbase(beta = b_vals_myo,                      ### returns b-vals adjusted for 5 main cellpopulations identified
                               arraytype = "EPIC")

myRefbase_homo <- champ.refbase(beta = b_vals_homo,                      ### returns b-vals adjusted for 5 main cellpopulations identified
                               arraytype = "EPIC")

# not satisfied with how this works

################################################################################

### identify GSEA

##############################################################################

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





