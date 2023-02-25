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

library(ENmix); library(FedData)

### set working directory to filepath that contains the data of interest
### this folder should contain the raw red and green signal .idat files and a CSV file of the targets

### save directory for later! 

mypath <- getwd()

#############################################################
### unzipping of zip files when downloading from GEO database (change file paths)


gunzip("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/GC-AS-10179.tar.gz", remove=FALSE)

untar("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/GC-AS-10179.tar", exdir = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/raw idat/")



#############################################################
#############################################################

### load data from idat files and filter

# load raw data to compare CPGs

myData <- champ.import(directory = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/raw idat/", 
                       arraytype = "EPIC")

myLoad <- champ.load(directory = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/raw idat/", 
                     arraytype = "EPIC", method = "minfi")

length(rownames(myData$beta))
length(rownames(myLoad$beta))
length(rownames(getBeta(myNorm)))

# notes from running champ load
####################################################################
# Failed CpG Fraction.
# 4BH          0.0009978820
# 2BH          0.0007463926
# 5PM          0.0006425668
# 5BM          0.0007337028
# 12BM         0.0006367987
# 1PM          0.0007452390
# 5BH          0.0009851921
# 8PH          0.0020326798
# 2BM          0.0006517957
# 12BH         0.0006010364
# 12PM         0.0006887116
# 7PM          0.0004995178
# 1PH          0.0005272047
# 12PH         0.0006310306
# 7PH          0.0005848857
# 6PH          0.0008536794
# 2PM          0.0006887116
# 8BM          0.0008409895
# 7BH          0.0005468162
# 4BM          0.0005341264
# 4PM          0.0005041323
# 6BH          0.0004672164
# 8PM          0.0013566580
# 1BM          0.0012909016
# 5PH          0.0008871344
# 7BM          0.0006241088
# 1BH          0.0005987292
# 6BM          0.0005191293
# 4PH          0.0005618133
# 6PM          0.0007152449
# 8BH          0.0012885944
# 2PH          0.0012897480

# Filtering probes with a detection p-value above 0.01 in one or more samples has removed 5517 probes from the analysis. If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples.
# << Filter DetP Done. >>
# 
# 
# There is no NA values in your matrix, there is no need to imputation.
# 
# Filtering probes with a beadcount <3 in at least 5% of samples, has removed 15774 from the analysis.
# << Filter Beads Done. >>
# 
# Filtering non-cg probes, has removed 2896 from the analysis.
# << Filter NoCG Done. >>
# 
# Using general EPIC SNP list for filtering.
# Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper, 2016, has removed 95783 from the analysis.
# << Filter SNP Done. >>
# 
# Filtering probes that align to multiple locations as identified in Nordlund et al, has removed 23 from the analysis.
# << Filter MultiHit Done. >>
# 
# Filtering probes on the X or Y chromosome has removed 16404 from the analysis.
# << Filter XY chromosome Done. >>
# 
# [Beta value is selected as output.]
# 
# Zeros in your dataset have been replaced with smallest positive value.
# 
# One in your dataset have been replaced with largest value below 1.
# 
# The analysis will be proceed with 730439 probes and 32 samples.
# 
# Current Data Set contains 0 NA in [Beta] Matrix.


####################################################################

CpG.GUI(arraytype = "EPIC")

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

png(file="C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/Figures/median beta.png",
    units="in",res=reso,height=length,width=length*2)


long_beta %>% 
        mutate(FP = ifelse(nchar(FP) < 4, paste(0, FP, sep = ""), FP)) %>% 
        mutate(NR = substr(FP,1,2)) %>% 
        ggplot(aes(x = FP, y = beta, color = NR))+
        geom_violin()+
        geom_boxplot(color = "black", width = 0.2)+
        theme(axis.text.x = element_text(angle = -90))+
        labs(title = "Raw beta values")+
        geom_point(data = beta_mean, aes(x = FP, y = m), inherit.aes = FALSE, size = 3)

dev.off()

long_beta %>% 
        group_by(FP) %>% 
        summarise(m = mean(beta),
                  s = sd(beta))%>% 
        ggplot(aes(x = FP, y = m))+
        geom_point()+
        geom_errorbar(aes(ymin = m-s, ymax = m+s), width = 0.1)





### functional normalization

myNorm <- champ.norm(beta = myLoad$beta,
                     rgSet = myLoad$rgSet,
                     method ="BMIQ",
                     arraytype = "EPIC")

# did not work, have to do it manually

myNorm_fun <- preprocessFunnorm(myLoad$rgSet, nPCs=2, sex = NULL, bgCorr = TRUE,
                  dyeCorr = TRUE, keepCN = FALSE, ratioConvert = TRUE,
                  verbose = TRUE)

myNorm_fun@colData@listData[["Sample_Group"]] <- substr_right(myLoad$pd$Sample_Name,2)
myData$pd$Sample_Group <- substr_right(myLoad$pd$Sample_Name,2)
myImport = myData


myImport$beta <- getBeta(myNorm_fun)

length(rownames(myImport$beta))

densityPlot(myImport$beta, sampGroups=myLoad$pd$Sample_Group,
            main="Normalized", legend=TRUE)
        
flt_beta <- champ.filter(beta = myImport$beta,
                         pd = myImport$pd,
                         detP = myImport$detP[rownames(myImport$beta),],
                         beadcount = myImport$beadcount[rownames(myImport$beta),],
                         Meth = myImport$Meth[rownames(myImport$beta),],
                         UnMeth = myImport$UnMeth[rownames(myImport$beta),],
                         arraytype = "EPIC")

# notes from running champ norm
####################################################################
# [preprocessFunnorm] Background and dye bias correction with noob
# [preprocessFunnorm] Mapping to genome
# [preprocessFunnorm] Quantile extraction
# [preprocessFunnorm] Normalization
# Warning message:
#         In .getSex(CN = CN, xIndex = xIndex, yIndex = yIndex, cutoff = cutoff) :
#         An inconsistency was encountered while determining sex. One possibility is that only one sex is present. We recommend further checks, for example with the plotSex function.

####################################################################

champ.QC(beta = flt_beta$beta)
### add group in myLoad$pd



myLoad$pd[,4] <- substr_right(myLoad$pd$Sample_Name,2)


### visualize nomralization and save density plots


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

library(limma)
library(RColorBrewer)

pal <- brewer.pal(8,"Dark2")
plotMDS(getM(myNorm), top=1000, gene.selection="common",
        col=pal[factor(myNorm$Sample_Group)])



### PCA analysis

# https://bioinfo4all.wordpress.com/2021/01/31/tutorial-6-how-to-do-principal-component-analysis-pca-in-r/

install.packages(c("factoextra", "FactoMineR"))

library(factoextra)
library(FactoMineR)

flt_beta$beta

pca.data <- PCA(myData$beta, scale.unit = FALSE, graph = FALSE)


fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 100))




champ.QC()      # default QC on myLoad

QC.GUI()



getwd()

setwd("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/")


### save temporary files

saveRDS(EPIC.manifest.hg19, file = "EPIC.manifest.hg19.RDATA")

saveRDS(multi.hit, file = "multi.hit.RDATA")

saveRDS(myLoad, file = "myLoad.RDATA")

saveRDS(probe.features, file = "probe.features.RDATA")

saveRDS(myNorm, file = "myNorm.RDATA")


### reload temporary files (not all are neaded for downstream analysis!)

EPIC.manifest.hg19 <- readRDS("EPIC.manifest.hg19.RDATA")

multi.hit <- readRDS("multi.hit.RDATA")

myLoad <- readRDS("myLoad.RDATA")

probe.features <- readRDS("probe.features.RDATA")

myNorm <- readRDS("myNorm.RDATA")

### set working directory back to original


setwd(mypath)


################################################################################################################################
################################################################################################################################
################################################################################################################################

### identify DMPs

myDMP <- champ.DMP(beta = getBeta(myNorm),
          arraytype = "EPIC")                   ### significant DMPs with BH adjusted p val of 0.05


# BH adjusted p val of 0.05 gave no significant DMPs between BH to PH, and BM to PM

# run these comparisons with unadjusted p val

myDMP2<- champ.DMP(beta = getBeta(myNorm),
                   arraytype = "EPIC",
          adjust.method = "none",
          adjPVal = 0.05)               ### worked, no p val adjustement


# run DMPs on M vals instead of B vals


mvals <- B2M(getBeta(myNorm))

saveRDS(mvals, file = "mvals.RDATA")


myDMP3 <- champ.DMP(beta = mvals, 
                    arraytype = "EPIC")


myDMP4 <- champ.DMP(beta = mvals, 
                    arraytype = "EPIC",
                    adjust.method = "none",
                    adjPVal = 0.05)

### visualize DMPs

DMP.GUI(myDMP[[2]])


saveRDS(myDMP, "myDMP.RDATA")

saveRDS(myDMP2, "myDMP2.RDATA")


myDMP <- readRDS("myDMP.RDATA")

myDMP2 <- readRDS("myDMP2.RDATA")

#### check number of dmps at different unadjusted pvals 0.05, 0.01, 0.001

myDMP4$BH_to_PH %>% 
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        filter(P.Value < 0.05) %>% 
        rownames() %>% 
        length()

# done to here

####

# to do
# - DMR
# - GSEA
# - global test etc.



################################################################################################################################

library(pheatmap)
library(viridis)



# create heatmap of dmps

# add annotation of condition

dfh <- data.frame(sample = as.character(colnames(b_vals)), condition = "Timepoint") %>% 
        column_to_rownames("sample") %>% 
        mutate(condition = substr_right(rownames(.),2))


# filter bvals for only conditions of interest

dfh_1 <- dfh %>% 
        filter(condition %in% c("BH", "BM"))

dmp_list <- rownames(myDMP$BH_to_BM)

b_vals <- getBeta(myNorm)

b_vals_1 <- b_vals[, rownames(dfh_1)]
        

dmp_bvals <- b_vals_1[dmp_list,]   # filter f_bVals for only significant DMPs

df_bvals <- head(dmp_bvals, n = 1000)



pheatmap(t(df_bvals), annotation_row = dfh_1, cutree_rows = 2, 
         show_colnames = FALSE, annotation_names_row = FALSE, 
         color = viridis(n = 100), scale = "none")
# all datapoints

pheatmap(t(dmp_bvals), annotation_row = dfh_1, cutree_rows = 2, 
         show_colnames = FALSE, annotation_names_row = FALSE, 
         color = viridis(n = 100), scale = "none")

# cant plot all DMPs
# isolate DMPs within islands

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

dmp_isl_change <- myDMP2$BH_to_PH %>%
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

DMP_BM_to_PM <- myDMP2$PM_to_BM %>%
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        mutate(logFC = logFC*-1) %>%
        arrange(logFC)

write.csv(DMP_BM_to_PM, "/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/DMPs/BM_to_PM_unadj. P 0.05.csv", row.names = FALSE)

DMP_BH_to_PH <- myDMP2$BH_to_PH %>%
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        arrange(logFC)

write.csv(DMP_BH_to_PH, "/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/DMPs/BH_to_PH_unadj. P 0.05.csv", row.names = FALSE)

DMP_BH_to_BM <- myDMP$BH_to_BM %>%
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        arrange(logFC)

write.csv(DMP_BH_to_BM, "/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/DMPs/BH_to_BM_FDR P 0.05.csv", row.names = FALSE)

DMP_PM_to_PH <- myDMP$PM_to_PH %>%
        rownames_to_column(var = "cpg") %>% 
        as.data.frame() %>% 
        mutate(logFC = logFC*-1) %>%
        arrange(logFC)

write.csv(DMP_PM_to_PH, "/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/DMPs/PH_to_PM_FDR P 0.05.csv", row.names = FALSE)




###

# cell pop correction

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










################################################################################################################################

### identify DMRs

myDMR <- champ.DMR(arraytype = "EPIC",
                   compare.group = c("BH", "BM"))          #no significant DMPs with BH adjusted p val of 0.05

myDMR <- champ.DMR(arraytype = "EPIC",
                   method = "Bumphunter",
                   adjPvalDmr = "none",
                   compare.group = c("BH", "BM"),
                   cores = 4)

### DMR did not work

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

doGT(pheno.v = pd1, 
     data.m = f_bVals,
     array = "850k")


################################################################################################################################
# manually do global test

doGT <- function(pheno.v,data.m,model=c("linear"),array=c("450k","850k"),ncores=4){
        if(array=="450k"){
                message(" Mapping 450k probes to genes... ")
                data("dualmap450kEID");
                subsets <- lapply(mapEIDto450k.lv,intersect,rownames(data.m),mc.cores = ncores);
                message(" Done ")
        }else {
                message(" Mapping EPIC probes to genes... ")
                data("dualmap850kEID");
                subsets <- mclapply(mapEIDto850k.lv,intersect,rownames(data.m),mc.cores = 1);
                message(" Done ")
        }
        nrep.v <- unlist(lapply(subsets,length));
        selG.idx <- which(nrep.v>0);
        message(" Running Global Test... ")
        gt.o <- gt(response=pheno.v,alternative=t(data.m),model=model,directional = FALSE, standardize = FALSE, permutations = 0, subsets=subsets[selG.idx],trace=F);
        message(" Done ")
        resGT.m <- as.matrix(result(gt.o));
        tmp.s <- sort(resGT.m[,2],decreasing=TRUE,index.return=TRUE);
        sresGT.m <- resGT.m[tmp.s$ix,];
        return(sresGT.m);
}

BiocManager::install("globaltest")

library(globaltest)


