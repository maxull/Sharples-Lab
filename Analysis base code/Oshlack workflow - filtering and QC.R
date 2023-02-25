###################################################################
###                                                             ###
### Methylation analysis - Oshlack workflow  filtering and QC   ###
###                                                             ###
###################################################################

### https://dockflow.org/workflow/methylation-array-analysis/#content






library(ChAMP); library(RColorBrewer); library(limma); library(maxprobes)



### direct R to your data


### description of how to load your data: https://bioconductor.org/packages/devel/bioc/vignettes/minfi/inst/doc/minfi.html#3_Reading_data


testDir=system.file("extdata", package = "ChAMPdata")           ### here you need to change to where your data is stored on your computer

### load data

targets <- read.metharray.sheet(testDir)

rgSet <- read.metharray.exp(targets = targets)

### merge targets and rgSet

targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID



#################################
###
### filtering
###
#################################


### filter out samples with mean P-values > 0.05

detP <- detectionP(rgSet)


keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]


targets <- targets[keep,]


detP <- detP[,keep]



###############################
###
### functional normalization
###
##################################


mSetSq <- preprocessFunnorm(rgSet)


detP <- detP[match(featureNames(mSetSq), rownames(detP)),]


### filter p_values <0.01


keep <- rowSums(detP < 0.01) == ncol(mSetSq)


mSetSqFlt <- mSetSq[keep,]


### if male and female samples, filter out XY chromosomes


ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in%  c("chrX","chrY")])



### filter out SNPs

mSetSqFlt2 <- dropLociWithSnps(mSetSqFlt)


### filter out cross reactive probes from Chen et al. 2013, Benton et al. 2015 for 450k data
### and Pidsley et al. 2016 and McCartney et al. 2016 for 850k


MsetExProbes <- dropXreactiveLoci(mSetSqFlt2)
