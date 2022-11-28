################################################################
###                                                          ###
### Global test -> ebGSEA                                    ###
###                                                          ###
################################################################

### Genome analysis and gene set enrichment analysis based on:

### Gene set enrichment analysis for genome-wide DNA methylation data
### https://doi.org/10.1101/2020.08.24.265702

### ChAMP: updated methylation analysis pipeline for Illumina BeadChips
### doi: 10.1093/bioinformatics/btx513


### extra functions https://github.com/YuanTian1991/ChAMP-Script


BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")


library(wateRmelon); library(methylumi);library(FDb.InfiniumMethylation.hg19);library(minfi); library(tidyverse);library(ggplot2)
library(tidyverse); library(pathview); library(gage); library(gageData); library(ChAMP); library(org.Hs.eg.db); library(AnnotationDbi);
library(reshape2)

mypath <- getwd()

setwd("/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/")


myLoad <- readRDS("myLoad.RDATA")

myNorm <- champ.norm(beta = myLoad$beta,
                     rgSet = myLoad$rgSet,
                     method ="BMIQ",
                     arraytype = "EPIC")                        ### bmiq worked, dont know why

myDMP<- champ.DMP(arraytype = "EPIC",
                  adjust.method = "none",
                  adjPVal = 0.05)               ### worked, no p val adjustement

df <- myDMP[["Baseline_to_7wk_Loading"]]


library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19); library(missMethyl)

annEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

CGtoGENE <- as.data.frame(annEPIC) %>% 
        select(Gene_Name = UCSC_RefGene_Name)

CGtoGENE$ENTREZ = as.character(mapIds(org.Hs.eg.db, 
                    key = as.character(CGtoGENE$Gene_Name), 
                    column = "ENTREZID", 
                    keytype ="ALIAS",
                    multiVals = first))

CGtoGENE <- CGtoGENE %>% 
        na.omit() %>% 
        filter(!ENTREZ == "NULL") %>% 
        select(2)

uniquecg <- as.data.frame(unique(CGtoGENE$ENTREZ))



###################################################################################
### convert gene ids to entrez ids

CGtoGENE2 <- as.data.frame(probe.features) %>% 
        select(gene)


CGtoGENE2$ENTREZ <- as.character(mapIds(org.Hs.eg.db,
                                        key = as.character(probe.features$gene),
                                        column = "ENTREZID",
                                        keytype ="ALIAS",
                                        multiVals = first))

CGtoGENE2 <- CGtoGENE2 %>% 
        na.omit() %>% 
        filter(!ENTREZ == "NULL") %>% 
        select(2)

uniquecg <- as.data.frame(unique(CGtoGENE2$ENTREZ)) %>% 
        na.omit()  ### 23329 unique entrez IDs instead of 16k from annEPIC file

#uniquegene <- as.data.frame(unique(probe.features$gene)) ### 26650 total unique gene ids in probe features
#uniquegene2 <- as.data.frame(unique(annEPIC$UCSC_RefGene_Name)) ### 66088 here

###

subsets <- lapply(CGtoGENE$ENTREZ, intersect, rownames(f_bVals))

nrep.v <- unlist(lapply(subsets,length))

print(nrep.v)

summarise(d)

################3
### custom large subset list generation for global test

a <- uniquecg[1,]
list(rownam)

f_bVals %>% 
        rownames_to_column() -> df

gene1 <- CGtoGENE2 %>% 
        filter(ENTREZ == a) %>% 
        rownames_to_column()

left_join(gene1, df, by = "rowname") -> g1


g1 %>% 
        na.omit() -> g1  
        g1[,-1:-2] -> g1
        
rownames(g1) <- g1[,1] 
g1[,-1]

g1 = as.data.frame(t(g1))

d <- rownames_to_column(d, var = "colnames")


################### trying global test
library(globaltest);library(AnnotationDbi)

gt(response = pd_gt, alternative = ~"Baseline" +"7wk_Loading", data = d)

###
#change pd1 to only contain rownames and group
pd1 %>% 
        select(colnames, Sample_Group) -> pd_gt

rownames(pd_gt) <- pd_gt[,1]

pd_gt <- pd_gt[,-1]







#################################################################################

# ivar input til for loop

array <- rownames_to_column(CGtoGENE2)

array2 <- left_join(array, df, by = "rowname")

array2 <- array2 %>% 
        na.omit()

#df %>% 
#        filter(rowname == "cg00001510")      # missing values in array 2... why are they missing? follow up



#for (i in 1: length(uniquecg$`unique(CGtoGENE2$ENTREZ)`)){
#index <- which(array2$ENTREZ == uniquecg[i,])
#z <- array2[index,]
#mygenes[[i]] <- z[,-1:-2]
#}

### funker men mangler å gi listene navnet til genene

mygenes <- list()



for (i in 1: length(uniquecg$`unique(CGtoGENE2$ENTREZ)`)){
        index <- which(array2$ENTREZ == uniquecg[i,])
        z <- array2[index,]
        mygenes[[(as.character(uniquecg[i,]))]] <- t(z[,-1:-2])
        print(i)
}

### works, however it returns some empty lists -> datafiltering beforehand should ensure that there are no empty lists! 

gene1 <- CGtoGENE2 %>% 
        filter(ENTREZ == "1947") %>% 
        rownames_to_column()

left_join(gene1, df, by = "rowname") -> g1


### i made new for loop that made the correct subsets, but its not saved... make again, should make output: mygenes2

colnames(f_bVals) <- a



########################################################################
### create mygenes2 list again


mygenes3 <- list()

for (i in 1:length(uniquecg$`unique(CGtoGENE2$ENTREZ)`)) {
        index <- which(array2$ENTREZ == uniquecg[i,])
        z <- array2[index,]
        mygenes3[[(as.character(uniquecg[i,]))]] <- z$rowname
        print(i)
        
}


# for loop successful
# now it is based on array, wich still includes NAs
# array2 contains no NAs
# run loop on array 2, and remove all lists that contain less than 3 cpg sites

getwd()
setwd("C:/Users/maxul/Documents/Skole/Master 21-22/Master/")
save.image(file = "gt_workspace.RData")

mygenes4 = mygenes3


for (i in 1:(length(mygenes3))) {
        x <- length(mygenes3[[i]])
        if(x<3){
                mygenes3[i] <- NULL
        }
        print(i)
}


# save removed lists in new list
shortgenes <- list()
shortgenes2 <- list()


for (i in 1:(length(mygenes3))) {
        x <- length(mygenes3[[i]])
        if(x<3){
                shortgenes[[(as.character(uniquecg[i,]))]] <- mygenes3[[i]]
                mygenes3[i] <- NULL
        }
        print(i)
}


for (i in 1:(length(mygenes4))) {
        x <- length(mygenes4[[i]])
        if(x<3){
                shortgenes[[i]] <- mygenes4[[i]]
                mygenes4[i] <- NULL
        }
        print(i)
}




shortgenes["100616112"]

mygenes3[["100616112"]]

####################################################################################






globalt <- gt(a ~., subsets = submygenes, data = t(f_bVals)) # doesn't work




# test one gene

one_gene <- as.data.frame(mygenes[[1]])
one_gene$timepoint <- a

###doesnt work


globalt <- gt(one_gene$timepoint~1, timepoint ~. , data = one_gene, model = "linear")   # doesnt work
### test gt
submygenes <- mygenes[3:12]


nrep.v <- unlist(submygenes, length)
selG.idx <- which(nrep.v>0)
a <- pd_gt$Sample_Group
b<- t(f_bVals)


gt.o <- gt(response = a, alternative = b, subsets = mygenes2)

#funker ikke

