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




d <- as.data.frame(d)



d <- t(g1)



d2 <- left_join(d, pd_gt, by = "colnames")

d2 %>% 
        aov(cg00000029  ~Sample_Group, data = .) %>% 
        summary()




cpgs <-rownames(g1)



d3 <- d2[,-1]

lm(.~as.factor(Sample_Group), data = d3)
anova(mod)


class(d3$Sample_Group)

d3$Sample_Group <- as.factor(d3$Sample_Group)



