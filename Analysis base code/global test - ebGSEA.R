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

BiocManager::install("ebGSEA")


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

# array of cpgs and gene names

array <- rownames_to_column(CGtoGENE2)

array2 <- left_join(array, df, by = "rowname")

array2 <- array2 %>% 
        na.omit()


########################################################################
### create mygenes2 list again


mygenes <- list()

for (i in 1:length(uniquecg$`unique(CGtoGENE2$ENTREZ)`)) {
        index <- which(array2$ENTREZ == uniquecg[i,])
        z <- array2[index,]
        mygenes[[(as.character(uniquecg[i,]))]] <- z$rowname
        print(i)
        
}


# for loop successful
# now it is based on array, wich still includes NAs
# array2 contains no NAs
# run loop on array 2, and remove all lists that contain less than 3 cpg sites

getwd()
setwd("C:/Users/maxul/Documents/Skole/Master 21-22/Master/")
save.image(file = "gt_workspace.RData")

mygenes2 = mygenes


# save removed lists in new list
shortgenes <- list()


for (i in rev(1:(length(mygenes2)))) {
        x <- length(mygenes2[[i]])
        if(x<3){                                                                ### filter out genes with less than 3 cpgs
                shortgenes[[(as.character(uniquecg[i,]))]] <- mygenes2[[i]]     ### list of filtered 
                mygenes2[i] <- NULL
        }
        print(i)
}



# example of global test output, run via ebGSEA function. however it is missing direction of the change/statistic

gtResults <- myGSEA3[["gtResult"]]


#convert f_bVals to expressionSet



library(Biobase)

bVals<-new("ExpressionSet", exprs=as.matrix(f_bVals))


####################################################################################






globalt <- gt(a ~., subsets = mygenes2, data = bVals) # doesn't work

    






gt_KEGG <- gtKEGG(a, bVals, probe2entrez = "eg", annotation = "org.Hs.eg.db")


gtKEGG(a, bVals, annotation = "org.Hs.eg.db", p.adjust = "BH")




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



##########################################################################
#######################################################################
#############################################################################

### chatGPT løsning


# prepare one gene dataset for chatGPT solution

gene1 <- as.data.frame(mygenes2[["5934"]]) %>% 
        select(cpgs = 1)


df <- rownames_to_column(f_bVals, var = "cpgs")
head(df)

gene_df <- merge(gene1, df, by.x = "cpgs", by.y = "cpgs", all.x = TRUE)

rownames(gene_df) <- gene_df$cpgs
        
gene_df1<- gene_df %>% 
        select(!1) %>% 
        t() %>% 
        as.data.frame() %>% 
        rownames_to_column()

gene_df2 <- merge(gene_df1, pd_gt, by.x = "rowname", by.y = "colnames")

data <- gene_df2 %>% 
        select(group = Sample_Group, 2:29)


data2 <- data %>% 
        select(2:29)

rownames(data2) <- data$group



# Load libraries
library(tidyverse)
library(limma)

# Load the data
data <- read.csv("methylation_data.csv")

# Create the model matrix
design <- model.matrix(~group, data)

# Fit the model
fit <- lmFit(t(data2), design)

# Perform the global test
globalTest <- eBayes(fit)

# Extract the p-values
pvals <- topTable(globalTest, coef=2, n=nrow(data))$P.Value

# Print the p-values
print(pvals)



eBayes(fit)


topTable(globalTest, n = Inf)

dim(globalTest)
names(globalTest)
volcanoplot(globalTest)
plotMD(globalTest)


# run gt on gene 1

gt(response = (as.numeric(as.factor(data$group))-1), alternative = data2, model = "logistic", direction = FALSE, permutations = 0)

# works

data3 = data2

data3$group = (as.numeric(as.factor(data$group))-1)

x = colnames(data2)

y = data3 %>% 
        mutate(average = rowMeans(data3)) %>% 
        group_by(group) %>% 
        summarize(mean = mean(average))

z <- y$mean        

# 0 = 7wk_loading, 
# 1 = baseline

# 0-1
z[1]-z[2]



# run gt on all subsets

gt_results <- gt(response = (as.numeric(as.factor(data$group))-1), alternative = t(f_bVals), model = "logistic", direction = FALSE, permutations = 0, subsets = mygenes2)

gt_results_df <- as.data.frame(gt_results@result)

# worked wohoooo


summary(gt_results)

z.score(gt_results)

covariates(gt_results[[2]])

subjects(gt_results[[2]])
sort(gt_results)

