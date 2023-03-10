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


#install package to install bioconductor packages


install.packages("BiocManager")


# install required packages

install.packages("tidyverse")


BiocManager::install("AnnotationDbi")

BiocManager::install("globaltest")



library(globaltest);library(AnnotationDbi); library(tidyverse)



#  set working direcotry to where you have the RData file stored





# file descriptions


# f_bVals       <-      full dataset
# mygenes       <-      unfiltered subset of dataset
# mygenes2      <-      filtered subset of dataset (use this)
# a             <-      list of column conditions in f_bVals
# gtResults     <-      example output from gt() function



# skip next part


####################################################################################################

# create subset list

mygenes <- list()

for (i in 1:length(uniquecg$`unique(CGtoGENE2$ENTREZ)`)) {
        index <- which(array2$ENTREZ == uniquecg[i,])
        z <- array2[index,]
        mygenes[[(as.character(uniquecg[i,]))]] <- z$rowname
        print(i)
        
}


# filter elements with less than 3 element

mygenes2 = mygenes
shortgenes = list()

for (i in rev(1:(length(mygenes2)))) {
        x <- length(mygenes2[[i]])
        if(x<3){                                                                ### filter out genes with less than 3 cpgs
                shortgenes[[(as.character(uniquecg[i,]))]] <- mygenes2[[i]]     ### list of filtered 
                mygenes2[i] <- NULL
        }
        print(i)
}




#########################################################################################################

### Det er fra her jeg trenger hjelpen din Oliver <3


# hensikten jeg prøver å oppnå er å objektivt estimere enring i genuttrykket basert på reguleringen av det genuttrykket. 
#
# Vi har ikke et faktisk mål på genuttrykket, men vi har målt en mengde punkter på DNAet som sammen regulerer hvert enkelt gen. problemet er at 
#
# vi har forskjellig mengde punkter som styrer uttrykket av hvert gen, og tradisjonelle tester er derfor biased mot gener som har mange målepunkter. 
#
# Jeg prøver derfor å komme med en samlet score, uavhengig av antallet målepunkter for å estimere denne reguleringen av genuttrykk. Kanskje det finnes en 
#
# bedre statistisk test jeg kan gjøre? men dette er det beste som har blitt brukt i litteraturen sålangt! 

### all info om pakken

# https://www.bioconductor.org/packages/devel/bioc/vignettes/globaltest/inst/doc/GlobalTest.pdf

# artikkel som bruker den og forklarer matten

# https://academic.oup.com/bioinformatics/article/20/1/93/229017






### jeg ønsker å bruke funksjonen:

gt()

# og muligens, men disse er mindre viktig! Begynne med bare gt()

gtKEGG()

gtGO()





# mitt forsøk på funksjonen

globalt <- gt(a ~., subsets = mygenes2, data = t(f_bVals)) # doesn't work/runs for hours with no end


# derfor prøvde jeg på subset med bare 3, men det funka fortsatt ikke! 

subset <- head(mygenes2, n = 3)

globalt <- gt(a ~., subsets = subset, data = t(f_bVals)) 

