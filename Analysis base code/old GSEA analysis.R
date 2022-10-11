################################################################
###                                                          ###
### Gene Set Enrichment analysis                             ###
###                                                          ###
################################################################


### if running for the first time you need to install several packages first, othervise skip straight to the loading og the libraries

install.packages("BiocManager")

BiocManager::install("minfi")

BiocManager::install("FDb.InfiniumMethylation.hg19")

BiocManager::install("wateRmelon")

BiocManager::install("methylumi")

BiocManager::install("maxprobes")

install.packages("tidyverse")

install.packages("ggplot2")

BiocManager::install("pathview")

BiocManager::install("gage")

BiocManager::install("gageData")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")


### now load the libraries


library(wateRmelon); library(methylumi);library(FDb.InfiniumMethylation.hg19);library(minfi); library(maxprobes); library(tidyverse);library(ggplot2)
library(tidyverse); library(pathview); library(gage); library(gageData); library(ChAMP); library(org.Hs.eg.db); library(AnnotationDbi)


### continue from "Oshlack workflow filtering and QC"


mVals <- getM(MsetExProbes)

myDMP <- champ.DMP(beta = mVals)                ### adjusted p value correction method = "Benjamini-Hochberg"
                                                ### works on mVals and bVals

df2 = as.data.frame(myDMP$C_to_T)               ### C_to_T = conditions between which we have identified DMPs

df2 = df2[order(abs(df2$logFC), decreasing = TRUE),]



#######################################################
#######################################################
#######################################################
###
### visualizations before doing GSEA, can skip if not interested! 
###
###

### volcano plot of log FC just to visualize DMPs

ggplot(df2, aes(logFC, -log10(adj.P.Val)))+
        geom_point()

#########################################################
###
### tornado plot of DMP in shores, islands etc.

df2 <- df2 %>% 
        mutate(direction = ifelse(logFC < 0, "hypo (-)", "hyper (+)"))

numb <- df2 %>%
        group_by(cgi) %>%
        summarise(count=n()) 

numb = numb[order(abs(numb$count), decreasing = TRUE),]

positions <- c("shelf","island","shore","opensea")     #### change this order based on decreasing order of "numb"

### hyper and hypo counts for annotation

numb2 <- df2 %>%
        group_by(cgi,direction) %>%
        summarise(count="n()")

a = numb2[3,3]
b = numb2[4,3]
c = numb2[7,3]
d = numb2[8,3]
e = numb2[1,3]
f = numb2[2,3]
g = numb2[5,3]
h = numb2[6,3]



ggplot(df2, aes(x = logFC, y = cgi, fill = direction))+
        geom_bar(stat = "identity")+
        theme_classic()+
        scale_fill_manual(values = c("red3", "springgreen4"))+          
        scale_y_discrete(limits = positions)+
        annotate(geom = "label", x = 8000, y = 4, label = a)+
        annotate(geom = "label", x = -9700, y = 4, label = b)+
        annotate(geom = "label", x = 6300, y = 3, label = c)+
        annotate(geom = "label", x = -2950, y = 3, label = d)+
        annotate(geom = "label", x = 7800, y = 2, label = e)+
        annotate(geom = "label", x = -1500, y = 2, label = f)+
        annotate(geom = "label", x = 2200, y = 1, label = g)+
        annotate(geom = "label", x = -1850, y = 1, label = h)+           ### labels are number of DMPs in each region, while x is the sum of the logFC (needs to be changed)
        labs(fill = "Methylation")


#######################################################
#######################################################
#######################################################


############################################################
###
### map ENTREZ ID to gene name
###
############################################################

columns(org.Hs.eg.db)

df2[!(is.na(df2$gene) | df2$gene==""), ]-> df3

df3$ENTREZ = mapIds(org.Hs.eg.db, 
                    key = as.character(df3$gene), 
                    column = "ENTREZID", 
                    keytype ="ALIAS",
                    multiVals = "first")

############################################################
###
### map out logFC to GO pathway
###
############################################################




foldchange = df3$logFC

names(foldchange) = df3$ENTREZ

data("go.sets.hs")
data("go.subs.hs")


gobpsets = go.sets.hs[go.subs.hs$BP]            ### BP is for biological processes

gobpdf3 = gage(exprs = foldchange, gsets = gobpsets, same.dir = TRUE)


view(gobpdf3$less)                              ### view less expressed GO pathways C vs. T
view(gobpdf3$greater)                           ### view more expressed GO pathways C vs. T




############################################################
###
### map out logFC to KEGG pathway
###
############################################################


data("kegg.sets.hs")
data("sigmet.idx.hs")   ### smaller subset dataset that contains only signalling and metabolic pathways

# data("sig.idx.hs")    ### which elements in kg.sets are signaling pathways
# data("met.idx.hs")    ### which elements in kg.sets are metabolism pathways
# data("dise.idx.hs")    ### which elements in kg.sets are disease pathways


kegg.subset = kegg.sets.hs[sigmet.idx.hs]       ### this now contains only the subset kegg pathwyas

### only signalling and metabolic pathways

keggres = gage(exprs = foldchange, gsets = kegg.subset, same.dir = TRUE)

view(keggres$less)                              ### view less expressed KEGG pathways C vs. T
view(keggres$greater)                           ### view more expressed KEGG pathways C vs. T

### all kegg pathways

keggres2 = gage(exprs = foldchange, gsets = kegg.sets.hs, same.dir = TRUE)

view(keggres2$less)                              ### view less expressed KEGG pathways C vs. T
view(keggres2$greater)                           ### view more expressed KEGG pathways C vs. T


### Map KEGG pathways to pdf

keggrespathways = data.frame(id = rownames(keggres$less), keggres$less) %>% 
        tibble::as.tibble() %>% 
        filter(row_number() <= 20) %>% 
        .$id %>% 
        as.character()
as.data.frame(keggrespathways)

keggresids = substr(keggrespathways, start = 1, stop = 8)
keggresids

tmp = sapply(keggresids, function(pid) pathview(gene.data = foldchange, pathway.id = pid, species = "hsa"))




###############################################################



keggres = gage(exprs = df$logFC, gsets = kegg.subset, same.dir = TRUE)






