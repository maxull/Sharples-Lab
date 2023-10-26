#
#
# Drugged dudes data


library(readxl)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gage)
library(pathview)
library(stringr)


RT_AS_vs_RT <- read_excel("/Users/maxul/Downloads/DGE_results_muscle_CoolMPS.xlsx", sheet = "Muscle_Group_Comparisons_4_of_6")


RT_AS_vs_RT %>% 
        filter(PValue < 0.05) %>% 
        dplyr::select("gene_names" = 1, 2:6) %>% 
        mutate(EntrezID = mapIds(org.Hs.eg.db, keys = gene_names, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")) %>% 
        na.omit() %>% 
        dplyr::select(EntrezID, logFC)  %>% 
        as.data.frame()-> exprsMat

KEGG_new <- kegg.gsets(species = "hsa", id.type = "kegg", check.new=TRUE)

rownames(exprsMat) <- exprsMat$EntrezID

exprsMat %>% 
        dplyr::select(logFC) %>% 
        as.matrix() -> exprsMat





kegg_res <- gage(exprs = exprsMat, gsets = KEGG_new$kg.sets, same.dir = TRUE, ref = NULL, samp = NULL)

kegg_up <- kegg_res$greater
kegg_down <- kegg_res$less



keggresids <- kegg_down %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "pathway") %>% 
        filter(p.val < 0.05) %>%
        mutate(id = substr(pathway, start = 1, stop = 8)) %>% 
        pull(id)



setwd("C:/Users/maxul/Documents/Skole/Master 21-22/")
tmp = sapply(keggresids, function(pid) pathview(gene.data = exprsMat, pathway.id = pid, species = "hsa", low = "blue", mid = "grey", high = "yellow"))


