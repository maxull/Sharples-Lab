###################################################################################################################
##################       MACS.  REACTOME GSEA.                                  ###################################
###################################################################################################################     


library(cowplot)
library(scales)
library(tidyverse)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)

BiocManager::install("signatureSearch")
library(signatureSearch)

library(ReactomePA)

BiocManager::install("DOSE")
library(DOSE)

BiocManager::install("enrichplot")
library(enrichplot)


########################################################################
#########               Load data                       ################
########################################################################


setwd("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/Skole/M.Sc/Master 21-22/Master/DATA/Epigenetics")

anno <- readRDS(file = "anno.RDATA")
DMPs_PH_vs_BH <- readRDS("DMPs_PH_vs_BH.RDATA")
DMPs_PM_vs_BM <- readRDS("DMPs_PM_vs_BM.RDATA")



########################################################################
# Convert gene symbols to Entrez Gene IDs

entrezIDs <- mapIds(org.Hs.eg.db, keys = anno$UCSC_RefGene_Name, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

anno$entrezIDs <- entrezIDs


########################################################################
# over representaion analysis with reactome
########################################################################

MYOINT_hypo_eids <- merge(DMPs_PH_vs_BH %>% 
                                  filter(p.value < 0.05,
                                         delta_M < 0),anno, by = "cpg") %>% 
        arrange(delta_M) %>% 
        na.omit(entrezIDs) %>% 
        pull(entrezIDs)

React_MYOINT_hypo <- enrichPathway(gene = MYOINT_hypo_eids, organism = "human")

MYOINT_hyper_eids <- merge(DMPs_PH_vs_BH %>% 
                                   filter(p.value < 0.05,
                                          delta_M > 0),anno, by = "cpg") %>% 
        arrange(delta_M) %>% 
        na.omit(entrezIDs) %>% 
        pull(entrezIDs)

React_MYOINT_hyper <- enrichPathway(gene = MYOINT_hyper_eids, organism = "human")

MYO_hypo_eids <- merge(DMPs_PM_vs_BM %>% 
                               filter(p.value < 0.05,
                                      delta_M < 0),anno, by = "cpg") %>% 
        arrange(delta_M) %>% 
        na.omit(entrezIDs) %>% 
        pull(entrezIDs)

React_MYO_hypo <- enrichPathway(gene = MYO_hypo_eids, organism = "human")

MYO_hyper_eids <- merge(DMPs_PM_vs_BM %>% 
                                filter(p.value < 0.05,
                                       delta_M > 0),anno, by = "cpg") %>% 
        arrange(delta_M) %>% 
        na.omit(entrezIDs) %>% 
        pull(entrezIDs)

React_MYO_hyper <- enrichPathway(gene = MYO_hyper_eids, organism = "human")

# get the results from reactome

React_MYO_hyper@result
React_MYO_hypo@result
React_MYOINT_hyper@result
React_MYOINT_hypo@result


########################################################################
# PLOTS: MYO+INT hyper
########################################################################

# change entrezIDs to gene symbol
React_MYOINT_hyper_readable <- setReadable(React_MYOINT_hyper, 'org.Hs.eg.db', 'ENTREZID')

# barplot

barplot(React_MYO_hypo, showCategory=20)+ scale_fill_gradient(low = "#F5C242", high = "#353F4F")+theme_classic()+scale_x_continuous(expand = c(0,0))

# dotplot
dotplot(React_MYO_hypo, showCategory=20)+ scale_fill_gradient(low = "#F5C242", high = "#353F4F")+theme_classic()

# interaction plot between Reactome terms
React_MYOINT_hyper_readable_pairwise <- pairwise_termsim(React_MYOINT_hyper_readable)

set.seed(123)
emapplot(React_MYOINT_hyper_readable_pairwise, cex.params = list(category_label = 0.6))+ scale_fill_gradient(low = "#F5C242", high = "#353F4F")+labs(fill = "p.adjust")

p <- treeplot(React_MYOINT_hyper_readable_pairwise, 
         showCategory = 30,
         hclust_method = "ward.D", 
         group_color = c("#5698AF", "#2F6B34", "#741B0C","#2B324F","#C3B29E"),
         nWords = 4,
         offset_tiplab = 1.5,
         offset = 80)+
        scale_color_gradient(low = "#F5C242", high = "#353F4F")+labs(color = "p.adjust")+
        theme(legend.position = "bottom")

p$layers[[4]]$aes_params$size <- 1.5
p$layers[[3]]$aes_params$size <- 4

p$layers[[3]]$mapping$x <- 70
p$layers[[4]]$mapping$x <- 68
p$layers[[4]]$mapping$xend <- 68
p

########################################################################
# PLOTS: MYO+INT hypo
########################################################################

# change entrezIDs to gene symbol
React_MYOINT_hypo_readable <- setReadable(React_MYOINT_hypo, 'org.Hs.eg.db', 'ENTREZID')

# barplot

barplot(React_MYOINT_hypo_readable, showCategory=20)+ scale_fill_gradient(low = "#F5C242", high = "#353F4F")+theme_classic()+scale_x_continuous(expand = c(0,0))

# dotplot
dotplot(React_MYOINT_hypo_readable, showCategory=20)+ scale_fill_gradient(low = "#F5C242", high = "#353F4F")+theme_classic()

# interaction plot between Reactome terms
React_MYOINT_hypo_readable_pairwise <- pairwise_termsim(React_MYOINT_hypo_readable)

set.seed(125)
emapplot(React_MYOINT_hypo_readable_pairwise, cex.params = list(category_label = 0.6))+ scale_fill_gradient(low = "#F5C242", high = "#353F4F")+labs(fill = "p.adjust")

p <- treeplot(React_MYOINT_hypo_readable_pairwise, 
              showCategory = 30,
              hclust_method = "ward.D", 
              group_color = c("#5698AF", "#2F6B34", "#741B0C","#2B324F","#C3B29E"),
              nWords = 4,
              offset_tiplab = 1.5,
              offset = 80)+
        scale_color_gradient(low = "#F5C242", high = "#353F4F")+labs(color = "p.adjust")+
        theme(legend.position = "bottom")

p$layers[[4]]$aes_params$size <- 1.5
p$layers[[3]]$aes_params$size <- 4

p$layers[[3]]$mapping$x <- 70
p$layers[[4]]$mapping$x <- 68
p$layers[[4]]$mapping$xend <- 68
p
########################################################################
# PLOTS: MYO hyper
########################################################################

# change entrezIDs to gene symbol
React_MYO_hyper_readable <- setReadable(React_MYO_hyper, 'org.Hs.eg.db', 'ENTREZID')

# barplot

barplot(React_MYO_hyper_readable, showCategory=20)+ scale_fill_gradient(low = "#F5C242", high = "#353F4F")+theme_classic()+scale_x_continuous(expand = c(0,0))

# dotplot
dotplot(React_MYO_hyper_readable, showCategory=20)+ scale_fill_gradient(low = "#F5C242", high = "#353F4F")+theme_classic()

# interaction plot between Reactome terms
React_MYO_hyper_readable_pairwise <- pairwise_termsim(React_MYO_hyper_readable)

set.seed(123)
emapplot(React_MYO_hyper_readable_pairwise, cex.params = list(category_label = 0.6))+ scale_fill_gradient(low = "#F5C242", high = "#353F4F")+labs(fill = "p.adjust")

p <- treeplot(React_MYO_hyper_readable_pairwise, 
              showCategory = 30,
              hclust_method = "ward.D", 
              group_color = c("#5698AF", "#2F6B34", "#741B0C","#2B324F","#C3B29E"),
              nWords = 4,
              offset_tiplab = 1.5,
              offset = 80)+
        scale_color_gradient(low = "#F5C242", high = "#353F4F")+labs(color = "p.adjust")+
        theme(legend.position = "bottom")

p$layers[[4]]$aes_params$size <- 1.5
p$layers[[3]]$aes_params$size <- 4

p$layers[[3]]$mapping$x <- 70
p$layers[[4]]$mapping$x <- 68
p$layers[[4]]$mapping$xend <- 68
p
########################################################################
# PLOTS: MYO hypo
########################################################################

# change entrezIDs to gene symbol
React_MYO_hypo_readable <- setReadable(React_MYO_hypo, 'org.Hs.eg.db', 'ENTREZID')

# barplot

barplot(React_MYO_hypo_readable, showCategory=20)+ scale_fill_gradient(low = "#F5C242", high = "#353F4F")+theme_classic()+scale_x_continuous(expand = c(0,0))

# dotplot
dotplot(React_MYO_hypo_readable, showCategory=20)+ scale_fill_gradient(low = "#F5C242", high = "#353F4F")+theme_classic()

# interaction plot between Reactome terms
React_MYO_hypo_readable_pairwise <- pairwise_termsim(React_MYO_hypo_readable)

set.seed(124)
emapplot(React_MYO_hypo_readable_pairwise, cex.params = list(category_label = 0.6))+ scale_fill_gradient(low = "#F5C242", high = "#353F4F")+labs(fill = "p.adjust")

p <- treeplot(React_MYO_hypo_readable_pairwise, 
              showCategory = 30,
              hclust_method = "ward.D", 
              group_color = c("#5698AF", "#2F6B34", "#741B0C","#2B324F","#C3B29E"),
              nWords = 4,
              offset_tiplab = 1.5,
              offset = 80)+
        scale_color_gradient(low = "#F5C242", high = "#353F4F")+labs(color = "p.adjust")+
        theme(legend.position = "bottom")

p$layers[[4]]$aes_params$size <- 1.5
p$layers[[3]]$aes_params$size <- 4

p$layers[[3]]$mapping$x <- 70
p$layers[[4]]$mapping$x <- 68
p$layers[[4]]$mapping$xend <- 68
p





########################################################################
# PLOTS: Pathway plot MYO
########################################################################

# calculate diff meth per gene

MYO_deltaM <- merge(DMPs_PM_vs_BM %>% 
                            filter(p.value < 0.05),anno, by = "cpg") %>% 
        arrange(delta_M) %>% 
        na.omit(entrezIDs) %>%
        distinct(entrezIDs, .keep_all = TRUE) %>% 
        pull(delta_M)

names(MYO_deltaM) <- merge(DMPs_PM_vs_BM %>% 
                                   filter(p.value < 0.05),anno, by = "cpg") %>% 
        arrange(-abs(delta_M)) %>% 
        na.omit(entrezIDs) %>%
        # head(20) %>% #pull(UCSC_RefGene_Name) -> filt_genes_myo
        distinct(entrezIDs, .keep_all = TRUE) %>% 
        pull(entrezIDs)


# visualize pathway without fold change
viewPathway(
        "Muscle contraction",
        organism = "human",
        readable = TRUE,
        foldChange = MYO_deltaM,
        keyType = "ENTREZID",
        layout = "kk"
)

# visualize pathway without fold change
viewPathway(
        "Protein-protein interactions at synapses",
        organism = "human",
        readable = TRUE,
        foldChange = MYO_deltaM,
        keyType = "ENTREZID",
        layout = "kk"
)

# visualize pathway without fold change
set.seed(123)
viewPathway(
        "Integration of energy metabolism",
        organism = "human",
        readable = TRUE,
        foldChange = MYO_deltaM,
        keyType = "ENTREZID",
        layout = "kk"
)

# visualize pathway without fold change
set.seed(123)
viewPathway(
        "Adipogenesis",
        organism = "human",
        readable = TRUE,
        foldChange = MYOINT_deltaM,
        keyType = "ENTREZID",
        layout = "kk"
)

merge(DMPs_PM_vs_BM %>% filter(p.value < 0.05), anno, by = "cpg") %>% filter(UCSC_RefGene_Name == "AR")

M_change[c("cg10164191","cg13873881","cg22908141"),]

########################################################################
# PLOTS: Pathway plot MYO+INT
########################################################################

# calculate diff meth per gene

MYOINT_deltaM <- merge(DMPs_PH_vs_BH %>% 
                            filter(p.value < 0.05),anno, by = "cpg") %>% 
        arrange(delta_M) %>% 
        na.omit(entrezIDs) %>%
        distinct(entrezIDs, .keep_all = TRUE) %>% 
        pull(delta_M)

names(MYOINT_deltaM) <- merge(DMPs_PH_vs_BH %>% 
                                   filter(p.value < 0.05),anno, by = "cpg") %>% 
        arrange(-abs(delta_M)) %>% 
        na.omit(entrezIDs) %>%
        # head(20) %>% #pull(UCSC_RefGene_Name) -> filt_genes_myo
        distinct(entrezIDs, .keep_all = TRUE) %>% 
        pull(entrezIDs)


# visualize pathway without fold change
viewPathway(
        "Muscle contraction",
        organism = "human",
        readable = TRUE,
        foldChange = MYOINT_deltaM,
        keyType = "ENTREZID",
        layout = "kk"
)

########################################################################
#       Get reactome DB
########################################################################

getDb <- function(organism) {
        if (organism == "worm") {
                organism = "celegans"
                warning("'worm' is deprecated, please use 'celegans' instead...")
        }
        
        annoDb <- switch(organism,
                         anopheles   = "org.Ag.eg.db",
                         arabidopsis = "org.At.tair.db",
                         bovine      = "org.Bt.eg.db",
                         canine      = "org.Cf.eg.db",
                         celegans    = "org.Ce.eg.db",
                         chicken     = "org.Gg.eg.db",
                         chimp       = "org.Pt.eg.db",
                         coelicolor  = "org.Sco.eg.db", 
                         ecolik12    = "org.EcK12.eg.db",
                         ecsakai     = "org.EcSakai.eg.db",
                         fly         = "org.Dm.eg.db",
                         gondii      = "org.Tgondii.eg.db",
                         human       = "org.Hs.eg.db",
                         malaria     = "org.Pf.plasmo.db",
                         mouse       = "org.Mm.eg.db",
                         pig         = "org.Ss.eg.db",
                         rat         = "org.Rn.eg.db",
                         rhesus      = "org.Mmu.eg.db",
                         xenopus     = "org.Xl.eg.db",
                         yeast       = "org.Sc.sgd.db",
                         zebrafish   = "org.Dr.eg.db",
        )
        return(annoDb)
}

get_Reactome_Env <- function() {
        if (!exists(".ReactomePA_Env", envir = .GlobalEnv)) {
                assign(".ReactomePA_Env", new.env(), .GlobalEnv)
        }
        get(".ReactomePA_Env", envir= .GlobalEnv)
}

getALLEG <- function(organism) {
        annoDb <- getDb(organism)
        require(annoDb, character.only = TRUE)
        annoDb <- eval(parse(text=annoDb))
        eg=keys(annoDb, keytype="ENTREZID")
        return(eg)
}

get_Reactome_DATA <- function(organism = "human") {
        ReactomePA_Env <- get_Reactome_Env()
        
        if (exists("organism", envir=ReactomePA_Env, inherits = FALSE)) {
                org <- get("organism", envir=ReactomePA_Env)
                if (org == organism &&
                    exists("PATHID2EXTID", envir = ReactomePA_Env) &&
                    exists("EXTID2PATHID", envir = ReactomePA_Env) &&
                    exists("PATHID2NAME",  envir = ReactomePA_Env)) {
                        
                        ## use_cached
                        return(ReactomePA_Env)
                }
        }
        
        ALLEG <- getALLEG(organism)
        
        EXTID2PATHID <- as.list(reactomeEXTID2PATHID)
        EXTID2PATHID <- EXTID2PATHID[names(EXTID2PATHID) %in% ALLEG]
        
        PATHID2EXTID <- as.list(reactomePATHID2EXTID) ## also contains reactions
        
        PATHID2NAME <- as.list(reactomePATHID2NAME)
        PI <- names(PATHID2NAME)
        ## > PATHID2NAME[['68877']]
        ## [1] "Homo sapiens: Mitotic Prometaphase" "Homo sapiens: Mitotic Prometaphase"
        PATHID2NAME <- lapply(PATHID2NAME, function(x) x[1])
        names(PATHID2NAME) <- PI
        
        PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% names(PATHID2NAME)]
        PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% unique(unlist(EXTID2PATHID))]
        PATHID2EXTID <- lapply(PATHID2EXTID, function(x) intersect(x, ALLEG))
        
        PATHID2NAME <- PATHID2NAME[names(PATHID2NAME) %in% names(PATHID2EXTID)]
        
        PATHID2NAME <- unlist(PATHID2NAME)
        PATHID2NAME <- gsub("^\\w+\\s\\w+:\\s+", "", PATHID2NAME) # remove leading spaces
        
        assign("PATHID2EXTID", PATHID2EXTID, envir=ReactomePA_Env)
        assign("EXTID2PATHID", EXTID2PATHID, envir=ReactomePA_Env)
        assign("PATHID2NAME", PATHID2NAME, envir=ReactomePA_Env)
        return(ReactomePA_Env)
}

# get reactome db
Reactome_DATA <- get_Reactome_DATA("human")


########################################################################
###             Heatplot of top 20 genes in pathways
########################################################################



########################################################################
###             betavalue pathway plot
########################################################################

# average the promoter methylation in my samples
average_promoter.m <- constAvBetaTSS(beta.m = beta, type = "850k")      # function from EpiSCORE package to calculate average methylation across promoter probes

average_promoter.m <- average_promoter.m %>% 
        as.data.frame() %>% 
        dplyr::select("1BH","2BH","4BH","5BH","6BH","7BH","8BH","12BH",
                      "1BM","2BM","4BM","5BM","6BM","7BM","8BM","12BM",
                      "1PH","2PH","4PH","5PH","6PH","7PH","8PH","12PH",
                      "1PM","2PM","4PM","5PM","6PM","7PM","8PM","12PM") %>% 
        mutate(BH = (`1BH`+`2BH`+`4BH`+`5BH`+`6BH`+`7BH`+`8BH`+`12BH`)/8,
               BM = (`1BM`+`2BM`+`4BM`+`5BM`+`6BM`+`7BM`+`8BM`+`12BM`)/8,
               PH = (`1PH`+`2PH`+`4PH`+`5PH`+`6PH`+`7PH`+`8PH`+`12PH`)/8,
               PM = (`1PM`+`2PM`+`4PM`+`5PM`+`6PM`+`7PM`+`8PM`+`12PM`)/8) %>% 
        # calculate change
        mutate(BM_vs_BH = BM-BH,
               PM_vs_PH = PM-PH,
               PH_vs_BH = PH-BH,
               PM_vs_BM = PM-BM)

# plot baseline homogenate

# visualize pathway without fold change

BH <- average_promoter.m %>% pull(BH)

names(BH) <- rownames(average_promoter.m)

set.seed(123)
viewPathway(
        "Adipogenesis",
        organism = "human",
        readable = TRUE,
        foldChange = BH,
        keyType = "ENTREZID",
        layout = "kk"
)

# plot post homogenate

PH <- average_promoter.m %>% pull(PH)

names(PH) <- rownames(average_promoter.m)

set.seed(123)
viewPathway(
        "Adipogenesis",
        organism = "human",
        readable = TRUE,
        foldChange = PH,
        keyType = "ENTREZID",
        layout = "kk"
)

# plot baseline homogenate

# visualize pathway without fold change

BM <- average_promoter.m %>% pull(BM)

names(BM) <- rownames(average_promoter.m)

set.seed(123)
viewPathway(
        "Adipogenesis",
        organism = "human",
        readable = TRUE,
        foldChange = BM,
        keyType = "ENTREZID",
        layout = "kk"
)

# plot post homogenate

PM <- average_promoter.m %>% pull(PM)

names(PM) <- rownames(average_promoter.m)

set.seed(123)
viewPathway(
        "Adipogenesis",
        organism = "human",
        readable = TRUE,
        foldChange = PM,
        keyType = "ENTREZID",
        layout = "kk"
)
