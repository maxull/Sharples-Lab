############################################################################################
###                                                                                      ###
###  MACS: PCA analysis                                                                  ###
###                                                                                      ###
############################################################################################


#This script contains the majority of the analysis code for the master thesis of Max Ullrich

BiocManager::install("missMethyl")

library(cowplot)
library(ggrepel)
library(scales)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(missMethyl)
# library(ChAMP)
# library(pathview)
library(viridis)
library(ENmix)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gage)
library(EpiSCORE)


setwd("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/Skole/M.Sc/Master 21-22/Master/DATA/Epigenetics")

beta <- readRDS(file = "beta.RDATA")
anno <- readRDS(file = "anno.RDATA")

########################################################################

### DMP data

########################################################################

# calculate group averages

M_change <- beta %>%
        B2M() %>% 
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



#DMPs_BM_vs_BH <- readRDS("DMPs_BM_vs_BH.RDATA")
DMPs_PH_vs_BH <- readRDS("DMPs_PH_vs_BH.RDATA")
DMPs_PM_vs_BM <- readRDS("DMPs_PM_vs_BM.RDATA")
#DMPs_PM_vs_PH <- readRDS("DMPs_PM_vs_PH.RDATA")

# add adjusted p value for within time comparison
#library(stats)

#DMPs_BM_vs_BH$adj.p.val <- p.adjust(p = DMPs_BM_vs_BH$p.value, method = "fdr", n = length(rownames(DMPs_BM_vs_BH)))
#DMPs_PM_vs_PH$adj.p.val <- p.adjust(p = DMPs_PM_vs_PH$p.value, method = "fdr", n = length(rownames(DMPs_PM_vs_PH)))        

#p.adjust(p = DMPs_PM_vs_PH$p.value, method = "bonferroni", n = length(rownames(DMPs_PM_vs_PH))) %>% 
#        as.data.frame() %>% 
#        filter(.<0.05) 

# filter all dmp lists for p.value <=0.05


#DMPs_BM_vs_BH <- DMPs_BM_vs_BH %>% 
#        filter(adj.p.val < 0.05)

#DMPs_PM_vs_PH <- DMPs_PM_vs_PH %>% 
#        filter(adj.p.val < 0.05) 

#DMPs_PM_vs_BM <- DMPs_PM_vs_BM %>% 
#        filter(p.value < 0.05)

#DMPs_PH_vs_BH <- DMPs_PH_vs_BH %>% 
#        filter(p.value < 0.05)


########################################################################

### PCA

################################################################


# PCA 1 correlates with cell population

pca.out <- prcomp(t(M_change[,1:32]), scale. = FALSE)

plot(pca.out$x[,1:2])

# plot nice PCA plot

pca_plt <- pca.out$x[,1:2] %>% 
        as.data.frame() %>% 
        rownames_to_column(var = 'ID') %>% 
        mutate(Sample = factor(ifelse(grepl('BH', ID), 'Baseline_MYOINT', 
                               ifelse(grepl('BM', ID), 'Baseline_MYO', 
                                      ifelse(grepl('PH', ID), 'Post_MYOINT', 'Post_MYO'))),
                               levels = c('Baseline_MYO','Post_MYO','Baseline_MYOINT','Post_MYOINT'))) %>% 
        ggplot(aes(x = PC1, y = PC2))+
        #geom_point(aes(size = 2+pca.out$sdev),shape = 3)+
        geom_point(aes(fill = Sample), size = 6, stroke = 3, color = 'Black', shape = 21)+
        scale_fill_manual(values = c("#440154FF","#B396B9","#5DC863FF","#B3D7B1"), labels = c('MYO (Baseline)', 'MYO (Post)', 'MYO+INT (Baseline)','MYO+INT (Post)'))+
        theme_classic(base_size = 20)+
        labs(x = 'Principal Component 1',
             y = 'Principal Component 2')+
        theme(legend.position = 'top')+
        guides(fill = guide_legend(override.aes = list(size = 6)))+
        theme(
                panel.background = element_rect(fill='transparent'), #transparent panel bg
                plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                panel.grid.major = element_blank(), #remove major gridlines
                panel.grid.minor = element_blank(), #remove minor gridlines
                legend.background = element_rect(fill='transparent'), #transparent legend bg
                legend.box.background = element_rect(fill='transparent')) #transparent legend panel


ggsave(plot = pca_plt, bg = "transparent", path = "Figures/", filename = "MACS_PCA_plot_transparent.png") 

ggsave(filename = "Figures/MACS_PCA_plot_transparent.png", plot = pca_plt, bg = "transparent", device = "png")


loadings <- pca.out$rotation[,1:2]

# look through top PC1 drivers

sorted_loadings <- sort((loadings), decreasing = TRUE) 
top_loadings <- names(sorted_loadings)[1:100]



as.data.frame(top_loadings) %>% 
        dplyr::select("cpg" = top_loadings) %>% 
        merge(., anno, by = "cpg")


# identify drivers towards left, i.e. homogenate by fintering the negative values and arranging


loadings %>% 
        as.data.frame() %>% 
        dplyr::select("PC1" = 1) %>% 
        arrange(PC1) %>% 
        filter(PC1 < 0) %>%     # nrow = 374432
        head(10000) %>% 
        rownames_to_column(var = "cpg") %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(UCSC_RefGene_Name != "NA") -> neg_PC1

loadings %>% 
        as.data.frame() %>% 
        dplyr::select("PC1" = 1) %>% 
        arrange(-PC1) %>% 
        filter(PC1 > 0) %>%     # nrow = 374432
        head(10000) %>% 
        rownames_to_column(var = "cpg") %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(UCSC_RefGene_Name != "NA") -> pos_PC1


# plot correlation between top 100 drivers of PC1 between MYO and MYO+INT

### run ORA on top negative PCA1 probes

# get updated go pathways

go.sets <- go.gsets(species = "human", pkg.name=NULL, id.type = "EG", keep.evidence=FALSE)

# Run ORA


GO_neg_PCA1 <- gsameth(sig.cpg = neg_PC1$cpg,
                       all.cpg = rownames(beta), 
                       collection = go.sets$go.sets,
                       array.type = "EPIC")

GO_neg_PCA1 <- GO_neg_PCA1 %>% filter(FDR<0.05)

GO_pos_PCA1 <- gsameth(sig.cpg = pos_PC1$cpg,
                       all.cpg = rownames(beta), 
                       collection = go.sets$go.sets,
                       array.type = "EPIC")

GO_pos_PCA1 <- GO_pos_PCA1 %>% filter(FDR<0.05)

### calculate mean methylations and plot

# add direction of methylation for GO pathways

# calculate mean delta M for each gene

mean_change_df <- M_change %>% dplyr::select(37)


# Assuming you have a data frame named 'mean_change_df' with columns 'PH_vs_BH' and 'PM_vs_BM', and rownames as probe names

# get all unique genes
unique_gene <- unique(anno %>% 
                              filter(UCSC_RefGene_Name != "NA") %>% 
                              mutate(UCSC_RefGene_Name = sub(";.*", "", UCSC_RefGene_Name)) %>% 
                              pull(UCSC_RefGene_Name))

gene_cpgs <- list()

# create subset of probes annotated to individual genes

for (i in 1:length(unique_gene)) {
        y <- anno %>% 
                filter(UCSC_RefGene_Name == unique_gene[i]) %>% 
                pull(cpg)
        gene_cpgs[[unique_gene[i]]] <- y
        print(i)
}



# Initialize an empty data frame to store the results
results_df <- data.frame(
        gene = character(length(unique_gene)),
        mean_change_BM_vs_BH = numeric(length(unique_gene)),
        stringsAsFactors = FALSE
)

# Loop through unique genes and calculate mean change
for (i in 1:length(unique_gene)) {
        gene <- unique_gene[i]
        gene_probes <- gene_cpgs[[gene]]
        
        # Filter the mean_change_df to keep only the rows corresponding to the gene's probes
        gene_mean_change <- mean_change_df[rownames(mean_change_df) %in% gene_probes, ]
        
        # Calculate the mean change for each contrast and store it in the results data frame
        results_df[i, "gene"] <- gene
        results_df[i, "mean_change_BM_vs_BH"] <- sum(gene_mean_change) / length(gene_probes)
        
        
        print(i)
}

# Convert gene symbols to Entrez Gene IDs
entrezIDs <- mapIds(org.Hs.eg.db, keys = results_df$gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

as.data.frame(entrezIDs) %>% 
        rownames_to_column(var = "gene") %>% 
        merge(results_df,., by = "gene") %>% 
        filter(entrezIDs != "<NA>") -> change_df


pathway_dir_go = data.frame()



for (i in 1:length(go.sets$go.sets)) {
        pathway = names(go.sets$go.sets[i])
        pathway_genes = go.sets$go.sets[[pathway]]
        
        gene_mean_change <- change_df[change_df$entrezIDs %in% pathway_genes,]
        
        pathway_dir_go[i, "pathway"] <- pathway
        pathway_dir_go[i, "mean_change_BM_vs_BH"] <- sum(gene_mean_change$mean_change_BM_vs_BH) / length(pathway_genes)
        
        print(i)
}


# Positive PC1: Hypermethylated in MYO compared to MYO+INT at baseline
GO_pos_PCA1 %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(mean_change_BM_vs_BH = round(mean_change_BM_vs_BH, digits = 4)) %>% 
        filter(FDR < 0.05) %>% 
        arrange(-mean_change_BM_vs_BH) %>% 
        head(20)


# Positive PC1: Hypomethylated in MYO compared to MYO+INT at baseline
GO_pos_PCA1 %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(mean_change_BM_vs_BH = round(mean_change_BM_vs_BH, digits = 4)) %>% 
        filter(FDR < 0.05) %>% 
        arrange(mean_change_BM_vs_BH) %>% 
        head(20)


# Negative PC1: Hypermethylated in MYO compared to MYO+INT at baseline
GO_neg_PCA1 %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(mean_change_BM_vs_BH = round(mean_change_BM_vs_BH, digits = 4)) %>% 
        filter(FDR < 0.05) %>% 
        arrange(-mean_change_BM_vs_BH) %>% 
        head(20)


# Negative PC1: Hypomethylated in MYO compared to MYO+INT at baseline
GO_neg_PCA1 %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(mean_change_BM_vs_BH = round(mean_change_BM_vs_BH, digits = 4)) %>% 
        filter(FDR < 0.05) %>% 
        arrange(mean_change_BM_vs_BH) %>% 
        head(20)



#############################################################################################################
###############################   check main drivers of PC1, that are the most exercise responsive   ########
#############################################################################################################

# isolate the abs arranged list DMPs of PC1 drivers




as.data.frame(loadings) %>% 
        rownames_to_column(var = "cpg") %>% 
        # merge with significant DMPs post vs pre MYO
        merge(.,DMPs_PM_vs_BM, by = "cpg") %>% 
        mutate( loadings = round(loadings, digits = 5)) %>% 
        
        
        # merge with annotation dataset
        
        merge(., anno, by = "cpg") %>% 
        arrange(-abs(loadings), -abs(delta_M)) %>% 
        
        # extract top 10 list
        
        head(10) %>% 
        pull(cpg) -> top10_PC1_MYO_exercise_DMPs


# plot M_values for these genes

M_change %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% top10_PC1_MYO_exercise_DMPs[1]) %>% 
        dplyr::select(1:33) %>% 
        pivot_longer(cols = 2:33) %>% 
        mutate(sample = as.factor(substr(name, nchar(name) -0, nchar(name))),
               timepoint = as.factor(substr(name, nchar(name) -1, nchar(name)))) %>% 
        mutate(timepoint = as.factor(substr(timepoint, 1,1))) %>%
        mutate(timepoint = ifelse(timepoint == "B", "Baseline", "Post"),
               sample = ifelse(sample == "H", "MYO+INT", "MYO")) %>% 
        ggplot(aes(x = timepoint, y = value, fill = sample)) +
        geom_boxplot(position = position_dodge(width = 0.8)) +
        geom_point(position = position_dodge(width = 0.8), size = 3) +
        labs(fill = "Sample",
             y = "M_Value",
             x = "Timepoint", 
             title = (anno %>% filter(cpg == top10_PC1_MYO_exercise_DMPs[1]) %>% pull(UCSC_RefGene_Name))) + 
        theme_classic(base_size = 20)


# for loop top 10 PC1 drivers

for (i in 1:length(top10_PC1_MYO_exercise_DMPs)) {
        
        p1 <- M_change %>% 
                rownames_to_column(var = "cpg") %>% 
                filter(cpg %in% top10_PC1_MYO_exercise_DMPs[i]) %>% 
                dplyr::select(1:33) %>% 
                pivot_longer(cols = 2:33) %>% 
                mutate(sample = as.factor(substr(name, nchar(name) -0, nchar(name))),
                       timepoint = as.factor(substr(name, nchar(name) -1, nchar(name)))) %>% 
                mutate(timepoint = as.factor(substr(timepoint, 1,1))) %>%
                mutate(timepoint = ifelse(timepoint == "B", "Baseline", "Post"),
                       sample = ifelse(sample == "H", "MYO+INT", "MYO")) %>% 
                ggplot(aes(x = timepoint, y = value, fill = sample)) +
                geom_boxplot(position = position_dodge(width = 0.8)) +
                geom_point(position = position_dodge(width = 0.8), size = 3) +
                labs(fill = "Sample",
                     y = "M_Value",
                     x = "Timepoint", 
                     title = (anno %>% filter(cpg == top10_PC1_MYO_exercise_DMPs[i]) %>% pull(UCSC_RefGene_Name))) + 
                theme_classic(base_size = 20)+
                theme(legend.key.size = unit(1.5, "cm"))
        
        
        ggsave(p1, file = paste0("./", "PC1_figures/",anno %>% filter(cpg == top10_PC1_MYO_exercise_DMPs[i]) %>% pull(UCSC_RefGene_Name),".PNG"))
        
        print(i)
        
}


for (i in 1:length(top10_PC1_MYO_exercise_DMPs  )) {
        x <- anno %>% filter(cpg == top10_PC1_MYO_exercise_DMPs[i]) %>% pull(UCSC_RefGene_Name)
        print(x)
}    






# check if top10 exist in MYO+INT DMPs

DMPs_PH_vs_BH %>% 
        filter(cpg %in% top10_PC1_MYO_exercise_DMPs)



# repeat for MYO+INT significant genes

as.data.frame(loadings) %>% 
        rownames_to_column(var = "cpg") %>% 
        # merge with significant DMPs post vs pre MYO
        merge(.,DMPs_PH_vs_BH, by = "cpg") %>% 
        mutate( loadings = round(loadings, digits = 5)) %>% 
        
        
        # merge with annotation dataset
        
        merge(., anno, by = "cpg") %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        arrange(-abs(loadings), -abs(delta_M)) %>% 
        
        # extract top 10 list
        
        head(10) %>% 
        pull(cpg) -> top10_PC1_MYOINT_exercise_DMPs


# for loop top 10 PC1 drivers

for (i in 1:length(top10_PC1_MYOINT_exercise_DMPs)) {
        
        p1 <- M_change %>% 
                rownames_to_column(var = "cpg") %>% 
                filter(cpg %in% top10_PC1_MYOINT_exercise_DMPs[i]) %>% 
                dplyr::select(1:33) %>% 
                pivot_longer(cols = 2:33) %>% 
                mutate(sample = as.factor(substr(name, nchar(name) -0, nchar(name))),
                       timepoint = as.factor(substr(name, nchar(name) -1, nchar(name)))) %>% 
                mutate(timepoint = as.factor(substr(timepoint, 1,1))) %>%
                mutate(timepoint = ifelse(timepoint == "B", "Baseline", "Post"),
                       sample = ifelse(sample == "H", "MYO+INT", "MYO")) %>% 
                ggplot(aes(x = timepoint, y = value, fill = sample)) +
                geom_boxplot(position = position_dodge(width = 0.8)) +
                geom_point(position = position_dodge(width = 0.8), size = 3) +
                labs(fill = "Sample",
                     y = "M_Value",
                     x = "Timepoint", 
                     title = (anno %>% filter(cpg == top10_PC1_MYOINT_exercise_DMPs[i]) %>% pull(UCSC_RefGene_Name))) + 
                theme_classic(base_size = 20)+
                theme(legend.key.size = unit(1.5, "cm"))
        
        
        ggsave(p1, file = paste0("./", "PC1_figures/MYOINT/",anno %>% filter(cpg == top10_PC1_MYOINT_exercise_DMPs[i]) %>% pull(UCSC_RefGene_Name),".PNG"))
        
        print(i)
        
}


for (i in 1:length(top10_PC1_MYO_exercise_DMPs  )) {
        x <- anno %>% filter(cpg == top10_PC1_MYOINT_exercise_DMPs[i]) %>% pull(UCSC_RefGene_Name)
        print(x)
}    

DMPs_PM_vs_BM %>% 
        filter(cpg %in% top10_PC1_MYOINT_exercise_DMPs) %>% 
        merge(.,anno, by = "cpg")




###############################################################################################

###            PC1 plot

###############################################################################################

# PLot top 5 PC1 positive and negative, hypo and hyper GO pathways

# using FDR of GO GSEA and mean_M_change for size of circle and intensity of color

# get name of top 5 go pathways in all 4 groups

# Positive PC1: Hypermethylated in MYO compared to MYO+INT at baseline
GO_pos_PCA1 %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(mean_change_BM_vs_BH = round(mean_change_BM_vs_BH, digits = 4)) %>% 
        filter(FDR < 0.05) %>% 
        arrange(-mean_change_BM_vs_BH) %>% 
        head(5) %>% 
        pull(pathway) -> Positive_hyper_PC1


# Positive PC1: Hypomethylated in MYO compared to MYO+INT at baseline
GO_pos_PCA1 %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(mean_change_BM_vs_BH = round(mean_change_BM_vs_BH, digits = 4)) %>% 
        filter(FDR < 0.05) %>% 
        arrange(mean_change_BM_vs_BH) %>% 
        head(5) %>% 
        pull(pathway) -> Positive_hypo_PC1


# Negative PC1: Hypermethylated in MYO compared to MYO+INT at baseline
GO_neg_PCA1 %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(mean_change_BM_vs_BH = round(mean_change_BM_vs_BH, digits = 4)) %>% 
        filter(FDR < 0.05) %>% 
        arrange(-mean_change_BM_vs_BH) %>% 
        head(5) %>% 
        pull(pathway) -> Negative_hyper_PC1


# Negative PC1: Hypomethylated in MYO compared to MYO+INT at baseline
GO_neg_PCA1 %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(mean_change_BM_vs_BH = round(mean_change_BM_vs_BH, digits = 4)) %>% 
        filter(FDR < 0.05) %>% 
        arrange(mean_change_BM_vs_BH) %>% 
        head(5) %>% 
        pull(pathway) -> Negative_hypo_PC1


# combine df of top 5 go pathways

# create ordering vector

order <- rbind(GO_pos_PCA1 %>% 
          rownames_to_column(var = "pathway") %>% 
          merge(.,pathway_dir_go, by = "pathway") %>% 
          mutate(BM_vs_BH = round(mean_change_BM_vs_BH, digits = 4)) %>% 
          filter(pathway %in% Positive_hyper_PC1 | pathway %in% Positive_hypo_PC1) %>% 
          mutate(PC1 = "positive",
                 methylation = ifelse(pathway %in% Positive_hyper_PC1, "Hyper", "Hypo")) %>% 
          arrange(BM_vs_BH),
      
      
      GO_neg_PCA1 %>% 
          rownames_to_column(var = "pathway") %>% 
          merge(.,pathway_dir_go, by = "pathway") %>% 
          mutate(BM_vs_BH = round(mean_change_BM_vs_BH, digits = 4)) %>% 
          filter(pathway %in% Negative_hyper_PC1 | pathway %in% Negative_hypo_PC1) %>% 
          mutate(PC1 = "negative",
                 methylation = ifelse(pathway %in% Negative_hyper_PC1, "Hyper", "Hypo")) %>% 
          arrange(BM_vs_BH)) %>% 
    
    # remove GO code
    mutate(pathway = substr(pathway, 12,100)) %>% 
    pull(pathway) 



rbind(GO_pos_PCA1 %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(BM_vs_BH = round(mean_change_BM_vs_BH, digits = 4)) %>% 
        filter(pathway %in% Positive_hyper_PC1 | pathway %in% Positive_hypo_PC1) %>% 
        mutate(PC1 = "positive",
               methylation = ifelse(pathway %in% Positive_hyper_PC1, "Hyper", "Hypo")) %>% 
              arrange(BM_vs_BH),


        GO_neg_PCA1 %>% 
                rownames_to_column(var = "pathway") %>% 
                merge(.,pathway_dir_go, by = "pathway") %>% 
                mutate(BM_vs_BH = round(mean_change_BM_vs_BH, digits = 4)) %>% 
                filter(pathway %in% Negative_hyper_PC1 | pathway %in% Negative_hypo_PC1) %>% 
                mutate(PC1 = "negative",
                       methylation = ifelse(pathway %in% Negative_hyper_PC1, "Hyper", "Hypo")) %>% 
                arrange(BM_vs_BH)) %>% 
                
                # remove GO code
                mutate(pathway = substr(pathway, 12,100),
                       group = paste0(PC1, "_", methylation),
                       pathway = factor(pathway, levels = c(order))) %>% 
                mutate(BM_vs_BH = ifelse(BM_vs_BH < 0, 
                                         rescale(BM_vs_BH, to = c(-1,0), from = c(min(BM_vs_BH),0)),
                                         rescale(BM_vs_BH, to = c(0,1), from = c(0,max(BM_vs_BH))))) %>% 
                
                # calculate % differentially methylated probes of terms
                mutate(Percent_DMGs_of_GO_term = round((DE/N)*100), 0) %>% 
                ggplot(aes(y = pathway, x = group, color = (BM_vs_BH)))+
                #geom_point(aes( size = -log10(FDR)))+ # use % of probes in term to determine size of point
                geom_point(aes( size = Percent_DMGs_of_GO_term, alpha = -log10(FDR)))+
                scale_color_viridis()+        
                theme_classic(base_size = 20)+
                theme(axis.text.x = element_text(angle = 30, hjust = 1))+
                theme(panel.grid.major.y = element_line(linewidth=0.2, color="grey"))+
                geom_vline(xintercept = c(0.5, 1.5, 2.5, 3.5,4.5), color = "white", linewidth = 10)+
                scale_x_discrete(expand = c(0.12,0.12))+
                labs(color = "HYPO - HYPER")+
                theme(axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      #legend.position = "top",
                      #legend.text = element_text(size = 10, color = c("white","white","white","white","white", "Black")),
                      legend.title = element_text(size = 12, face = "bold", hjust = -0.5, vjust = 1),
                      axis.text.y = element_text(size = 12, face = "bold"))
        







###############################################################################################

###            dot plot exercise responisive genes in MYO vs. MYO+INT

###############################################################################################

# plot dot plot of exercise DMPs, wuth MYO on x-axis and MYO+INT on y-axis

# shared dots = black
# MYO = pink
# MYO+INT = grey


# create list of all significant probes

DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        pull(cpg) -> MYOINT

DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>%
        pull(cpg) -> MYO
        
all_dmps <- c(MYOINT, MYO)

shared_dmps <- all_dmps[duplicated(all_dmps)]

# merge exercise repsponsive DMPs

merge(DMPs_PH_vs_BH, DMPs_PM_vs_BM, by = "cpg") %>% 
        filter(cpg %in% all_dmps) %>% 
        # add significant in what group column
        mutate(signif = ifelse(p.value.x < 0.05, "MYOINT", "MYO")) %>% 
        mutate(signif = ifelse(cpg %in% shared_dmps, "Shared", signif)) %>% 
        merge(., anno, by = "cpg")%>% 
        mutate(UCSC_RefGene_Name = ifelse(UCSC_RefGene_Name == "NA", "", UCSC_RefGene_Name))-> exercise_responsive_dmps

# get labels of all points with delta_M more than ABS 1

exercise_responsive_dmps %>% 
        filter(abs(delta_M.x) > 1) -> x
exercise_responsive_dmps %>% 
        filter(abs(delta_M.y) > 1) -> y

rbind(x,y) -> labs
labs <- labs[!duplicated(labs$cpg),]

exercise_dmps_plot <- exercise_responsive_dmps %>% 
     #   head(10000) %>% 
        ggplot(aes(x = as.numeric(delta_M.y), y = as.numeric(delta_M.x), fill = signif))+
        geom_hline(yintercept = 0, alpha = 0.5)+
        geom_vline(xintercept = 0, alpha = 0.5)+
        geom_point(aes(color = signif)) +
        geom_label_repel(data = labs, aes(x = delta_M.y, y = delta_M.x, label = UCSC_RefGene_Name ) ,box.padding = 1.5, size = 3.5, point.padding = 1)+
        theme_classic(base_size = 20)+
        labs(x = "MYO delta_M",
             y = "MYO+INT delta_M")+
        scale_color_manual(values = c("#440154FF","#5DC863FF", "RED"))+
        scale_fill_manual(values = c("#B396B9","#B3D7B1", "Red"))+
        scale_x_continuous(n.breaks = 6)+
        scale_y_continuous(n.breaks = 6)

# exports interactive plot

library(htmlwidgets)
library(plotly)        

saveWidget(ggplotly(exercise_dmps_plot), file = "./Figures/exercise_responsive_DMPs.html")


# plot the DMPs with the highest delta m -> labs df

labs_flt <- labs %>% 
        mutate(absolute_delta = abs(delta_M.x)+abs(delta_M.y)) %>% 
        arrange(-absolute_delta) %>% 
        filter(UCSC_RefGene_Name != "") 

CPG =   labs_flt[5,1]    


p1 <- M_change[CPG,1:32] %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 1:32) %>% 
        mutate(Timepoint = ifelse(grepl("B", FP), "Baseline", "Post"),
               Sample = ifelse(grepl("H", FP), "MYO+INT", "MYO"),
               FP = paste0("FP",as.numeric(gsub("[^0-9]", "", FP)))) %>% 
        ggplot(aes(y = M_value, x = Timepoint, fill  = Timepoint))+
        geom_boxplot(width = 0.6)+
        geom_point(size = 3, alpha = 0.6)+
        geom_line(aes(group = FP), alpha = 0.6)+
        facet_grid(~Sample)+
        theme_classic(base_size = 20)+
        labs(title = paste(labs_flt %>% filter(cpg == CPG) %>% pull(UCSC_RefGene_Name), 
                           "  Signif in:", labs_flt %>% filter(cpg == CPG) %>% pull(signif)))


p2 <- M_change[CPG,1:32] %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 1:32) %>% 
        mutate(Timepoint = ifelse(grepl("B", FP), "Baseline", "Post"),
               Sample = ifelse(grepl("H", FP), "MYO+INT", "MYO"),
               FP = paste0("FP",as.numeric(gsub("[^0-9]", "", FP)))) %>% 
        pivot_wider(names_from = Sample, values_from = M_value) %>% 
        ggplot(aes(y = `MYO+INT`, x = MYO, color = Timepoint))+
        geom_point(size = 3)+
        theme_classic(base_size = 20)+
        labs(title = paste(labs_flt %>% filter(cpg == CPG) %>% pull(UCSC_RefGene_Name), 
                           "  Signif in:", labs_flt %>% filter(cpg == CPG) %>% pull(signif)))

ggarrange(p1,p2, ncol = 1, common.legend = TRUE, legend = "bottom")



###############################################################################################
##########      Exercise related DMPs in MYO                            #######################
##########      Promoter and island related DMPs                        #######################
###############################################################################################

# filter MYO dmps for CpG islands within promoters

# get hypermethylated probes
MYO_hyper <- merge(DMPs_PM_vs_BM, anno %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island"), 
        by = "cpg") %>% 
        filter(p.value < 0.05) %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        filter(delta_M > 0) %>% 
        pull(cpg)

MYO_hypo <- merge(DMPs_PM_vs_BM, anno %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island"), 
        by = "cpg") %>% 
        filter(p.value < 0.05) %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        filter(delta_M < 0) %>% 
        pull(cpg)
        
        
# get updated go pathways

go.sets <- go.gsets(species = "human", pkg.name=NULL, id.type = "EG", keep.evidence=FALSE)

# create lists of the different subsets og GO terms

go.BP <- names(go.sets$go.sets[go.sets$go.subs$BP])
go.CC <- names(go.sets$go.sets[go.sets$go.subs$CC])
go.MF <- names(go.sets$go.sets[go.sets$go.subs$MF])



# Run ORA

MYO_hyper_GO.all <- gsameth(sig.cpg = MYO_hyper,
                       all.cpg = rownames(beta), 
                       collection = go.sets$go.sets,
                       array.type = "EPIC")

MYO_hyper_GO.all %>% 
        filter(FDR < 0.05) %>% 
        mutate(percent_DE = DE/N) %>% 
        arrange(-percent_DE) %>% 
        rownames_to_column(var = "GO_term") %>% 
        mutate(Subset = ifelse(GO_term %in% go.BP, 'BP', ifelse(GO_term %in% go.CC, 'CC', 'MF'))) %>% 
        write_csv(., 'MYO_RT_HYPER_ISLAND_PROMOTER_GOTERMs.csv')

MYO_hyper_GO.BP <- gsameth(sig.cpg = MYO_hyper,
                            all.cpg = rownames(beta), 
                            collection = go.sets$go.sets[go.sets$go.subs$BP],
                            array.type = "EPIC") %>% 
        filter(FDR < 0.05) %>% 
        mutate(percent_DE = DE/N) %>% 
        arrange(-percent_DE)

MYO_hyper_GO.CC <- gsameth(sig.cpg = MYO_hyper,
                            all.cpg = rownames(beta), 
                            collection = go.sets$go.sets[go.sets$go.subs$CC],
                            array.type = "EPIC") %>% 
        filter(FDR < 0.05) %>% 
        mutate(percent_DE = DE/N) %>% 
        arrange(-percent_DE)

MYO_hyper_GO.MF <- gsameth(sig.cpg = MYO_hyper,
                            all.cpg = rownames(beta), 
                            collection = go.sets$go.sets[go.sets$go.subs$MF],
                            array.type = "EPIC") %>% 
        filter(FDR < 0.05) %>% 
        mutate(percent_DE = DE/N) %>% 
        arrange(-percent_DE)


MYO_hypo_GO.all <- gsameth(sig.cpg = MYO_hypo,
                            all.cpg = rownames(beta), 
                            collection = go.sets$go.sets,
                            array.type = "EPIC")

MYO_hypo_GO.BP <- gsameth(sig.cpg = MYO_hypo,
                           all.cpg = rownames(beta), 
                           collection = go.sets$go.sets[go.sets$go.subs$BP],
                           array.type = "EPIC") %>% 
        filter(FDR < 0.1) %>% 
        mutate(percent_DE = DE/N) %>% 
        arrange(-percent_DE)

MYO_hypo_GO.CC <- gsameth(sig.cpg = MYO_hypo,
                           all.cpg = rownames(beta), 
                           collection = go.sets$go.sets[go.sets$go.subs$CC],
                           array.type = "EPIC") %>% 
        filter(FDR < 0.1) %>% 
        mutate(percent_DE = DE/N) %>% 
        arrange(-percent_DE)

MYO_hypo_GO.MF <- gsameth(sig.cpg = MYO_hypo,
                           all.cpg = rownames(beta), 
                           collection = go.sets$go.sets[go.sets$go.subs$MF],
                           array.type = "EPIC") %>% 
        filter(FDR < 0.1) %>% 
        mutate(percent_DE = DE/N) %>% 
        arrange(-percent_DE)


# find tropomyosin GO term

go.sets$go.sets["GO:0005862 muscle thin filament tropomyosin"] %>% as.data.frame

anno[go.sets$go.sets[["GO:0005862 muscle thin filament tropomyosin"]],]


# Convert gene symbols to Entrez Gene IDs
entrezIDs <- mapIds(org.Hs.eg.db, keys = go.sets$go.sets[["GO:0005862 muscle thin filament tropomyosin"]], column ="SYMBOL", keytype = "ENTREZID", multiVals = "first")

TPM <- DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        filter(UCSC_RefGene_Name %in% c("TPM1","TPM2", "TPM3", "TPM4"))

M_change[TPM$cpg,]


####
#       Plot horizontal histogram of N probes + DMP within the pathway
####



MYO_hypo_plot <- rbind(MYO_hyper_GO.BP %>% 
              rownames_to_column(var = "GO_term") %>% 
              mutate(GO_term = substr(GO_term, 12, nchar(GO_term)),
                     N_DE = N-DE) %>% 
              pivot_longer(cols = c("DE", "N_DE")) %>% 
              mutate(GO_cat = "BP"),
      MYO_hyper_GO.CC %>% 
              rownames_to_column(var = "GO_term") %>% 
              mutate(GO_term = substr(GO_term, 12, nchar(GO_term)),
                     N_DE = N-DE) %>% 
              pivot_longer(cols = c("DE", "N_DE")) %>% 
              mutate(GO_cat = "CC")) %>% 
        rbind(., 
              MYO_hyper_GO.MF %>% 
                rownames_to_column(var = "GO_term") %>% 
                mutate(GO_term = substr(GO_term, 12, nchar(GO_term)),
                        N_DE = N-DE) %>% 
                pivot_longer(cols = c("DE", "N_DE")) %>%
                mutate(GO_cat = "MF")) %>% 
        ggplot(aes(x = value, y = GO_term, fill = factor(name, levels = rev(unique(name)))))+
        geom_bar(stat = "identity", width = 0.5)+
        facet_wrap(~GO_cat, nrow = 3, strip.position = "right")+
        labs(x = "N genes", 
             fill = "Genes in pathway")+
        theme_classic(base_size = 16)+
        theme(axis.title.y = element_blank(), 
              legend.position = "top", 
              axis.title.x = element_blank())+
        scale_x_continuous(expand = c(0,0))+
        scale_fill_manual(values = c("#440154FF","#5DC863FF"),
                          labels = c("All genes in GO term", "DMGs in GO term"))



p1_BP <- MYO_hyper_GO.BP %>% 
        rownames_to_column(var = "GO_term") %>% 
        mutate(GO_term = substr(GO_term, 12, nchar(GO_term)),
               N_DE = N-DE) %>% 
        pivot_longer(cols = c("DE", "N_DE")) %>% 
        ggplot(aes(x = value, y = GO_term, fill = factor(name, levels = rev(unique(name)))))+
                geom_bar(stat = "identity", width = 0.5)+
        labs(x = "N genes", 
             fill = "Genes in pathway")+
        theme_classic(base_size = 16)+
        theme(axis.title.y = element_blank(), 
              legend.position = "top", 
              axis.title.x = element_blank())+
        scale_x_continuous(expand = c(0,0))+
        scale_fill_manual(values = c("#440154FF","#5DC863FF"),
                          labels = c("All genes in GO term", "DMGs in GO term"))

p2_CC <- MYO_hyper_GO.CC %>% 
        rownames_to_column(var = "GO_term") %>% 
        mutate(GO_term = substr(GO_term, 12, nchar(GO_term)),
               N_DE = N-DE) %>% 
        pivot_longer(cols = c("DE", "N_DE")) %>% 
        ggplot(aes(x = value, y = GO_term, fill = factor(name, levels = rev(unique(name)))))+
        geom_bar(stat = "identity", width = 0.5)+
        labs(x = "N genes", 
             fill = "Genes in pathway")+
        theme_classic(base_size = 16)+
        theme(axis.title.y = element_blank(), 
              legend.position = "none",
              axis.title.x = element_blank())+
        scale_x_continuous(expand = c(0,0))+
        scale_fill_manual(values = c("#440154FF","#5DC863FF"),
                          labels = c("All genes in GO term", "DMGs in GO term"))

p3_MF <- MYO_hyper_GO.MF %>% 
        rownames_to_column(var = "GO_term") %>% 
        mutate(GO_term = substr(GO_term, 12, nchar(GO_term)),
               N_DE = N-DE) %>% 
        pivot_longer(cols = c("DE", "N_DE")) %>% 
        ggplot(aes(x = value, y = GO_term, fill = factor(name, levels = rev(unique(name)))))+
        geom_bar(stat = "identity", width = 0.5)+
        labs(x = "N genes", 
             fill = "Genes in pathway")+
        theme_classic(base_size = 16)+
        theme(axis.title.y = element_blank(), 
              legend.position = "none")+
        scale_x_continuous(expand = c(0,0))+
        scale_fill_manual(values = c("#440154FF","#5DC863FF"),
                          labels = c("All genes in GO term", "DMGs in GO term"))

plot_grid(p1_BP, p2_CC, p3_MF, ncol = 1)


####
#               Check ORA of MYO+INT island and promoter hypermethylation
####

# filter MYO+INT dmps for CpG islands within promoters

# get hypermethylated probes
MYOINT_hyper <- merge(DMPs_PH_vs_BH, anno %>% 
                           filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island"), 
                   by = "cpg") %>% 
        filter(p.value < 0.05) %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        filter(delta_M > 0) %>% 
        pull(cpg)


MYOINT_hyper_GO.all <- gsameth(sig.cpg = MYOINT_hyper,
                           all.cpg = rownames(beta), 
                           collection = go.sets$go.sets,
                           array.type = "EPIC")

MYOINT_hyper_GO.all %>% 
        filter(FDR < 0.05) %>% 
        mutate(percent_DE = DE/N) %>% 
        arrange(-percent_DE) %>% 
        rownames_to_column(var = "GO_term") %>% 
        mutate(Subset = ifelse(GO_term %in% go.BP, 'BP', ifelse(GO_term %in% go.CC, 'CC', 'MF'))) %>% 
        write_csv(., 'MYO+INT_RT_HYPER_ISLAND_PROMOTER_GOTERMs.csv')

# get hyp0methylated probes
MYOINT_hypo <- merge(DMPs_PH_vs_BH, anno %>% 
                              filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island"), 
                      by = "cpg") %>% 
        filter(p.value < 0.05) %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        filter(delta_M < 0) %>% 
        pull(cpg)


MYOINT_hypo_GO.all <- gsameth(sig.cpg = MYOINT_hypo,
                               all.cpg = rownames(beta), 
                               collection = go.sets$go.sets,
                               array.type = "EPIC")

MYOINT_hypo_GO.all %>% 
        filter(FDR < 0.1) %>% 
        mutate(percent_DE = DE/N) %>% 
        arrange(-percent_DE)


# check genes

anno %>% 
        filter(cpg %in% MYO_hyper) %>% 
        distinct(UCSC_RefGene_Name, .keep_all = TRUE)



###############################################################################################

###            correlations

###############################################################################################

# plot absoulte correlation between MYO vs. MYO+INT

M_change %>% 
        ggplot(aes(x = BH, y = BM)) +
        geom_point() +
        theme_classic(base_size = 20)+
        labs(x = "Baseline MYO+INT",
             y = "Baseline MYO")

M_change %>% 
        ggplot(aes(x = PH, y = PM)) +
        geom_point() +
        theme_classic(base_size = 20)+
        labs(x = "Post MYO+INT",
             y = "Post MYO")
        
baseline_cor <- cor(x = M_change [1:8], y = M_change[9:16])

# mean cor baseline

((baseline_cor[1,1]+baseline_cor[2,2]+baseline_cor[3,3]+baseline_cor[4,4]+baseline_cor[5,5]+baseline_cor[6,6]+baseline_cor[7,7]+baseline_cor[8,8])/8)^2

post_cor <- cor(x = M_change [17:24], y = M_change[25:32])

((post_cor[1,1]+post_cor[2,2]+post_cor[3,3]+post_cor[4,4]+post_cor[5,5]+post_cor[6,6]+post_cor[7,7]+post_cor[8,8])/8)^2

###
# plot correlation between top exercise responsive genes
###



M_change[rownames(M_change) %in% labs_flt$cpg,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Timepoint = ifelse(grepl("B", FP), "Baseline", "Post"),
               Sample = ifelse(grepl("H", FP), "MYO+INT", "MYO"),
               FP = paste0("FP",as.numeric(gsub("[^0-9]", "", FP)))) %>% 
        pivot_wider(names_from = Sample, values_from = M_value) %>% 
        ggplot(aes(y = `MYO+INT`, x = MYO, color = Timepoint))+
        geom_point(size = 3)+
        theme_classic(base_size = 20)


M_change[rownames(M_change) %in% labs_flt$cpg,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Timepoint = ifelse(grepl("B", FP), "Baseline", "Post"),
               Sample = ifelse(grepl("H", FP), "MYO+INT", "MYO"),
               FP = paste0("FP",as.numeric(gsub("[^0-9]", "", FP)))) %>% 
        pivot_wider(names_from = Sample, values_from = M_value) -> cor_df

# correlation matrix


cor(M_change[rownames(M_change) %in% labs_flt$cpg,c(1:8,17:24)], 
    M_change[rownames(M_change) %in% labs_flt$cpg,c(9:16,25:32)])

pheatmap::pheatmap(cor(M_change[rownames(M_change) %in% labs_flt$cpg,c(1:8,17:24)], 
                       M_change[rownames(M_change) %in% labs_flt$cpg,c(9:16,25:32)]), 
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers = TRUE, fontsize_number = 10, number_color = "black", viridis(10))

# correlation baseline

mean(0.67,0.78,0.78,0.61,0.68,0.64,0.5,0.49)

# correlation post

mean(0.61,0.56,0.75,0.76,0.74,0.71,0.68,0.67)

# variation matrix


var(M_change[rownames(M_change) %in% labs_flt$cpg,c(1:8,17:24)], 
    M_change[rownames(M_change) %in% labs_flt$cpg,c(9:16,25:32)])

pheatmap::pheatmap(var(M_change[rownames(M_change) %in% labs_flt$cpg,c(1:8,17:24)], 
                       M_change[rownames(M_change) %in% labs_flt$cpg,c(9:16,25:32)]), 
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers = TRUE, fontsize_number = 10, number_color = "black", viridis(10))





###
# plot correlation between top MYO+INT exercise responsive genes
###

# get list

ex.resp_myoint <- DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        arrange(-abs(delta_M)) %>% 
        head(100) %>% 
        pull(cpg)

# plot heatmap

pheatmap::pheatmap(cor(M_change[rownames(M_change) %in% ex.resp_myoint,c(1:8,17:24)], 
                       M_change[rownames(M_change) %in% ex.resp_myoint,c(9:16,25:32)]), 
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers = TRUE, fontsize_number = 10, number_color = "black", viridis(10))

# plot var

pheatmap::pheatmap(var(M_change[rownames(M_change) %in% ex.resp_myoint,c(1:8,17:24)], 
                       M_change[rownames(M_change) %in% ex.resp_myoint,c(9:16,25:32)]), 
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers = TRUE, fontsize_number = 10, number_color = "black", viridis(10))


###
# plot correlation between top MYO exercise responsive genes
###

# get list

ex.resp_myo <- DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        arrange(-abs(delta_M)) %>% 
        head(100) %>% 
        pull(cpg)

# plot heatmap

pheatmap::pheatmap(cor(M_change[rownames(M_change) %in% ex.resp_myo,c(1:8,17:24)], 
                       M_change[rownames(M_change) %in% ex.resp_myo,c(9:16,25:32)]), 
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers = TRUE, fontsize_number = 10, number_color = "black", viridis(10))

pheatmap::pheatmap(var(M_change[rownames(M_change) %in% ex.resp_myo,c(1:8,17:24)], 
                       M_change[rownames(M_change) %in% ex.resp_myo,c(9:16,25:32)]), 
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers = TRUE, fontsize_number = 10, number_color = "black", viridis(10))





###################################################################################################################
###################################################################################################################

###            plot hypo and hyper methylation of all genes relative to transcriptional start site

###################################################################################################################
###################################################################################################################


# workflow

# 1. get genomic position of all annotated probes
# 2. get transcriptional start site for the gene
# 3. calculate distance from tss
# 4. plot y axis = M_value, x axis = distance to TSS, color = annotated region 


###
# get TSS
### 


# add TSS to figure

BiocManager::install("biomaRt")

library(biomaRt)

listMarts()

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mart_attributes <- listAttributes(mart = ensembl)

region_data <- getBM(attributes = c("ensembl_gene_id", 
                                 "external_gene_name", 
                                 "start_position",
                                 "end_position",
                                 "5_utr_start",
                                 "5_utr_end",
                                 "3_utr_start",
                                 "3_utr_end",
                                 "transcription_start_site", 
                                 "chromosome_name", 
                                 "strand"),
                  mart = ensembl, verbose = TRUE)


# filter out duplicated Enseml IDs


region_data <- region_data %>%
        distinct(ensembl_gene_id, .keep_all = TRUE)


# merge annotation files and get the probeID, and position

anno_data <- region_data %>% 
        as.data.frame() %>% 
        unique() %>% 
        filter(external_gene_name != "") %>% 
        mutate(UCSC_RefGene_Name = external_gene_name) 

        
        merge(.,anno,by = "UCSC_RefGene_Name")
        pull(transcription_start_site) %>% mean()


illumina_anno <-readRDS("./EPIC.manifest.hg19.RDATA")


merged_anno <- illumina_anno %>% 
        dplyr::select(CpG_chrm, CpG_beg, CpG_end) %>% 
        rownames_to_column(var = "cpg") %>% 
        merge(anno, ., by = "cpg") %>% 
        left_join(., anno_data, by = "UCSC_RefGene_Name") %>% 
        distinct(cpg, .keep_all = TRUE)


anno_flt <- merged_anno %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        mutate(dist_cpg_to_TSS = transcription_start_site-CpG_beg,
               gene_body = ifelse(CpG_beg > start_position & CpG_beg < end_position, "Gene_Body", "")) 



M_change[,1:32] %>% 
        head(100000) %>% 
        rownames_to_column(var = "cpg") %>% 
        merge(., anno_flt, by = "cpg" ) %>% 
        dplyr::select(1:33, dist_cpg_to_TSS, gene_body, Regulatory_Feature_Group) %>% 
        filter(dist_cpg_to_TSS != "") %>% 
        mutate(TSS = ifelse(dist_cpg_to_TSS > -200 & dist_cpg_to_TSS < 0, "TSS",""), 
        
               TSS = ifelse(gene_body == "Gene_Body", "gene_body", TSS)) %>% 
        pull(TSS) %>% unique()
        pivot_longer(cols = 2:33, names_to = "FP") %>%
        mutate(Timepoint = ifelse(grepl("B", FP), "Baseline", "Post"),
               Sample = ifelse(grepl("H", FP), "MYO+INT", "MYO"),
               FP = paste0("FP",as.numeric(gsub("[^0-9]", "", FP)))) %>% 
        
        
        ggplot(aes(x = dist_cpg_to_TSS, y = value, color = Sample))+
        geom_point()+
        geom_vline(xintercept = 0)+
        scale_x_continuous(n.breaks = 20)


        
###################################################################################################
##########       Fit mixed effects model on aggregate promoter methylation     ####################
###################################################################################################
        
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

average_promoter.m %>% dim()

# create model matrix

promoter_model <- data.frame(ID = colnames(average_promoter.m[1:32]), 
                            row.names = colnames(average_promoter.m[1:32])) 

promoter_model$Time[grep('B', promoter_model$ID)] <- "Baseline"
promoter_model$Time[grep('P', promoter_model$ID)] <- 'Post'
promoter_model$Sample[grep('M', promoter_model$ID)] <- 'MYO'
promoter_model$Sample[grep('H', promoter_model$ID)] <- 'MYOINT'
promoter_model$FP <- paste0("FP_",(unlist(str_extract_all(promoter_model$ID, "\\d+"))))

promoter_model <- lapply(promoter_model, as.factor)

str(promoter_model)

# for loop that merges the annotation and methylation data and then runs mixed effects model

library(lme4)
library(lmerTest)
library(emmeans)
library(parallel)
library(pbapply)


model_results <- list()

for (i in 1:nrow(average_promoter.m)) {
        
        methylation_levels <- as.numeric(average_promoter.m[i,1:32])
        
        data <- cbind(methylation_levels, as.data.frame(promoter_model))
        
        model <- lmer(methylation_levels ~ Time * Sample + (1 | FP), data = data)
        
        model_results[[i]] <- model
        
        names(model_results)[i] <- rownames(average_promoter.m)[i]
        
        print(i)
        
}


summary(model_results[[2]])

# Obtain estimated marginal means
emm <- emmeans(model_results[[2]], ~ Time * Sample)

# Contrast between Baseline and Post for each Sample
contrast_time <- contrast(emm, method = "pairwise", by = "Sample")
print(contrast_time)

# Contrast between MYO and MYOINT for each Time
contrast_sample <- contrast(emm, method = "pairwise", by = "Time")
print(contrast_sample)        

# Interaction contrasts (simple main effects)
interaction_contrasts <- contrast(emm, interaction = "pairwise")
print(interaction_contrasts)


####
#      Extract emm and contrasts for all models
####

# Initialize a data frame to store contrast results
contrast_results <- data.frame(
        gene = character(),
        time_contrast_myo = numeric(),
        se_time_myo = numeric(),
        df_time_myo = numeric(),
        t_ratio_time_myo = numeric(),
        p_value_time_myo = numeric(),
        time_contrast_myoint = numeric(),
        se_time_myoint = numeric(),
        df_time_myoint = numeric(),
        t_ratio_time_myoint = numeric(),
        p_value_time_myoint = numeric(),
        sample_result_baseline = numeric(),
        se_sample_baseline = numeric(),
        df_sample_baseline = numeric(),
        t_ratio_sample_baseline = numeric(),
        p_value_sample_baseline = numeric(),
        sample_result_post = numeric(),
        se_sample_post = numeric(),
        df_sample_post = numeric(),
        t_ratio_sample_post = numeric(),
        p_value_sample_post = numeric(),
        estimate_interaction = numeric(),
        se_interaction = numeric(),
        df_interaction = numeric(),
        t_ratio_interaction = numeric(),
        p_value_interaction = numeric(),
        stringsAsFactors = FALSE
)


# Function to extract contrasts for a given model
extract_contrasts <- function(model, gene_name) {
        emm <- emmeans(model, ~ Time * Sample)
        # Contrast between Baseline and Post for each Sample
        contrast_time <- contrast(emm, method = "pairwise", by = "Sample")
        # Contrast between MYO and MYOINT for each Time
        contrast_sample <- contrast(emm, method = "pairwise", by = "Time")
        # Interaction contrasts (simple main effects)
        interaction_contrasts <- contrast(emm, interaction = "pairwise")
        # Extract results for time contrast (Baseline - Post) within each sample
        time_result_myo <- as.data.frame(contrast_time)[1, ]
        time_result_myoint <- as.data.frame(contrast_time)[2, ]
        # Extract results for sample contrast (MYO - MYOINT) within each time
        sample_result_baseline <- as.data.frame(contrast_sample)[1, ]
        sample_result_post <- as.data.frame(contrast_sample)[2, ]
        # Extract results for interaction contrast
        interaction_result <- as.data.frame(interaction_contrasts)[1, ]
        
        # Consolidate results into a single row for this gene
        result_row <- data.frame(
                gene = gene_name,
                time_contrast_myo = time_result_myo$estimate,
                se_time_myo = time_result_myo$SE,
                df_time_myo = time_result_myo$df,
                t_ratio_time_myo = time_result_myo$t.ratio,
                p_value_time_myo = time_result_myo$p.value,
                time_contrast_myoint = time_result_myoint$estimate,
                se_time_myoint = time_result_myoint$SE,
                df_time_myoint = time_result_myoint$df,
                t_ratio_time_myoint = time_result_myoint$t.ratio,
                p_value_time_myoint = time_result_myoint$p.value,
                sample_result_baseline = sample_result_baseline$estimate,
                se_sample_baseline = sample_result_baseline$SE,
                df_sample_baseline = sample_result_baseline$df,
                t_ratio_sample_baseline = sample_result_baseline$t.ratio,
                p_value_sample_baseline = sample_result_baseline$p.value,
                sample_result_post = sample_result_post$estimate,
                se_sample_post = sample_result_post$SE,
                df_sample_post = sample_result_post$df,
                t_ratio_sample_post = sample_result_post$t.ratio,
                p_value_sample_post = sample_result_post$p.value,
                estimate_interaction = interaction_result$estimate,
                se_interaction = interaction_result$SE,
                df_interaction = interaction_result$df,
                t_ratio_interaction = interaction_result$t.ratio,
                p_value_interaction = interaction_result$p.value,
                stringsAsFactors = FALSE
        )
        
        return(result_row)
}

# Number of cores to use for parallel processing
num_cores <- detectCores() - 1

# Use pblapply for parallel processing instead of mclapply to include progressbar and estimated time to completion
contrast_results <- pblapply(1:length(model_results), function(i) {
        gene_name <- names(model_results)[i]
        model <- model_results[[i]]
        extract_contrasts(model, gene_name)
}, cl = num_cores)

# Combine the list of data frames into a single data frame
contrast_results_df <- do.call(rbind, contrast_results)


######################################################################################################################
#######          save                                                                   ##############################
######################################################################################################################

write_csv(contrast_results_df, file = './aggregate_promoter_contrasts.csv')
write_csv(average_promoter.m, file = './average_promoter.meth.csv')

######################################################################################################################
#######          load data directly                                                     ##############################
######################################################################################################################

contrast_results_df <- read.csv('./aggregate_promoter_contrasts.csv')
average_promoter.m  <- read.csv('./average_promoter.meth.csv')




####
#       check for significant differential methylation of aggregate promoter methylation
####

contrast_results_df %>% 
        mutate(time_contrast_myo_FDR = p.adjust(p_value_time_myo)) %>% 
        arrange(time_contrast_myo_FDR) %>% 
        head(5)

contrast_results_df %>% 
        mutate(time_contrast_myoint_FDR = p.adjust(p_value_time_myoint),
               time_contrast_myo_FDR = p.adjust(p_value_time_myo),
               interaction_FDR = p.adjust(p_value_interaction)) %>% 
        arrange(abs(estimate_interaction)) %>% 
        rownames_to_column(var = "idx") %>% 
        filter(gene == "54762")
        head(20)


# one significant FDR < 0.05 = 54762 (GRAMD1C)

# plot GRAMD1C
emm_54762 <- emmeans(model_results[["54762"]], ~ Time * Sample)
a = 0.053
b = 0.053
c = 0.053
d = 0.053

# GRAMD1C <-
        data.frame(Time = summary(emm_54762)$Time,
                      Sample = summary(emm_54762)$Sample,
                      emmean = summary(emm_54762)$emmean,
                      lower.CL = summary(emm_54762)$lower.CL,
                      upper.CL = summary(emm_54762)$upper.CL) %>% 
        ggplot(aes(x = Time, y = emmean, fill = Sample))+
        geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.5, linewidth = 1.2)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.5), linewidth = 0.8)+
        scale_fill_manual(values = c("#B396B9","#B3D7B1"), labels = c("MYO", "MYO+INT"))+
        scale_y_continuous(expand = c(0,0), limits = c(0, 0.05))+
        theme_classic(base_size = 20)+
        labs(y = "Aggregate promoter methylation (\u03b2-value)", 
             title = "GRAMD1C")+
        geom_point(data = t(average_promoter.m["54762",1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(ID, "GRAMD1C" = "54762") %>% merge(., promoter_model, by = "ID"), 
                   aes(x= Time, y = GRAMD1C, fill = Sample), position = position_dodge(width = 0.5), size = 2)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 20, face = "bold", colour = "black"))+
        # Baseline comparison
                geom_segment(aes(x = 0.875, xend = 1.125, y = 0.78*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the left side of the bracket
                geom_segment(aes(x = 0.875, xend = 0.875, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the right side of the bracket
                geom_segment(aes(x = 1.125, xend = 1.125, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
                # Adding label for the bracket
                geom_text(aes(x = 1, y = 0.805*a, label = paste("p =", round(contrast_results_df %>% filter(gene == "54762") %>% pull(p_value_sample_baseline),3))), size = 5)+
        # Post comparison
                geom_segment(aes(x = 1.875, xend = 2.125, y = 0.78*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the left side of the bracket
                geom_segment(aes(x = 1.875, xend = 1.875, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the right side of the bracket
                geom_segment(aes(x = 2.125, xend = 2.125, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
                # Adding label for the bracket
                geom_text(aes(x = 2, y = 0.805*b, label = paste("p =", round(contrast_results_df %>% filter(gene == "54762") %>% pull(p_value_sample_post),4))), size = 5)+
        # MYO comparison
                geom_segment(aes(x = 0.875, xend = 1.875, y = 0.88*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the left side of the bracket
                geom_segment(aes(x = 0.875, xend = 0.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the right side of the bracket
                geom_segment(aes(x = 1.875, xend = 1.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
                # Adding label for the bracket
                geom_text(aes(x = 1.37, y = 0.905*c, label = paste("p =", round(contrast_results_df %>% filter(gene == "54762") %>% pull(p_value_time_myo),3))), size = 5)+
        # MYO+INT comparison
                geom_segment(aes(x = 1.125, xend = 2.125, y = 0.83*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the left side of the bracket
                geom_segment(aes(x = 1.125, xend = 1.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the right side of the bracket
                geom_segment(aes(x = 2.125, xend = 2.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
                # Adding label for the bracket
                geom_text(aes(x = 1.67, y = 0.855*d, label = paste("p =", round(contrast_results_df %>% filter(gene == "54762") %>% pull(p_value_time_myoint),8))), size = 5)
        

# check which genes were significant in both MYo and MYO+INT with a lower threshold

contrast_results_df %>% 
        filter(p_value_time_myo < 0.05 & p_value_time_myoint < 0.05) %>% 
        arrange(-abs(estimate_interaction)) 



# plot MYEOV
emm_26579 <- emmeans(model_results[["26579"]], ~ Time * Sample)

MYEOV <- data.frame(Time = summary(emm_26579)$Time,
                      Sample = summary(emm_26579)$Sample,
                      emmean = summary(emm_26579)$emmean,
                      lower.CL = summary(emm_26579)$lower.CL,
                      upper.CL = summary(emm_26579)$upper.CL) %>% 
        ggplot(aes(x = Time, y = emmean, fill = Sample))+
        geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.5, linewidth = 1.2)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.5), linewidth = 0.8)+
        scale_fill_manual(values = c("#B396B9","#B3D7B1"), labels = c("MYO", "MYO+INT"))+
        scale_y_continuous(expand = c(0,0), limits = c(0, 0.6))+
        theme_classic(base_size = 20)+
        labs(y = "Aggregate promoter methylation (\u03b2-value)", 
             title = "MYEOV")+
        geom_point(data = t(average_promoter.m["26579",1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(ID, "MYEOV" = "26579") %>% merge(., promoter_model, by = "ID"), 
                   aes(x= Time, y = MYEOV, fill = Sample), position = position_dodge(width = 0.5), size = 2)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 20, face = "bold", colour = "black"))


# plot MYH16
emm_84176 <- emmeans(model_results[["84176"]], ~ Time * Sample)

MYH16 <- data.frame(Time = summary(emm_84176)$Time,
                    Sample = summary(emm_84176)$Sample,
                    emmean = summary(emm_84176)$emmean,
                    lower.CL = summary(emm_84176)$lower.CL,
                    upper.CL = summary(emm_84176)$upper.CL) %>% 
        ggplot(aes(x = Time, y = emmean, fill = Sample))+
        geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.5, linewidth = 1.2)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.5), linewidth = 0.8)+
        scale_fill_manual(values = c("#B396B9","#B3D7B1"), labels = c("MYO", "MYO+INT"))+
        scale_y_continuous(expand = c(0,0), limits = c(0, 0.7))+
        theme_classic(base_size = 20)+
        labs(y = "Aggregate promoter methylation (\u03b2-value)", 
             title = "MYH16")+
        geom_point(data = t(average_promoter.m["84176",1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(ID, "MYH16" = "84176") %>% merge(., promoter_model, by = "ID"), 
                   aes(x= Time, y = MYH16, fill = Sample), position = position_dodge(width = 0.5), size = 2)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 20, face = "bold", colour = "black"))


# plot MYBPH
summary(model_results[["4608"]])

emm_4608 <- emmeans(model_results[["4608"]], ~ Time * Sample)

# MYBPH <- 
        
        data.frame(Time = summary(emm_4608)$Time,
                    Sample = summary(emm_4608)$Sample,
                    emmean = summary(emm_4608)$emmean,
                    lower.CL = summary(emm_4608)$lower.CL,
                    upper.CL = summary(emm_4608)$upper.CL) %>% 
        ggplot(aes(x = Time, y = emmean, fill = Sample))+
        geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.5, linewidth = 1.2)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.5), linewidth = 0.8)+
        scale_fill_manual(values = c("#B396B9","#B3D7B1"), labels = c("MYO", "MYO+INT"))+
        scale_y_continuous(expand = c(0,0), limits = c(0, 1))+
        theme_classic(base_size = 20)+
        labs(y = "Aggregate promoter methylation (\u03b2-value)", 
             title = "MYBPH")+
        geom_point(data = t(average_promoter.m["4608",1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(ID, "MYBPH" = "4608") %>% merge(., promoter_model, by = "ID"), 
                   aes(x= Time, y = MYBPH, fill = Sample), position = position_dodge(width = 0.5), size = 2)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 20, face = "bold", colour = "black"))+
        # Baseline comparison
                geom_segment(aes(x = 0.875, xend = 1.125, y = 0.78, yend = 0.78), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the left side of the bracket
                geom_segment(aes(x = 0.875, xend = 0.875, y = 0.76, yend = 0.78), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the right side of the bracket
                geom_segment(aes(x = 1.125, xend = 1.125, y = 0.76, yend = 0.78), linetype = "solid", linewidth = 0.6) +
                # Adding label for the bracket
                geom_text(aes(x = 1, y = 0.805, label = paste("p =", round(contrast_results_df %>% filter(gene == "4608") %>% pull(p_value_sample_baseline),3))), size = 5)+
        # Post comparison
                geom_segment(aes(x = 1.875, xend = 2.125, y = 0.78, yend = 0.78), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the left side of the bracket
                geom_segment(aes(x = 1.875, xend = 1.875, y = 0.76, yend = 0.78), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the right side of the bracket
                geom_segment(aes(x = 2.125, xend = 2.125, y = 0.76, yend = 0.78), linetype = "solid", linewidth = 0.6) +
                # Adding label for the bracket
                geom_text(aes(x = 2, y = 0.805, label = paste("p =", round(contrast_results_df %>% filter(gene == "4608") %>% pull(p_value_sample_post),3))), size = 5)+
        # MYO comparison
                geom_segment(aes(x = 0.875, xend = 1.875, y = 0.88, yend = 0.88), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the left side of the bracket
                geom_segment(aes(x = 0.875, xend = 0.875, y = 0.86, yend = 0.88), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the right side of the bracket
                geom_segment(aes(x = 1.875, xend = 1.875, y = 0.86, yend = 0.88), linetype = "solid", linewidth = 0.6) +
                # Adding label for the bracket
                geom_text(aes(x = 1.37, y = 0.905, label = paste("p =", round(contrast_results_df %>% filter(gene == "4608") %>% pull(p_value_time_myo),3))), size = 5)+
        # MYO+INT comparison
                geom_segment(aes(x = 1.125, xend = 2.125, y = 0.83, yend = 0.83), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the left side of the bracket
                geom_segment(aes(x = 1.125, xend = 1.125, y = 0.81, yend = 0.83), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the right side of the bracket
                geom_segment(aes(x = 2.125, xend = 2.125, y = 0.81, yend = 0.83), linetype = "solid", linewidth = 0.6) +
                # Adding label for the bracket
                geom_text(aes(x = 1.67, y = 0.855, label = paste("p =", round(contrast_results_df %>% filter(gene == "4608") %>% pull(p_value_time_myoint),3))), size = 5)


# most changes in MYO
        
contrast_results_df %>% 
        mutate(time_contrast_myo_FDR = p.adjust(p_value_time_myo)) %>% 
        arrange(-abs(time_contrast_myo)) %>% 
        head(5)        

# plot LINC00327 long intergenic non-protein coding RNA 327
summary(model_results[["100506697"]])

emm_100506697 <- emmeans(model_results[["100506697"]], ~ Time * Sample)
a = 1.05
b = 1.05
c = 1.08
d = 1.05

# LINC00327 <-
data.frame(Time = summary(emm_100506697)$Time,
           Sample = summary(emm_100506697)$Sample,
           emmean = summary(emm_100506697)$emmean,
           lower.CL = summary(emm_100506697)$lower.CL,
           upper.CL = summary(emm_100506697)$upper.CL) %>% 
        ggplot(aes(x = Time, y = emmean, fill = Sample))+
        geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.5, linewidth = 1.2)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.5), linewidth = 0.8)+
        scale_fill_manual(values = c("#B396B9","#B3D7B1"), labels = c("MYO", "MYO+INT"))+
        scale_y_continuous(expand = c(0,0), limits = c(0, 1))+
        theme_classic(base_size = 20)+
        labs(y = "Aggregate promoter methylation (\u03b2-value)", 
             title = "LINC00327")+
        geom_point(data = t(average_promoter.m["100506697",1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(ID, "LINC00327" = "100506697") %>% merge(., promoter_model, by = "ID"), 
                   aes(x= Time, y = LINC00327, fill = Sample), position = position_dodge(width = 0.5), size = 2)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 20, face = "bold", colour = "black"))+
        # # Baseline comparison
        # geom_segment(aes(x = 0.875, xend = 1.125, y = 0.78*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the left side of the bracket
        # geom_segment(aes(x = 0.875, xend = 0.875, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the right side of the bracket
        # geom_segment(aes(x = 1.125, xend = 1.125, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # # Adding label for the bracket
        # geom_text(aes(x = 1, y = 0.805*a, label = paste("p =", round(contrast_results_df %>% filter(gene == "100506697") %>% pull(p_value_sample_baseline),3))), size = 5)+
        # # Post comparison
        # geom_segment(aes(x = 1.875, xend = 2.125, y = 0.78*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the left side of the bracket
        # geom_segment(aes(x = 1.875, xend = 1.875, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the right side of the bracket
        # geom_segment(aes(x = 2.125, xend = 2.125, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # # Adding label for the bracket
        # geom_text(aes(x = 2, y = 0.805*b, label = paste("p =", round(contrast_results_df %>% filter(gene == "100506697") %>% pull(p_value_sample_post),4))), size = 5)+
        # MYO comparison
        geom_segment(aes(x = 0.875, xend = 1.875, y = 0.88*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 0.875, xend = 0.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 1.875, xend = 1.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 1.37, y = 0.905*c, label = paste("p =", round(contrast_results_df %>% filter(gene == "100506697") %>% pull(p_value_time_myo),3))), size = 5)
        # # MYO+INT comparison
        # geom_segment(aes(x = 1.125, xend = 2.125, y = 0.83*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the left side of the bracket
        # geom_segment(aes(x = 1.125, xend = 1.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the right side of the bracket
        # geom_segment(aes(x = 2.125, xend = 2.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Adding label for the bracket
        # geom_text(aes(x = 1.67, y = 0.855*d, label = paste("p =", round(contrast_results_df %>% filter(gene == "100506697") %>% pull(p_value_time_myoint),8))), size = 5)


  

# plot GPR65 G protein-coupled receptor 65 (GPR65) 
summary(model_results[["8477"]])

emm_8477 <- emmeans(model_results[["8477"]], ~ Time * Sample)
a = 0.75
b = 0.75
c = 0.72
d = 0.7

# LINC00327 <-
data.frame(Time = summary(emm_8477)$Time,
           Sample = summary(emm_8477)$Sample,
           emmean = summary(emm_8477)$emmean,
           lower.CL = summary(emm_8477)$lower.CL,
           upper.CL = summary(emm_8477)$upper.CL) %>% 
        ggplot(aes(x = Time, y = emmean, fill = Sample))+
        geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.5, linewidth = 1.2)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.5), linewidth = 0.8)+
        scale_fill_manual(values = c("#B396B9","#B3D7B1"), labels = c("MYO", "MYO+INT"))+
        scale_y_continuous(expand = c(0,0), limits = c(0, 0.7))+
        theme_classic(base_size = 20)+
        labs(y = "Aggregate promoter methylation (\u03b2-value)", 
             title = "GPR65")+
        geom_point(data = t(average_promoter.m["8477",1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(ID, "GPR65" = "8477") %>% merge(., promoter_model, by = "ID"), 
                   aes(x= Time, y = GPR65, fill = Sample), position = position_dodge(width = 0.5), size = 2)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 20, face = "bold", colour = "black"))+
        # Baseline comparison
        geom_segment(aes(x = 0.875, xend = 1.125, y = 0.78*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 0.875, xend = 0.875, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 1.125, xend = 1.125, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 1, y = 0.805*a, label = paste("p =", round(contrast_results_df %>% filter(gene == "8477") %>% pull(p_value_sample_baseline),3))), size = 5)+
        # Post comparison
        geom_segment(aes(x = 1.875, xend = 2.125, y = 0.78*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 1.875, xend = 1.875, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 2.125, xend = 2.125, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 2, y = 0.805*b, label = paste("p =", round(contrast_results_df %>% filter(gene == "8477") %>% pull(p_value_sample_post),7))), size = 5)+
        # MYO comparison
        geom_segment(aes(x = 0.875, xend = 1.875, y = 0.88*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 0.875, xend = 0.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 1.875, xend = 1.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 1.37, y = 0.905*c, label = paste("p =", round(contrast_results_df %>% filter(gene == "8477") %>% pull(p_value_time_myo),3))), size = 5)
        # # MYO+INT comparison
        # geom_segment(aes(x = 1.125, xend = 2.125, y = 0.83*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the left side of the bracket
        # geom_segment(aes(x = 1.125, xend = 1.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the right side of the bracket
        # geom_segment(aes(x = 2.125, xend = 2.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Adding label for the bracket
        # geom_text(aes(x = 1.67, y = 0.855*d, label = paste("p =", round(contrast_results_df %>% filter(gene == "8477") %>% pull(p_value_time_myoint),8))), size = 5)


# plot SATB2 SATB homeobox 2
emm_23314 <- emmeans(model_results[["23314"]], ~ Time * Sample)
a = 0.8
b = 0.92
c = 0.88
d = 1.1

# SATB2 <-
data.frame(Time = summary(emm_23314)$Time,
           Sample = summary(emm_23314)$Sample,
           emmean = summary(emm_23314)$emmean,
           lower.CL = summary(emm_23314)$lower.CL,
           upper.CL = summary(emm_23314)$upper.CL) %>% 
        ggplot(aes(x = Time, y = emmean, fill = Sample))+
        geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.5, linewidth = 1.2)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.5), linewidth = 0.8)+
        scale_fill_manual(values = c("#B396B9","#B3D7B1"), labels = c("MYO", "MYO+INT"))+
        scale_y_continuous(expand = c(0,0), limits = c(0, 0.85))+
        theme_classic(base_size = 20)+
        labs(y = "Aggregate promoter methylation (\u03b2-value)", 
             title = "SATB2")+
        geom_point(data = t(average_promoter.m["23314",1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(ID, "SATB2" = "23314") %>% merge(., promoter_model, by = "ID"), 
                   aes(x= Time, y = SATB2, fill = Sample), position = position_dodge(width = 0.5), size = 2)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 20, face = "bold", colour = "black"))+
        # # Baseline comparison
        # geom_segment(aes(x = 0.875, xend = 1.125, y = 0.78*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the left side of the bracket
        # geom_segment(aes(x = 0.875, xend = 0.875, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the right side of the bracket
        # geom_segment(aes(x = 1.125, xend = 1.125, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # # Adding label for the bracket
        # geom_text(aes(x = 1, y = 0.805*a, label = paste("p =", round(contrast_results_df %>% filter(gene == "23314") %>% pull(p_value_sample_baseline),3))), size = 5)+
        # Post comparison
        geom_segment(aes(x = 1.875, xend = 2.125, y = 0.78*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 1.875, xend = 1.875, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 2.125, xend = 2.125, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 2, y = 0.805*b, label = paste("p =", round(contrast_results_df %>% filter(gene == "23314") %>% pull(p_value_sample_post),4))), size = 5)+
        # MYO comparison
        geom_segment(aes(x = 0.875, xend = 1.875, y = 0.88*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 0.875, xend = 0.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 1.875, xend = 1.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 1.37, y = 0.905*c, label = paste("p =", round(contrast_results_df %>% filter(gene == "23314") %>% pull(p_value_time_myo),3))), size = 5)
        # # MYO+INT comparison
        # geom_segment(aes(x = 1.125, xend = 2.125, y = 0.83*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the left side of the bracket
        # geom_segment(aes(x = 1.125, xend = 1.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the right side of the bracket
        # geom_segment(aes(x = 2.125, xend = 2.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Adding label for the bracket
        # geom_text(aes(x = 1.67, y = 0.855*d, label = paste("p =", round(contrast_results_df %>% filter(gene == "23314") %>% pull(p_value_time_myoint),8))), size = 5)


# Identify genes with high baseline difference, which increases with RT. 

contrast_results_df %>%
        arrange(-abs(sample_result_baseline)) %>% 
        filter(abs(sample_result_baseline) > 0.1) %>% 
        filter(abs(sample_result_post) > abs(sample_result_baseline)) %>% 
        filter(p_value_time_myo < 0.05 | p_value_time_myoint < 0.05)
        mutate(diff = abs(sample_result_baseline)+abs(sample_result_post)) %>% 
        arrange(-abs(diff)) 

# plot MYPN
# Not interesting
        
gene <- "84665"
gene_name = "MYPN"
emm <- emmeans(model_results[[gene]], ~ Time * Sample)
a = 0.8
b = 0.92
c = 0.88
d = 1.1

# SATB2 <-
data.frame(Time = summary(emm)$Time,
           Sample = summary(emm)$Sample,
           emmean = summary(emm)$emmean,
           lower.CL = summary(emm)$lower.CL,
           upper.CL = summary(emm)$upper.CL) %>% 
        ggplot(aes(x = Time, y = emmean, fill = Sample))+
        geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.5, linewidth = 1.2)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.5), linewidth = 0.8)+
        scale_fill_manual(values = c("#B396B9","#B3D7B1"), labels = c("MYO", "MYO+INT"))+
        scale_y_continuous(expand = c(0,0), limits = c(0, 0.6))+
        theme_classic(base_size = 20)+
        labs(y = "Aggregate promoter methylation (\u03b2-value)", 
             title = gene_name)+
        geom_point(data = t(average_promoter.m[gene,1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(ID, gene_name = gene) %>% merge(., promoter_model, by = "ID"), 
                   aes(x= Time, y = gene_name, fill = Sample), position = position_dodge(width = 0.5), size = 2)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 20, face = "bold", colour = "black"))+
        # # Baseline comparison
        # geom_segment(aes(x = 0.875, xend = 1.125, y = 0.78*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the left side of the bracket
        # geom_segment(aes(x = 0.875, xend = 0.875, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the right side of the bracket
        # geom_segment(aes(x = 1.125, xend = 1.125, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # # Adding label for the bracket
        # geom_text(aes(x = 1, y = 0.805*a, label = paste("p =", round(contrast_results_df %>% filter(gene == "gene_name") %>% pull(p_value_sample_baseline),3))), size = 5)+
        # Post comparison
        geom_segment(aes(x = 1.875, xend = 2.125, y = 0.78*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 1.875, xend = 1.875, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 2.125, xend = 2.125, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 2, y = 0.805*b, label = paste("p =", round(contrast_results_df %>% filter(gene == "gene_name") %>% pull(p_value_sample_post),4))), size = 5)+
        # MYO comparison
        geom_segment(aes(x = 0.875, xend = 1.875, y = 0.88*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 0.875, xend = 0.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 1.875, xend = 1.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 1.37, y = 0.905*c, label = paste("p =", round(contrast_results_df %>% filter(gene == gene_name) %>% pull(p_value_time_myo),3))), size = 5)
        # # MYO+INT comparison
        # geom_segment(aes(x = 1.125, xend = 2.125, y = 0.83*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the left side of the bracket
        # geom_segment(aes(x = 1.125, xend = 1.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the right side of the bracket
        # geom_segment(aes(x = 2.125, xend = 2.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Adding label for the bracket
        # geom_text(aes

# plot MYPN
# Not interesting

gene <- "9119"
gene_name = "KRT75"
emm <- emmeans(model_results[[gene]], ~ Time * Sample)
a = 0.8
b = 0.92
c = 0.88
d = 1.1

# SATB2 <-
data.frame(Time = summary(emm)$Time,
           Sample = summary(emm)$Sample,
           emmean = summary(emm)$emmean,
           lower.CL = summary(emm)$lower.CL,
           upper.CL = summary(emm)$upper.CL) %>% 
        ggplot(aes(x = Time, y = emmean, fill = Sample))+
        geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.5, linewidth = 1.2)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.5), linewidth = 0.8)+
        scale_fill_manual(values = c("#B396B9","#B3D7B1"), labels = c("MYO", "MYO+INT"))+
        scale_y_continuous(expand = c(0,0), limits = c(0, 0.6))+
        theme_classic(base_size = 20)+
        labs(y = "Aggregate promoter methylation (\u03b2-value)", 
             title = gene_name)+
        geom_point(data = t(average_promoter.m[gene,1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(ID, gene_name = gene) %>% merge(., promoter_model, by = "ID"), 
                   aes(x= Time, y = gene_name, fill = Sample), position = position_dodge(width = 0.5), size = 2)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 20, face = "bold", colour = "black"))+
        # # Baseline comparison
        # geom_segment(aes(x = 0.875, xend = 1.125, y = 0.78*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the left side of the bracket
        # geom_segment(aes(x = 0.875, xend = 0.875, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the right side of the bracket
        # geom_segment(aes(x = 1.125, xend = 1.125, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # # Adding label for the bracket
        # geom_text(aes(x = 1, y = 0.805*a, label = paste("p =", round(contrast_results_df %>% filter(gene == "gene_name") %>% pull(p_value_sample_baseline),3))), size = 5)+
        # Post comparison
        geom_segment(aes(x = 1.875, xend = 2.125, y = 0.78*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 1.875, xend = 1.875, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 2.125, xend = 2.125, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 2, y = 0.805*b, label = paste("p =", round(contrast_results_df %>% filter(gene == "gene_name") %>% pull(p_value_sample_post),4))), size = 5)+
        # MYO comparison
        geom_segment(aes(x = 0.875, xend = 1.875, y = 0.88*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 0.875, xend = 0.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 1.875, xend = 1.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 1.37, y = 0.905*c, label = paste("p =", round(contrast_results_df %>% filter(gene == gene_name) %>% pull(p_value_time_myo),3))), size = 5)
        # # MYO+INT comparison
        # geom_segment(aes(x = 1.125, xend = 2.125, y = 0.83*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the left side of the bracket
        # geom_segment(aes(x = 1.125, xend = 1.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Vertical segment on the right side of the bracket
        # geom_segment(aes(x = 2.125, xend = 2.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # # Adding label for the bracket
        # geom_text(aes(x = 1.67, y = 0.855*d, label = paste("p =", round(contrast_results_df %>% filter(gene == gene_name %>% pull(p_value_time_myoint),8))), size = 5)





###################################################################################################################
##################       dot plot of exercise responsive genes on aggregate methylation.         ##################
###################################################################################################################   

signif_promoter_exercise <- contrast_results_df %>% 
        filter(p_value_time_myo < 0.05 | p_value_time_myoint < 0.05) %>% 
        mutate(signif = ifelse(p_value_time_myo < 0.05, "MYO", 
                               ifelse(p_value_time_myoint < 0.05, "MYOINT", "none"))) %>% 
        mutate(signif = as.factor(ifelse(p_value_time_myo < 0.05 & p_value_time_myoint < 0.05, "BOTH", signif))) %>% 
        merge(., average_promoter.m[,39:40] %>% rownames_to_column(var = "gene"), by = "gene") 
        
        
# Convert gene symbols to Entrez Gene IDs
signif_promoter_exercise$entrezIDs <- mapIds(org.Hs.eg.db, keys = signif_promoter_exercise$gene, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

shared_labs <- signif_promoter_exercise %>% 
        dplyr::select(28:30) %>% 
        filter(entrezIDs != "MYBPH") %>% 
        filter(entrezIDs != "NACA4P") %>% 
        filter(entrezIDs != "VARS1",
               entrezIDs != "MYEOV") 

Labelled_genes <- c("LINC00327", "GPR174", "SLC2A9", "OR1J2", "LINC01603", 
                    "ACSM1", "BLK", "LINC00613", "TAAR5", "CCDC168", "OR5A1", "LINC02877", 
                    "LARP1", "OR51S1", "GNG2", 
                    "NPPB", "OTOA", "LINC00332", "GPR65", "TFEC", "ROS1", "SATB2","ZNF430", 
                    "MYBPH", "NEUROD6", "MS4A13", "H1-8", "MYH16", "MYEOV", "RTN4", "THPO", "GSTM2", "CLCC1", "GRAMD1C")

Labelled_genes2 <- c("LINC00327", "GPR174", "SLC2A9", "OR1J2", 
                    "ACSM1", "LINC00613", "TAAR5", "CCDC168", "LINC02877", 
                    "LARP1", "OR51S1", "GNG2", 
                    "NPPB", "LINC00332", "GPR65", "TFEC", "ROS1", "SATB2", 
                    "MYBPH", "NEUROD6", "MYH16", "MYEOV", "THPO")

labs <- signif_promoter_exercise %>% 
        dplyr::select(28:30) %>% 
        filter(entrezIDs %in% Labelled_genes2) %>% 
        mutate(entrezIDs = factor(entrezIDs, levels = Labelled_genes2)) %>% 
        arrange(entrezIDs)



agg_plot <- signif_promoter_exercise %>% 
        ggplot(aes(x = PM_vs_BM, y = PH_vs_BH, fill = signif))+
        geom_hline(yintercept = 0, alpha = 0.5)+
        geom_vline(xintercept = 0, alpha = 0.5)+
        geom_point(aes(color = signif), size = 2, shape = 21) +
        #geom_label_repel(aes(label = entrezIDs ) ,box.padding = 0.3, size = 3.5, point.padding = 0.5, force = 1.6, min.segment.length = 0.1, max.time = 2, max.overlaps = 12, show.legend = FALSE)+
        theme_classic(base_size = 25)+
        labs(x = "MYO",
             y = "MYO+INT", fill = "Significant", color = "Significant")+
        scale_color_manual(values = c("RED",  "#440154FF","#5DC863FF"), labels = c("Both", "MYO", "MYO+INT"))+
        scale_fill_manual(values = c("Red","#B396B9" ,"#B3D7B1"), labels = c("Both", "MYO", "MYO+INT"))+
        scale_x_continuous(n.breaks = 5, limits = c(-0.1,0.08), expand = c(0,0.0011))+
        scale_y_continuous(n.breaks = 5, limits = c(-0.1,0.08), expand = c(0,0.0005))+
        geom_label_repel(data = labs, aes(x = PM_vs_BM, y = PH_vs_BH,label = entrezIDs ),inherit.aes = FALSE,
                         box.padding = 0.3, size = 5, point.padding = 0.5, force = 1.6, min.segment.length = 0.1, max.time = 2, max.overlaps = 12, 
                         fill = c(rep("#B396B9",4), rep("#B3D7B1",8)  ,rep("#B396B9",6)  ,rep("red",5)), show.legend = FALSE)+
        theme(legend.position = "top", legend.key.spacing.x = unit(4, "mm"))+
        guides(color = guide_legend(override.aes = list(size = 6)))+
        theme(axis.text.r = element_text(color = 'black'),
                panel.background = element_rect(fill='transparent'), #transparent panel bg
                plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                panel.grid.major = element_blank(), #remove major gridlines
                panel.grid.minor = element_blank(), #remove minor gridlines
                legend.background = element_rect(fill='transparent'), #transparent legend bg
                legend.box.background = element_rect(fill='transparent')) #transparent legend panel


ggsave(plot = agg_plot, bg = "transparent", path = "Figures/", filename = "MACS_agg_plot_transparent.png") 

# save plot 

ggsave(plot = agg_plot, path = "Figures/Aggregate_promoter_methylation/", filename = "exercise_responsive_aggregate_promoter_methylation.png")


# count and calculate % of shared promoter genes and unique

signif_promoter_exercise %>%  filter(signif == 'MYO') %>% nrow() /signif_promoter_exercise %>% nrow()
signif_promoter_exercise %>%  filter(signif == 'MYOINT') %>% nrow() /signif_promoter_exercise %>% nrow()
signif_promoter_exercise %>%  filter(signif == 'BOTH') %>% nrow() /signif_promoter_exercise %>% nrow()

# individual gene exporation

# plot NEUROD6
# Not interesting

gene <- "63974"
gene_name = "NEUROD6"
emm <- emmeans(model_results[[gene]], ~ Time * Sample)
a = 0.8
b = 0.92
c = 0.88
d = 1.1

# SATB2 <-
data.frame(Time = summary(emm)$Time,
           Sample = summary(emm)$Sample,
           emmean = summary(emm)$emmean,
           lower.CL = summary(emm)$lower.CL,
           upper.CL = summary(emm)$upper.CL) %>% 
        ggplot(aes(x = Time, y = emmean, fill = Sample))+
        geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.5, linewidth = 1.2)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.5), linewidth = 0.8)+
        scale_fill_manual(values = c("#B396B9","#B3D7B1"), labels = c("MYO", "MYO+INT"))+
        scale_y_continuous(expand = c(0,0), limits = c(0, 1))+
        theme_classic(base_size = 20)+
        labs(y = "Aggregate promoter methylation (\u03b2-value)", 
             title = gene_name)+
        geom_point(data = t(average_promoter.m[gene,1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(ID, gene_name = gene) %>% merge(., promoter_model, by = "ID"), 
                   aes(x= Time, y = gene_name, fill = Sample), position = position_dodge(width = 0.5), size = 2)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 20, face = "bold", colour = "black"))+
        # Baseline comparison
        geom_segment(aes(x = 0.875, xend = 1.125, y = 0.78*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 0.875, xend = 0.875, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 1.125, xend = 1.125, y = 0.76*a, yend = 0.78*a), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 1, y = 0.805*a, label = paste("p =", round(contrast_results_df %>% filter(gene == "gene_name") %>% pull(p_value_sample_baseline),3))), size = 5)+
        # Post comparison
        geom_segment(aes(x = 1.875, xend = 2.125, y = 0.78*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 1.875, xend = 1.875, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 2.125, xend = 2.125, y = 0.76*b, yend = 0.78*b), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 2, y = 0.805*b, label = paste("p =", round(contrast_results_df %>% filter(gene == "gene_name") %>% pull(p_value_sample_post),4))), size = 5)+
        # MYO comparison
        geom_segment(aes(x = 0.875, xend = 1.875, y = 0.88*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 0.875, xend = 0.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 1.875, xend = 1.875, y = 0.86*c, yend = 0.88*c), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 1.37, y = 0.905*c, label = paste("p =", round(contrast_results_df %>% filter(gene == gene_name) %>% pull(p_value_time_myo),3))), size = 5)+
        # MYO+INT comparison
        geom_segment(aes(x = 1.125, xend = 2.125, y = 0.83*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the left side of the bracket
        geom_segment(aes(x = 1.125, xend = 1.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # Vertical segment on the right side of the bracket
        geom_segment(aes(x = 2.125, xend = 2.125, y = 0.81*d, yend = 0.83*d), linetype = "solid", linewidth = 0.6) +
        # Adding label for the bracket
        geom_text(aes(x = 1.67, y = 0.855*d, label = paste("p =", round(contrast_results_df %>% filter(gene == gene_name) %>% pull(p_value_time_myoint),8))), size = 5)


####
#       Write a loop to plot and save all labelled exercise responsive aggregate methylation genes
####

Labelled_genes <- c("LINC00327", "GPR174", "SLC2A9", "OR1J2", "LINC01603", 
                       "ACSM1", "BLK", "LINC00613", "TAAR5", "CCDC168", "OR5A1", "LINC02877", 
                       "LARP1", "OR51S1", "GNG2", 
                       "NPPB", "OTOA", "LINC00332", "GPR65", "TFEC", "ROS1", "SATB2","ZNF430", 
                       "MYBPH", "NEUROD6", "MS4A13", "H1-8", "MYH16", "MYEOV", "RTN4", "THPO", "GSTM2", "CLCC1", "GRAMD1C"
                       )

i = 32
for (i in 1:length(Labelled_genes)) {
        
        gene_name = Labelled_genes[i]
        # Convert gene symbols to Entrez Gene IDs
        Gene <- mapIds(org.Hs.eg.db, keys = gene_name, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
        
        emm <- emmeans(model_results[[Gene]], ~ Time * Sample)
        
        upper = round(max(t(average_promoter.m[Gene,1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(gene_name = 2))*1.3, 2)
        
        a = max(t(average_promoter.m[Gene,1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(gene_name = 2))*1.025
        b = max(t(average_promoter.m[Gene,1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(gene_name = 2))*1.025
        c = max(t(average_promoter.m[Gene,1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(gene_name = 2))*1.025
        d = max(t(average_promoter.m[Gene,1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(gene_name = 2))*1.025
        
        
      plot <- data.frame(Time = summary(emm)$Time,
                   Sample = summary(emm)$Sample,
                   emmean = summary(emm)$emmean,
                   lower.CL = summary(emm)$lower.CL,
                   upper.CL = summary(emm)$upper.CL) %>% 
                ggplot(aes(x = Time, y = emmean, fill = Sample))+
                geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.5, linewidth = 1.2)+
                geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.5), linewidth = 0.8)+
                scale_fill_manual(values = c("#B396B9","#B3D7B1"), labels = c("MYO", "MYO+INT"))+
                scale_y_continuous(expand = c(0,0), limits = c(0, upper))+
                theme_classic(base_size = 20)+
                labs(y = "Aggregate promoter methylation (\u03b2-value)", 
                     title = gene_name)+
                geom_point(data = t(average_promoter.m[Gene,1:32]) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>%  dplyr::select(ID, gene_name = 2) %>% merge(., promoter_model, by = "ID"), 
                           aes(x= Time, y = gene_name, fill = Sample), position = position_dodge(width = 0.5), size = 2)+
                theme(axis.title.x = element_blank(),
                      axis.text.x = element_text(size = 20, face = "bold", colour = "black"))+
                # Baseline comparison
                geom_segment(aes(x = 0.875, xend = 1.125, y = 1.02*a, yend = 1.02*a), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the left side of the bracket
                geom_segment(aes(x = 0.875, xend = 0.875, y = 1*a, yend = 1.02*a), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the right side of the bracket
                geom_segment(aes(x = 1.125, xend = 1.125, y = 1*a, yend = 1.02*a), linetype = "solid", linewidth = 0.6) +
                # Adding label for the bracket
                geom_text(aes(x = 1, y = 1.045*a, label = paste("p =", round(contrast_results_df %>% filter(gene == Gene) %>% pull(p_value_sample_baseline),5))), size = 5)+
                # Post comparison
                geom_segment(aes(x = 1.875, xend = 2.125, y = 1.02*b, yend = 1.02*b), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the left side of the bracket
                geom_segment(aes(x = 1.875, xend = 1.875, y = 1*b, yend = 1.02*b), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the right side of the bracket
                geom_segment(aes(x = 2.125, xend = 2.125, y = 1*b, yend = 1.02*b), linetype = "solid", linewidth = 0.6) +
                # Adding label for the bracket
                geom_text(aes(x = 2, y = 1.045*b, label = paste("p =", round(contrast_results_df %>% filter(gene == Gene) %>% pull(p_value_sample_post),5))), size = 5)+
                # MYO comparison
                geom_segment(aes(x = 0.875, xend = 1.875, y = 1.12*c, yend = 1.12*c), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the left side of the bracket
                geom_segment(aes(x = 0.875, xend = 0.875, y = 1.1*c, yend = 1.12*c), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the right side of the bracket
                geom_segment(aes(x = 1.875, xend = 1.875, y = 1.1*c, yend = 1.12*c), linetype = "solid", linewidth = 0.6) +
                # Adding label for the bracket
                geom_text(aes(x = 1.37, y = 1.145*c, label = paste("p =", round(contrast_results_df %>% filter(gene == Gene) %>% pull(p_value_time_myo),5))), size = 5)+
                # MYO+INT comparison
                geom_segment(aes(x = 1.125, xend = 2.125, y = 1.07*d, yend = 1.07*d), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the left side of the bracket
                geom_segment(aes(x = 1.125, xend = 1.125, y = 1.05*d, yend = 1.07*d), linetype = "solid", linewidth = 0.6) +
                # Vertical segment on the right side of the bracket
                geom_segment(aes(x = 2.125, xend = 2.125, y = 1.05*d, yend = 1.07*d), linetype = "solid", linewidth = 0.6) +
                # Adding label for the bracket
                geom_text(aes(x = 1.67, y = 1.095*d, label = paste("p =", round(contrast_results_df %>% filter(gene == Gene) %>% pull(p_value_time_myoint),5))), size = 5)
        
      ggsave(plot = plot, path = "Figures/Aggregate_promoter_methylation/Labelled_genes/", filename = paste0(gene_name, ".png"))
}




#############################################################################################################################
##################       plot all probes from the most ingeresting gene of aggregate promoter genes        ##################
#############################################################################################################################   




_________________________________________________

# find probes in genes

_________________________________________________

Illumina_anno <- Illumina_anno %>% 
        mutate("cpg" = IlmnID) 

x <- Illumina_anno %>% 
        mutate(UCSC_RefGene_Name = str_split(UCSC_RefGene_Name, ";", simplify = TRUE)[, 1]) 

x$cpg = x$Name

m_change_df <- M_change[39:40] %>% 
        rownames_to_column(var = "cpg") %>% 
        merge(., x, by = "cpg") %>% 
        dplyr::select(1:3, UCSC_RefGene_Name)


gene_name <- "LINC00327"       # write gene name

Probes <- anno %>% 
        filter(UCSC_RefGene_Name == gene_name) 

x = c("cg25431166", "cg11866473") # extra probes if misidentified

Chromosome <- merge(Probes, Illumina_anno %>% mutate(cpg = IlmnID), by = "cpg") %>% 
        dplyr::select(CHR) %>% 
        unique()

m_change_df %>% 
        filter(UCSC_RefGene_Name == gene_name) %>%
        merge(.,Illumina_anno, by = "cpg") %>% 
        dplyr::select(cpg, PH_vs_BH, PM_vs_BM, MAPINFO, Relation_to_UCSC_CpG_Island, CHR) %>% 
        pivot_longer(names_to = "contrast", values_to = "mean_change", cols = 2:3) %>% 
        arrange(desc(MAPINFO)) %>% 
        ggplot(aes(x = MAPINFO, y = mean_change, color = contrast, group = contrast))+
        geom_hline(yintercept = 0)+
        geom_point(aes(shape = Relation_to_UCSC_CpG_Island), size = 3)+
        geom_smooth(method = "loess", alpha = 0.5, se = FALSE)+
        labs(title = paste(gene_name,";", "N Probes =", nrow(Probes), ";","CHR", Chromosome),
             y = "mean M-value change",
             x = "nucleotide")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90))+
        scale_x_continuous(n.breaks = 20)



####
#       Find TSS LNC00327
####



BiocManager::install("biomaRt")

library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

tss_data <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "transcription_start_site", "chromosome_name", "strand"),
                  mart = ensembl)

tss_data %>% 
        as.data.frame() %>% 
        filter(external_gene_name == gene_name) %>% 
        pull(transcription_start_site) %>% mean()

# calculate LINC00327 distance to tss

Illumina_anno %>% 
        filter(UCSC_RefGene_Name == "LINC00327") %>% 
        mutate(TSS_distance = Start_hg38-tss_data %>% 
                       as.data.frame() %>% 
                       filter(external_gene_name == gene_name) %>% 
                       pull(transcription_start_site) %>% mean()) %>% 
        arrange(TSS_distance)


# plot boxplot of LINC00327

LINC00327 <- Illumina_anno %>% 
        filter(UCSC_RefGene_Name == "LINC00327") %>% 
        mutate(TSS_distance = Start_hg38-tss_data %>% 
                       as.data.frame() %>% 
                       filter(external_gene_name == gene_name) %>% 
                       pull(transcription_start_site) %>% mean()) %>% 
        arrange(TSS_distance) %>% 
        merge(., M_change %>% rownames_to_column(var = "cpg"), by = "cpg") %>% 
        pivot_longer(cols = 55:86, names_to = "Sample", values_to = "M_val") %>% 
        dplyr::select(1, TSS_distance, Sample, M_val) %>% 
        mutate(Time = "na", 
               FP = "na")


LINC00327$Time[grep('BM', LINC00327$Sample)] <- "Baseline_MYO"
LINC00327$Time[grep('PM', LINC00327$Sample)] <- 'Post_MYO'
LINC00327$Time[grep('BH', LINC00327$Sample)] <- "Baseline_MYOINT"
LINC00327$Time[grep('PH', LINC00327$Sample)] <- 'Post_MYOINT'
# LINC00327$FP <- paste0("FP_",(unlist(str_extract_all(LINC00327$Sample, "\\d+"))))
# LINC00327$Sample[grep('M', LINC00327$Sample)] <- 'MYO'
# LINC00327$Sample[grep('H', LINC00327$Sample)] <- 'MYOINT'


# plot

c1 <- LINC00327 %>% 
        mutate(Time = factor(Time, levels = c("Baseline_MYO", "Post_MYO", "Baseline_MYOINT", "Post_MYOINT"))) %>% 
        filter(cpg == "cg19409956") %>% 
        ggplot(aes(x = Time, y = M_val, fill = Time))+
        geom_boxplot(outliers = FALSE, linewidth = 1.5) +
        geom_point(show.legend = FALSE, size =3) + 
        scale_y_continuous(limits = c(0,4), expand = c(0,0))+
        scale_fill_manual(values = c( "#440154FF","#B396B9","#5DC863FF","#B3D7B1"), labels = c("MYO (Baseline)", "MYO (Post)", "MYO+INT (Baseline)", "MYO+INT (Post)"))+
        theme_classic(base_size = 20) +
        theme(axis.title.x = element_blank(), 
              axis.text.x = element_blank(), 
              axis.title.y = element_blank(), 
              panel.grid.major.y = element_line(), 
              legend.position = "none")+
        theme(axis.text.r = element_text(color = 'black'),
              panel.background = element_rect(fill='transparent'), #transparent panel bg
              plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
              panel.grid.major = element_blank(), #remove major gridlines
              panel.grid.minor = element_blank(), #remove minor gridlines
              legend.background = element_rect(fill='transparent'), #transparent legend bg
              legend.box.background = element_rect(fill='transparent'))

c2 <- LINC00327 %>% 
        mutate(Time = factor(Time, levels = c("Baseline_MYO", "Post_MYO", "Baseline_MYOINT", "Post_MYOINT"))) %>% 
        filter(cpg == "cg22576398") %>% 
        ggplot(aes(x = Time, y = M_val, fill = Time))+
        geom_boxplot(outliers = FALSE, linewidth = 1.5) +
        geom_point(show.legend = FALSE, size =3) + 
        scale_y_continuous(limits = c(0,4), expand = c(0,0))+
        scale_fill_manual(values = c( "#440154FF","#B396B9","#5DC863FF","#B3D7B1"), labels = c("MYO (Baseline)", "MYO (Post)", "MYO+INT (Baseline)", "MYO+INT (Post)"))+
        theme_classic(base_size = 20) +
        theme(axis.title.x = element_blank(), 
              axis.text.x = element_blank(), 
              axis.title.y = element_blank(), 
              panel.grid.major.y = element_line(), 
              legend.position = "none", 
              axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())+
        theme(axis.text.r = element_text(color = 'black'),
              panel.background = element_rect(fill='transparent'), #transparent panel bg
              plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
              panel.grid.major = element_blank(), #remove major gridlines
              panel.grid.minor = element_blank(), #remove minor gridlines
              legend.background = element_rect(fill='transparent'), #transparent legend bg
              legend.box.background = element_rect(fill='transparent'))

c3 <- LINC00327 %>% 
        mutate(Time = factor(Time, levels = c("Baseline_MYO", "Post_MYO", "Baseline_MYOINT", "Post_MYOINT"))) %>% 
        filter(cpg == "cg15896012") %>% 
        ggplot(aes(x = Time, y = M_val, fill = Time))+
        geom_boxplot(outliers = FALSE, linewidth = 1.5) +
        geom_point(show.legend = FALSE, size =3) + 
        scale_y_continuous(limits = c(0,4), expand = c(0,0))+
        scale_fill_manual(values = c( "#440154FF","#B396B9","#5DC863FF","#B3D7B1"), labels = c("MYO (Baseline)", "MYO (Post)", "MYO+INT (Baseline)", "MYO+INT (Post)"))+
        theme_classic(base_size = 20) +
        theme(axis.title.x = element_blank(), 
              axis.text.x = element_blank(), 
              axis.title.y = element_blank(), 
              panel.grid.major.y = element_line(), 
              legend.position = "none", 
              axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())+
        theme(axis.text.r = element_text(color = 'black'),
              panel.background = element_rect(fill='transparent'), #transparent panel bg
              plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
              panel.grid.major = element_blank(), #remove major gridlines
              panel.grid.minor = element_blank(), #remove minor gridlines
              legend.background = element_rect(fill='transparent'), #transparent legend bg
              legend.box.background = element_rect(fill='transparent'))

c4 <- LINC00327 %>% 
        mutate(Time = factor(Time, levels = c("Baseline_MYO", "Post_MYO", "Baseline_MYOINT", "Post_MYOINT"))) %>% 
        filter(cpg == "cg21995219") %>% 
        ggplot(aes(x = Time, y = M_val, fill = Time))+
        geom_boxplot(outliers = FALSE, linewidth = 1.5) +
        geom_point(show.legend = FALSE, size =3) + 
        scale_y_continuous(limits = c(0,4), expand = c(0,0))+
        scale_fill_manual(values = c( "#440154FF","#B396B9","#5DC863FF","#B3D7B1"), labels = c("MYO (Baseline)", "MYO (Post)", "MYO+INT (Baseline)", "MYO+INT (Post)"))+
        theme_classic(base_size = 20) +
        theme(axis.title.x = element_blank(), 
              axis.text.x = element_blank(), 
              axis.title.y = element_blank(), 
              panel.grid.major.y = element_line(), 
              legend.position = "none", 
              axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())+
        theme(axis.text.r = element_text(color = 'black'),
              panel.background = element_rect(fill='transparent'), #transparent panel bg
              plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
              panel.grid.major = element_blank(), #remove major gridlines
              panel.grid.minor = element_blank(), #remove minor gridlines
              legend.background = element_rect(fill='transparent'), #transparent legend bg
              legend.box.background = element_rect(fill='transparent'))


linc00327 <- plot_grid(c1, c2, c3, c4, nrow = 1)+
        theme(axis.text.r = element_text(color = 'black'),
              panel.background = element_rect(fill='transparent'), #transparent panel bg
              plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
              panel.grid.major = element_blank(), #remove major gridlines
              panel.grid.minor = element_blank(), #remove minor gridlines
              legend.background = element_rect(fill='transparent'), #transparent legend bg
              legend.box.background = element_rect(fill='transparent')) #transparent legend panel


ggsave(plot = linc00327, bg = "transparent", path = "Figures/", filename = "MACS_linc00327_transparent.png") 


###################################################################################################################
##################       plot up and down regulated GO pathways in bubble plot     ################################
###################################################################################################################        

# get updated go pathways

go.sets <- go.gsets(species = "human", pkg.name=NULL, id.type = "EG", keep.evidence=FALSE)


# run GSEA on hypo and hyper, MYO and MYO+INT

# subset GO sets

go.sets$go.sets[go.sets$go.subs$BP]


# 

GO.BP_MYOINT_hypo <- gsameth(sig.cpg = DMPs_PH_vs_BH %>% 
                               filter(p.value < 0.05,
                                      delta_M < 0) %>% 
                               pull(cpg),
                       all.cpg = rownames(beta), 
                       collection = go.sets$go.sets[go.sets$go.subs$BP],
                       array.type = "EPIC")

GO.BP_MYOINT_hyper <- gsameth(sig.cpg = DMPs_PH_vs_BH %>% 
                                  filter(p.value < 0.05,
                                         delta_M > 0) %>% 
                                  pull(cpg),
                          all.cpg = rownames(beta), 
                          collection = go.sets$go.sets[go.sets$go.subs$BP],
                          array.type = "EPIC")

GO.BP_MYO_hypo <- gsameth(sig.cpg = DMPs_PM_vs_BM %>% 
                                  filter(p.value < 0.05,
                                         delta_M < 0) %>% 
                                  pull(cpg),
                          all.cpg = rownames(beta), 
                          collection = go.sets$go.sets[go.sets$go.subs$BP],
                          array.type = "EPIC")

GO.BP_MYO_hyper <- gsameth(sig.cpg = DMPs_PM_vs_BM %>% 
                                   filter(p.value < 0.05,
                                          delta_M > 0) %>% 
                                   pull(cpg),
                           all.cpg = rownames(beta), 
                           collection = go.sets$go.sets[go.sets$go.subs$BP],
                           array.type = "EPIC")


###################################################################################################################
##################       REACTOME GSEA.                                            ################################
###################################################################################################################        


BiocManager::install("signatureSearch")
library(signatureSearch)

# Convert gene symbols to Entrez Gene IDs
entrezIDs <- mapIds(org.Hs.eg.db, keys = anno$UCSC_RefGene_Name, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

anno$entrezIDs <- entrezIDs

library(ReactomePA)


# over representaion analysis with reactome

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


# visualize pathway without fold change
viewPathway(
        "Muscle contraction",
        organism = "human",
        readable = TRUE,
        foldChange = NULL,
        keyType = "ENTREZID",
        layout = "kk"
)


####
#       add diff methylation to reactome plots
####

# calculate diff meth per gene

unique(MYO_hyper_eids) %>% length()




####
#       Get reactome DB
####

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


###################################################################################################################
##################       density plot of differences                 ##############################################
###################################################################################################################        



# density plot of differences in M_value at the same timepoint
        

M_change[,37:38] %>% 
        pivot_longer(cols = 1:2, names_to = "Comparison", values_to = "M_difference") %>% 
        ggplot(aes(x = M_difference, group = Comparison, color = Comparison, fill = Comparison))+
        geom_density(linewidth = 1)+
        scale_y_sqrt()+
        facet_wrap(~Comparison, ncol = 1, nrow = 2)+
        theme_classic(base_size = 20)+
        theme(strip.text = element_blank())+
        scale_color_manual(values = c("#440154FF","#5DC863FF"))+
        scale_fill_manual(values = c("#B396B9","#B3D7B1"))+
        labs(title = "Difference between MYO+INT vs. MYO at Baseline and Post RT",
             x = "MYO vs. MYO+INT")+
        geom_vline(xintercept = 0)



# blant altman plot

library(BlandAltmanLeh)
library(ggExtra)

        
# BH vs. BM
bland.altman.plot(M_change[,33], M_change[,34], graph.sys = "ggplot2", conf.int = 0.95)

# stats

bland.altman.stats(M_change[,33], M_change[,34])

# same plot, just with distribution histograms

ggMarginal(bland.altman.plot(M_change[,33], M_change[,34], graph.sys = "ggplot2", type = "histogram", size=4))


# PH vs. PM

ggMarginal(bland.altman.plot(M_change[,35], M_change[,36], graph.sys = "ggplot2", type = "histogram", size=4))

# stats

bland.altman.stats(M_change[,35], M_change[,36])
