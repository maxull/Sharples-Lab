############################################################################################
###                                                                                      ###
###  MACS: PCA analysis                                                                  ###
###                                                                                      ###
############################################################################################


#This script contains the majority of the analysis code for the master thesis of Max Ullrich



library(cowplot)
library(ggrepel)
library(scales)
library(tidyverse)
library(ggplot2)
library(stringr)
library(missMethyl)
library(ChAMP)
library(pathview)
library(viridis)
library(ENmix)
library(AnnotationDbi)
library(org.Hs.eg.db)

setwd("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/")

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
        as.data.frame() %>% 
        filter(.<0.05) 

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

loadings <- pca.out$rotation[,1]

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

go.sets <- go.gsets(species = "human", pkg.name=NULL, id.type = "eg", keep.evidence=FALSE)

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

exercise_responsive_dmps %>% 
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


        




































