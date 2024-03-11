############################################################################################
###                                                                                      ###
###  Redoing analysis with manual DMPs                                                   ###
###                                                                                      ###
############################################################################################


#This script contains the majority of the analysis code for the master thesis of Max Ullrich


library(cowplot)
library(dendextend)
library(ggpubr)
library(ggrepel)
library(scales)
library(tidyverse)
library(ggplot2)
library(stringr)
library(missMethyl)
library(ChAMP)
library(pathview)
library(FedData)
library(lme4)
library(lmerTest)
library(ggvenn)
library(viridis)
library(ENmix)
library(gageData)

setwd("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/")

beta <- readRDS(file = "beta.RDATA")

########################################################################

### re-do DMPs manually

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

# calculate paired sample t-test

# Extract relevant columns from M_change
BH <- as.matrix(M_change[,1:8])
BM <- as.matrix(M_change[,9:16])
PH <- as.matrix(M_change[,17:24])
PM <- as.matrix(M_change[,25:32])

# Pre-allocate memory for DMPs_PM_vs_PH
DMPs_BM_vs_BH <- data.frame(
        cpg = rownames(M_change),
        t = numeric(length(rownames(M_change))),
        p.value = numeric(length(rownames(M_change))),
        conf.int_0.025 = numeric(length(rownames(M_change))),
        conf.int_0.975 = numeric(length(rownames(M_change))),
        delta_M = numeric(length(rownames(M_change)))
)

DMPs_PM_vs_PH <- data.frame(
        cpg = rownames(M_change),
        t = numeric(length(rownames(M_change))),
        p.value = numeric(length(rownames(M_change))),
        conf.int_0.025 = numeric(length(rownames(M_change))),
        conf.int_0.975 = numeric(length(rownames(M_change))),
        delta_M = numeric(length(rownames(M_change)))
)

DMPs_PM_vs_BM <- data.frame(
        cpg = rownames(M_change),
        t = numeric(length(rownames(M_change))),
        p.value = numeric(length(rownames(M_change))),
        conf.int_0.025 = numeric(length(rownames(M_change))),
        conf.int_0.975 = numeric(length(rownames(M_change))),
        delta_M = numeric(length(rownames(M_change)))
)

DMPs_PH_vs_BH <- data.frame(
        cpg = rownames(M_change),
        t = numeric(length(rownames(M_change))),
        p.value = numeric(length(rownames(M_change))),
        conf.int_0.025 = numeric(length(rownames(M_change))),
        conf.int_0.975 = numeric(length(rownames(M_change))),
        delta_M = numeric(length(rownames(M_change)))
)


# Loop through rows and perform t.test

for (i in 1:804857) {
        test <- t.test(x = BM[i,], y = BH[i,], paired = TRUE)
        DMPs_BM_vs_BH[i, c("t", "p.value", "conf.int_0.025", "conf.int_0.975", "delta_M")] <- 
                c(test$statistic, test$p.value, test$conf.int[1], test$conf.int[2], test$estimate)
        
        print(i)
}

for (i in 1:804857) {
        test <- t.test(x = PM[i,], y = PH[i,], paired = TRUE)
        DMPs_PM_vs_PH[i, c("t", "p.value", "conf.int_0.025", "conf.int_0.975", "delta_M")] <- 
                c(test$statistic, test$p.value, test$conf.int[1], test$conf.int[2], test$estimate)
        
        print(i)
}

for (i in 1:804857) {
        test <- t.test(x = PM[i,], y = BM[i,], paired = TRUE)
        DMPs_PM_vs_BM[i, c("t", "p.value", "conf.int_0.025", "conf.int_0.975", "delta_M")] <- 
                c(test$statistic, test$p.value, test$conf.int[1], test$conf.int[2], test$estimate)
        
        print(i)
}

for (i in 1:804857) {
        test <- t.test(x = PH[i,], y = BH[i,], paired = TRUE)
        DMPs_PH_vs_BH[i, c("t", "p.value", "conf.int_0.025", "conf.int_0.975", "delta_M")] <- 
                c(test$statistic, test$p.value, test$conf.int[1], test$conf.int[2], test$estimate)
        
        print(i)
}


# performed these loops/calculations on confocal pc


DMPs_BM_vs_BH <- readRDS("DMPs_BM_vs_BH.RDATA")
DMPs_PH_vs_BH <- readRDS("DMPs_PH_vs_BH.RDATA")
DMPs_PM_vs_BM <- readRDS("DMPs_PM_vs_BM.RDATA")
DMPs_PM_vs_PH <- readRDS("DMPs_PM_vs_PH.RDATA")

# add adjusted p value for within time comparison
library(stats)

DMPs_BM_vs_BH$adj.p.val <- p.adjust(p = DMPs_BM_vs_BH$p.value, method = "fdr", n = length(rownames(DMPs_BM_vs_BH)))
DMPs_PM_vs_PH$adj.p.val <- p.adjust(p = DMPs_PM_vs_PH$p.value, method = "fdr", n = length(rownames(DMPs_PM_vs_PH)))        

p.adjust(p = DMPs_PM_vs_PH$p.value, method = "bonferroni", n = length(rownames(DMPs_PM_vs_PH))) %>% 
        as.data.frame() %>% 
        filter(.<0.05) 

# filter all dmp lists for p.value <=0.05


DMPs_BM_vs_BH <- DMPs_BM_vs_BH %>% 
        filter(adj.p.val < 0.05)

DMPs_PM_vs_PH <- DMPs_PM_vs_PH %>% 
        filter(adj.p.val < 0.05) 

DMPs_PM_vs_BM <- DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05)

DMPs_PH_vs_BH <- DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05)


#############################################################################################

### # create horizontal barchart of DMPs in different positions

########################################################################


# baseline myonuclei vs. homogenate



# annotate and dd methylation direction

baseline_dmps <- DMPs_BM_vs_BH %>% 
        as.data.frame() %>% 
        merge(.,anno, by = "cpg") %>% 
        mutate("dir" = as.factor(ifelse(delta_M <0, "Negative_skew", "Positive_skew")))



dfh_1 <- dfh %>% 
        filter(condition %in% c("BH", "BM"))        

change_b_vals_baseline <- beta[,rownames(dfh_1)] %>% 
        as.data.frame() %>% 
        mutate(average_homo = (`1BH`+ `2BH`+ `4BH`+ `5BH`+ `6BH`+ `7BH`+ `8BH`+`12BH`)/8,
               average_myo = (`1BM`+ `2BM`+ `4BM`+ `5BM`+ `6BM`+ `7BM`+ `8BM`+`12BM`)/8,
               change = average_myo-average_homo) %>% 
        
        dplyr::select(change) %>% 
        rownames_to_column(var = "cpg")


n_count <- merge(baseline_dmps, change_b_vals_baseline, by.x = "cpg") %>% 
        filter(dir == "Negative_skew") %>% 
        mutate(Relation_to_Island = factor(Relation_to_Island, levels = c("S_Shelf","N_Shelf","Island", "S_Shore", "N_Shore", "OpenSea"))) %>% 
        summarise(dplyr::count(.,Relation_to_Island)) %>% 
        as.data.frame()

p_count <- merge(baseline_dmps, change_b_vals_baseline, by.x = "cpg") %>% 
        filter(dir == "Positive_skew") %>%
        mutate(Relation_to_Island = factor(Relation_to_Island, levels = c("S_Shelf","N_Shelf","Island", "S_Shore", "N_Shore", "OpenSea"))) %>% 
        summarise(dplyr::count(.,Relation_to_Island)) %>% 
        as.data.frame() 

df <- data.frame(
        category = n_count[6:1,1],
        neg_skew = (n_count[6:1,2])*-1,
        pos_skew = p_count[6:1,2]) %>% 
        mutate(total = (neg_skew)*-1+pos_skew)

p1 <- ggplot(df, aes(x = category)) +
        geom_bar(aes(y = neg_skew, fill = "neg_skew"), stat = "identity") +
        geom_bar(aes(y = pos_skew, fill = "pos_skew"), stat = "identity") +
        scale_fill_manual(values = c("#453781FF", "#DCE319FF")) + # set colors
        geom_label(aes(label = as.integer(neg_skew)*-1, x = category, y = -15000), size = 4, alpha = 0.8, nudge_y = -3000)+
        geom_label(aes(label = as.integer(pos_skew), x = category, y = 15000), size = 4, alpha = 0.8, nudge_y = 3000)+
        coord_flip()+
        theme_classic(base_size = 20)+
        theme(legend.title = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_text(hjust = 0.58))+
        scale_y_continuous(n.breaks = 6)+
        labs(y = "Number of DMPs")


# find number of significant cpgs that were altered more than <0.01 B < 0.1 B < 0.2 B <


change_b_vals_baseline %>% 
        mutate(change = abs(change)) %>% # change all numbers to positive
        merge(baseline_dmps,.,by.x = "cpg") %>% 
        # find values in ranges
        mutate(skew = ifelse(change <0.01, "<0.01", 
                             ifelse(change >0.01 & change < 0.1, "0.01<0.1", 
                                    ifelse(change > 0.1 & change < 0.2, "0.1<0.2", 
                                           ifelse(change > 0.2 & change < 0.3, "0.2<0.3",">0.3"))))) %>% 
        summarise(dplyr::count(.,skew)) -> df1

df2 <- data.frame(category = factor(df1[,1], levels = c("<0.01", "0.01<0.1", "0.1<0.2", "0.2<0.3", ">0.3")),
                  count = df1[,2])

library(scales)
p3 <- ggplot(data = df2, aes(x = category))+
        geom_bar(aes(y = count), stat = "identity")+
        geom_text(aes(label = count, y = count), vjust = -0.5)+
        theme_classic(base_size = 20)+
        labs(x = "B-value skew")+
        theme(axis.title.y = element_blank())+
        scale_y_continuous(trans = "log2", expand  = c(0, 1), n.breaks = 9)


# repeat for post samples

post_dmps <- DMPs_PM_vs_PH %>% 
        as.data.frame() %>% 
        merge(.,anno, by = "cpg") %>% 
        mutate("dir" = as.factor(ifelse(delta_M <0, "Negative_skew", "Positive_skew")))


# merge with newer annotation


dfh_2 <- dfh %>% 
        filter(condition %in% c("PH", "PM"))        

change_b_vals_post <- beta[,rownames(dfh_2)] %>% 
        as.data.frame() %>% 
        mutate(average_homo = (`1PH`+ `2PH`+ `4PH`+ `5PH`+ `6PH`+ `7PH`+ `8PH`+`12PH`)/8,
               average_myo = (`1PM`+ `2PM`+ `4PM`+ `5PM`+ `6PM`+ `7PM`+ `8PM`+`12PM`)/8,
               change = average_myo-average_homo) %>% 
        dplyr::select(change) %>% 
        rownames_to_column(var = "cpg")


n_count2 <- merge(post_dmps, change_b_vals_post, by.x = "cpg") %>% 
        filter(dir == "Negative_skew") %>%
        mutate(Relation_to_Island = factor(Relation_to_Island, levels = c("S_Shelf","N_Shelf","Island", "S_Shore", "N_Shore", "OpenSea"))) %>% 
        summarise(dplyr::count(.,Relation_to_Island)) %>% 
        as.data.frame()

p_count2 <- merge(post_dmps, change_b_vals_post, by.x = "cpg") %>% 
        filter(dir == "Positive_skew") %>% 
        mutate(Relation_to_Island = factor(Relation_to_Island, levels = c("S_Shelf","N_Shelf","Island", "S_Shore", "N_Shore", "OpenSea"))) %>% 
        summarise(dplyr::count(.,Relation_to_Island)) %>% 
        as.data.frame() 

df2 <- data.frame(
        category = n_count2[6:1,1],
        neg_skew = (n_count2[6:1,2])*-1,
        pos_skew = p_count2[6:1,2])

p2 <- ggplot(df2, aes(x = category)) +
        geom_bar(aes(y = neg_skew, fill = "neg_skew"), stat = "identity") +
        geom_bar(aes(y = pos_skew, fill = "pos_skew"), stat = "identity") +
        scale_fill_manual(values = c("#453781FF", "#DCE319FF")) + # set colors
        geom_label(aes(label = as.integer(neg_skew)*-1, x = category, y = -17000), size = 4, alpha = 0.8, nudge_y = -3000)+
        geom_label(aes(label = as.integer(pos_skew), x = category, y = 16000), size = 4, alpha = 0.8, nudge_y = 3000)+
        coord_flip()+
        theme_classic(base_size = 20)+
        theme(legend.title = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_text(hjust = 0.58))+
        scale_y_continuous(labels = label_comma(),
                           n.breaks = 7)+
        labs(y = "Number of DMPs")


# find number of significant cpgs that were altered more than <0.01 B < 0.1 B < 0.2 B <


change_b_vals_post %>% 
        mutate(change = abs(change)) %>% # change all numbers to positive
        merge(post_dmps,.,by.x = "cpg") %>% 
        # find values in ranges
        mutate(skew = ifelse(change <0.01, "<0.01", 
                             ifelse(change >0.01 & change < 0.1, "0.01<0.1", 
                                    ifelse(change > 0.1 & change < 0.2, "0.1<0.2", 
                                           ifelse(change > 0.2 & change < 0.3, "0.2<0.3",">0.3"))))) %>% 
        summarise(dplyr::count(.,skew)) -> df3

df4 <- data.frame(category = factor(df3[,1], levels = c("<0.01", "0.01<0.1", "0.1<0.2", "0.2<0.3", ">0.3")),
                  count = df3[,2])

p4 <- ggplot(data = df4, aes(x = category))+
        geom_bar(aes(y = count), stat = "identity")+
        geom_text(aes(label = count, y = count), vjust = -0.5)+
        theme_classic(base_size = 20)+
        labs(x = "B-value skew")+
        theme(axis.title.y = element_blank())+
        scale_y_continuous(trans = "log2", expand  = c(0, 1), n.breaks = 9)


plot_grid(p1,p3,p2,p4, rel_widths = c(1.5,1), labels = c("A", "B", "C", "D"), label_size = 20)



# count baseline DMPs for islands within promoters

DMPs_BM_vs_BH %>% 
        as.data.frame() %>% 
        merge(., anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island") %>% 
        nrow()


# count post dmps for islnds within promoters

DMPs_PM_vs_PH %>% 
        as.data.frame() %>% 
        merge(., anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" ) %>% 
        nrow()


###########################################################################################

### DMP venn diagrams within time, across cell populations

###########################################################################################

# extract probe names of DMPs that are negatively and positively skjewed at each timepoint


base_neg <- DMPs_BM_vs_BH %>% 
        as.data.frame() %>% 
        mutate(dir = ifelse(delta_M <0, "hypo", "hyper")) %>% 
        filter(dir == "hypo")

base_pos <- DMPs_BM_vs_BH %>% 
        as.data.frame() %>% 
        mutate(dir = ifelse(delta_M <0, "hypo", "hyper")) %>% 
        filter(dir == "hyper")

post_neg <- DMPs_PM_vs_PH %>% 
        as.data.frame() %>% 
        mutate(dir = ifelse(delta_M <0, "hypo", "hyper")) %>% 
        filter(dir == "hypo")

post_pos <- DMPs_PM_vs_PH %>% 
        as.data.frame() %>% 
        mutate(dir = ifelse(delta_M <0, "hypo", "hyper")) %>% 
        filter(dir == "hyper")

DMPs <- list(
        base_neg = base_neg$cpg,
        Baseline = base_pos$cpg,
        Post = post_neg$cpg,
        post_pos = post_pos$cpg
)

venn <- ggvenn(DMPs, set_name_size = 10, stroke_size = 1, 
       fill_color = c("#453781FF", "#DCE319FF","#453781FF", "#DCE319FF"),
       text_size = 10,
       stroke_alpha = 0.8,
       set_name_color = c("White", "Black","Black","White"))

plot_grid(venn, labels = "C", label_size = 27, hjust = -4.5)

# get the cpgs of the differentially overlapping genes

merge(base_pos, post_neg, by = "cpg") %>% 
        merge(.,anno, by = "cpg") %>% 
        dplyr::select(cpg, UCSC_RefGene_Name)

merge(base_neg, post_pos, by = "cpg")%>% 
        merge(.,anno, by = "cpg") %>% 
        dplyr::select(cpg, UCSC_RefGene_Name)

# find position of the cpgs in homogenate and myonuclei time DMPs

DMPs_PH_vs_BH %>% 
        rownames_to_column(var = "number") %>% 
        dplyr::select(number, cpg) %>% 
        merge(merge(base_pos, post_neg, by = "cpg"),.,by = "cpg") %>% 
        select(cpg, number) %>% 
        dplyr::arrange(number) -> x

DMPs_PM_vs_BM %>% 
        rownames_to_column(var = "number") %>% 
        dplyr::select(number, cpg) %>% 
        merge(merge(base_pos, post_neg, by = "cpg"),.,by = "cpg") %>% 
        select(cpg, number) %>% 
        dplyr::arrange(number) -> y

merge(x,y, by = "cpg") %>% 
        pull(cpg) -> z

DMPs_PH_vs_BH %>% 
        filter(cpg %in% z)

DMPs_PM_vs_BM %>% 
        filter(cpg %in% z)

DMPs_PH_vs_BH %>% 
        rownames_to_column(var = "number") %>% 
        dplyr::select(number, cpg) %>% 
        merge(merge(base_neg, post_pos, by = "cpg"),.,by = "cpg") %>% 
        select(cpg, number) %>% 
        dplyr::arrange(number) -> x

DMPs_PM_vs_BM %>% 
        rownames_to_column(var = "number") %>% 
        dplyr::select(number, cpg) %>% 
        merge(merge(base_neg, post_pos, by = "cpg"),.,by = "cpg") %>% 
        select(cpg, number) %>% 
        dplyr::arrange(number) -> y

merge(x,y, by = "cpg") %>% 
        pull(cpg) -> z

DMPs_PH_vs_BH %>% 
        filter(cpg %in% z)

DMPs_PM_vs_BM %>% 
        filter(cpg %in% z)

# create supplementary csv of overlapping 20 DMPs in venn diagram figure 6A

DMPs_BM_vs_BH %>% 
        mutate()

merge(base_neg, post_pos, by = "cpg") %>% pull(cpg)
merge(base_pos, post_neg, by = "cpg") %>% pull(cpg)


########################################################################################

### Heatmaps DMPs in islands at baseline and at post

########################################################################################

# isolate beta values of DMPs in islands

# get baseline DMP names

baseline_DMPs <- DMPs_BM_vs_BH %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Relation_to_Island == "Island") %>% 
        pull(cpg)

post_DMPs <- DMPs_PM_vs_PH %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Relation_to_Island == "Island") %>% 
        pull(cpg)

DMPs <- list(baseline_isl = baseline_DMPs,
             post_isl = post_DMPs)


ggvenn(DMPs, set_name_size = 10, stroke_size = 1, 
      # fill_color = c("#453781FF", "#DCE319FF","#453781FF", "#DCE319FF"),
       text_size = 10,
       stroke_alpha = 0.8,
      # set_name_color = c("White", "Black","Black","White"),
      )


# take all DMPs from baseline and post and plot in heatmaps

heatmap_beta <- beta %>% 
        as.data.frame() %>% 
        dplyr::select("1BH","2BH","4BH","5BH","6BH","7BH","8BH","12BH",
                      "1BM","2BM","4BM","5BM","6BM","7BM","8BM","12BM",
                      "1PH","2PH","4PH","5PH","6PH","7PH","8PH","12PH",
                      "1PM","2PM","4PM","5PM","6PM","7PM","8PM","12PM") %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% baseline_DMPs | cpg %in% post_DMPs) 

# plot baseline heatmap         


baseline_heatmap <- pheatmap(t(heatmap_beta[,2:17]), annotation_row = dfh_1, cluster_rows = FALSE, gaps_row = 8, 
                          show_colnames = FALSE, annotation_names_row = FALSE, 
                          color = viridis(n = 100), scale = "none")

# baseline heatmap looks good

post_heatmap <- pheatmap(t(heatmap_beta[,18:33]), annotation_row = dfh_2, cluster_rows = FALSE, gaps_row = 8, 
                             show_colnames = FALSE, annotation_names_row = FALSE, 
                             color = viridis(n = 100), scale = "none")

#need to flip dendrogram in post heatmap



col_dend <- post_heatmap[[2]]

col_dend$order <- rev(col_dend$order)

post_heatmap <- pheatmap(t(heatmap_beta[,18:33]), annotation_row = dfh_2, cluster_rows = FALSE, gaps_row = 8, 
                      show_colnames = FALSE, annotation_names_row = FALSE, 
                      color = viridis(n = 100), scale = "none", cluster_cols = as.hclust(col_dend))

plot_heat <- plot_grid(baseline_heatmap[[4]], post_heatmap[[4]], ncol=1, labels = c("A", "B"))   

ggsave2(plot_heat, filename = "heatmap_homo_vs_myo.pdf",units = "cm", width = 19, height = 21, bg = "white")


# remake heatmap with the baseline dendrogram for both plots


col_dend <- baseline_heatmap[[2]]

post_heatmap2 <- pheatmap(t(heatmap_beta[,18:33]), annotation_row = dfh_2, cluster_rows = FALSE, gaps_row = 8, 
                         show_colnames = FALSE, annotation_names_row = FALSE, 
                         color = viridis(n = 100), scale = "none", cluster_cols = as.hclust(col_dend))

plot_heat <- plot_grid(baseline_heatmap[[4]], post_heatmap2[[4]], ncol=1, labels = c("A", "B"))


######################################################################################

### k means plots

#######################################################################################

som_data <-  B2M(beta)%>% 
        as.data.frame() %>% 
        mutate(avg_BH = round((`1BH`+`2BH`+`4BH`+`5BH`+`6BH`+`7BH`+`8BH`+`12BH`)/8,2),
               avg_BM = round((`1BM`+`2BM`+`4BM`+`5BM`+`6BM`+`7BM`+`8BM`+`12BM`)/8,2),
               avg_PH = round((`1PH`+`2PH`+`4PH`+`5PH`+`6PH`+`7PH`+`8PH`+`12PH`)/8,2),
               avg_PM = round((`1PM`+`2PM`+`4PM`+`5PM`+`6PM`+`7PM`+`8PM`+`12PM`)/8,2)) %>% 
        dplyr::select(avg_BH, avg_BM, avg_PH, avg_PM) %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::filter(cpg %in% c(DMPs_BM_vs_BH$cpg,DMPs_PM_vs_PH$cpg))


rownames(som_data) <- som_data[,1]
som_data <- som_data[,-1]

### run kmeans clustering on som_data



set.seed(1)


xyss <- vector()

for (i in 1:10) {
        xyss[i] <- sum(kmeans(som_data,i)$withinss)
}
plot(1:10, xyss, type = "b", main = "clusters of M-value profiles", xlab = "number of clusters", ylab = "XYSS")

# elbow at ~4 clusters

# make cluster with identified elbow
set.seed(1)
kmeans <- kmeans(som_data, 4, iter.max = 300, nstart = 10)

# add cluster number to som data and replot

k_means <-as.data.frame(kmeans$cluster)

som_data[,"kmeans"] <- k_means[,1]


# remake som plot with individual y axis scaling

# count cluster DMPs

som_data %>% 
        filter(kmeans == "4") %>% nrow()



k1 <- som_data %>% 
        pivot_longer(names_to = "condition", values_to = "mean_M", cols = 1:4) %>% 
        mutate(condition = factor(condition, levels = c("avg_BH", "avg_BM", "avg_PH", "avg_PM")),
               kmeans = factor(kmeans))%>% 
        dplyr::group_by(kmeans,condition) %>% 
        dplyr::summarize(m = mean(mean_M),
                         s = sd(mean_M)) %>% 
        mutate(sample = substr_right(as.character(condition), 1),
               timepoint = substr_right(as.character(condition), 2)) %>% 
        mutate(sample = ifelse(sample == "M", "Myonuclei", "Homogenate"),
               timepoint = substr(timepoint, 1,1)) %>% 
        mutate(timepoint = ifelse(timepoint == "B", "Baseline", "Post")) %>% 
        mutate(timepoint = factor(timepoint, levels = c("Baseline", "Post")),
               sample = factor(sample, levels = c("Myonuclei","Homogenate"))) %>% 
        filter(kmeans == "2") %>% 
        ggplot(aes(x = timepoint, y = m, color = sample))+
        geom_line(aes(group = sample), size = 1.5, position = position_dodge(width = 0.2))+
        geom_errorbar(aes(ymin = m-s,
                          ymax = m+s, x = timepoint, group = sample),inherit.aes = FALSE, width = 0.2, size = 1, position = position_dodge(width = 0.2), alpha = 0.5)+
        labs(y = "M-value",
             title = "Cluster 1: N = 90250")+
        theme_bw(base_size = 20)+
        theme(axis.text.x = element_text(size = 15),
              legend.position = "none",
              axis.title.x = element_blank(),
              panel.grid.minor = element_blank())+
        scale_y_continuous(n.breaks = 4)+
        geom_text(aes(label = round(m,2), y = m),position = position_dodge(width = 1), vjust = 1.5,fontface = "bold", size = 8)


k2 <- som_data %>% 
        pivot_longer(names_to = "condition", values_to = "mean_M", cols = 1:4) %>% 
        mutate(condition = factor(condition, levels = c("avg_BH", "avg_BM", "avg_PH", "avg_PM")),
               kmeans = factor(kmeans))%>% 
        dplyr::group_by(kmeans,condition) %>% 
        dplyr::summarize(m = mean(mean_M),
                         s = sd(mean_M)) %>% 
        mutate(sample = substr_right(as.character(condition), 1),
               timepoint = substr_right(as.character(condition), 2)) %>% 
        mutate(sample = ifelse(sample == "M", "Myonuclei", "Homogenate"),
               timepoint = substr(timepoint, 1,1)) %>% 
        mutate(timepoint = ifelse(timepoint == "B", "Baseline", "Post")) %>% 
        mutate(timepoint = factor(timepoint, levels = c("Baseline", "Post")),
               sample = factor(sample, levels = c("Myonuclei","Homogenate"))) %>% 
        filter(kmeans == "4") %>% 
        ggplot(aes(x = timepoint, y = m, color = sample))+
        geom_line(aes(group = sample), size = 1.5, position = position_dodge(width = 0.2))+
        geom_errorbar(aes(ymin = m-s,
                          ymax = m+s, x = timepoint, group = sample),inherit.aes = FALSE, width = 0.2, size = 1, position = position_dodge(width = 0.2), alpha = 0.5)+
        labs(y = "M-value",
             title = "Cluster 2: N = 83018")+
        theme_bw(base_size = 20)+
        theme(axis.text.x = element_text(size = 15),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.minor = element_blank())+
        geom_text(aes(label = round(m,2), y = m),position = position_dodge(width = 1), vjust = 1.5,fontface = "bold", size = 8)


k3 <- som_data %>% 
        pivot_longer(names_to = "condition", values_to = "mean_M", cols = 1:4) %>% 
        mutate(condition = factor(condition, levels = c("avg_BH", "avg_BM", "avg_PH", "avg_PM")),
               kmeans = factor(kmeans))%>% 
        dplyr::group_by(kmeans,condition) %>% 
        dplyr::summarize(m = mean(mean_M),
                         s = sd(mean_M)) %>% 
        mutate(sample = substr_right(as.character(condition), 1),
               timepoint = substr_right(as.character(condition), 2)) %>% 
        mutate(sample = ifelse(sample == "M", "Myonuclei", "Homogenate"),
               timepoint = substr(timepoint, 1,1)) %>% 
        mutate(timepoint = ifelse(timepoint == "B", "Baseline", "Post")) %>% 
        mutate(timepoint = factor(timepoint, levels = c("Baseline", "Post")),
               sample = factor(sample, levels = c("Myonuclei","Homogenate"))) %>% 
        filter(kmeans == "3") %>% 
        ggplot(aes(x = timepoint, y = m, color = sample))+
        geom_line(aes(group = sample), size = 1.5, position = position_dodge(width = 0.2))+
        geom_errorbar(aes(ymin = m-s,
                          ymax = m+s, x = timepoint, group = sample),inherit.aes = FALSE, width = 0.2, size = 1, position = position_dodge(width = 0.2), alpha = 0.5)+
        labs(y = "M-value",
             title = "Cluster 3: N = 67620")+
        theme_bw(base_size = 20)+
        theme(axis.text.x = element_text(size = 15),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.minor = element_blank())+
        geom_text(aes(label = round(m,2), y = m),position = position_dodge(width = 1), vjust = 1.5,fontface = "bold", size = 8)


k4 <- som_data %>% 
        pivot_longer(names_to = "condition", values_to = "mean_M", cols = 1:4) %>% 
        mutate(condition = factor(condition, levels = c("avg_BH", "avg_BM", "avg_PH", "avg_PM")),
               kmeans = factor(kmeans))%>% 
        dplyr::group_by(kmeans,condition) %>% 
        dplyr::summarize(m = mean(mean_M),
                         s = sd(mean_M)) %>% 
        mutate(sample = substr_right(as.character(condition), 1),
               timepoint = substr_right(as.character(condition), 2)) %>% 
        mutate(sample = ifelse(sample == "M", "Myonuclei", "Homogenate"),
               timepoint = substr(timepoint, 1,1)) %>% 
        mutate(timepoint = ifelse(timepoint == "B", "Baseline", "Post")) %>% 
        mutate(timepoint = factor(timepoint, levels = c("Baseline", "Post")),
               sample = factor(sample, levels = c("Myonuclei","Homogenate"))) %>% 
        filter(kmeans == "1") %>% 
        ggplot(aes(x = timepoint, y = m, color = sample))+
        geom_line(aes(group = sample), size = 1.5, position = position_dodge(width = 0.2))+
        geom_errorbar(aes(ymin = m-s,
                          ymax = m+s, x = timepoint, group = sample),inherit.aes = FALSE, width = 0.2, size = 1, position = position_dodge(width = 0.2), alpha = 0.5)+
        labs(y = "M-value",
             title = "Cluster 4: N = 35842")+
        theme_bw(base_size = 20)+
        theme(axis.text.x = element_text(size = 15),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.minor = element_blank())+
        geom_text(aes(label = round(m,2), y = m),position = position_dodge(width = 1), vjust = 1.5,fontface = "bold", size = 8)





# plot kmeans
kmeans_plot <- ggarrange(k1,k2,k3,k4,common.legend = TRUE, legend = "right", nrow = 1,widths = c(1,0.9,0.9,0.9), labels = c("D", "","",""), font.label = list(size = 27))

# plot venn and kmeans together
plot <- plot_grid(venn, kmeans_plot, ncol = 1, rel_heights = c(2,1), labels = c("C",""), label_size = 27)

ggsave2(plot, filename = "venn_kmeans.pdf",units = "cm", width = 19, height = 21, bg = "white")
# plot with heatmap (doesnt wirk)
plot_grid(plot_heat, plot)

#############################################################################

### Kmeans plot of 901 overlapping DMPs in myo+int and myo

##########################################################################

# isolate beta values of 901 cpgs

# identify cpgs

common_dmps <- merge(DMPs_PH_vs_BH, DMPs_PM_vs_BM, by = "cpg") %>% 
        pull(cpg) 

som_data <-B2M(beta)%>% 
        as.data.frame() %>% 
        mutate(avg_BH = round((`1BH`+`2BH`+`4BH`+`5BH`+`6BH`+`7BH`+`8BH`+`12BH`)/8,2),
               avg_BM = round((`1BM`+`2BM`+`4BM`+`5BM`+`6BM`+`7BM`+`8BM`+`12BM`)/8,2),
               avg_PH = round((`1PH`+`2PH`+`4PH`+`5PH`+`6PH`+`7PH`+`8PH`+`12PH`)/8,2),
               avg_PM = round((`1PM`+`2PM`+`4PM`+`5PM`+`6PM`+`7PM`+`8PM`+`12PM`)/8,2)) %>% 
        dplyr::select(avg_BH, avg_BM, avg_PH, avg_PM) %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::filter(cpg %in% common_dmps) %>% 
        mutate(change_homo = avg_PH-avg_BH, 
               change_myo = avg_PM-avg_BM) %>% 
        pivot_longer(names_to = "sample", values_to = "average", cols = 2:7) %>% 
        filter(sample %in% c("change_homo","change_myo")) %>% 
        pivot_wider(names_from = sample, values_from = average) %>% 
        as.data.frame()

 


rownames(som_data) <- som_data[,1]
som_data <- som_data[,-1]


### run kmeans clustering on som_data



set.seed(1)


xyss <- vector()

for (i in 1:10) {
        xyss[i] <- sum(kmeans(som_data,i)$withinss)
}
plot(1:10, xyss, type = "b", main = "clusters of M-value profiles", xlab = "number of clusters", ylab = "XYSS")

# elbow at ~4 clusters

# make cluster with identified elbow
set.seed(1)
kmeans <- kmeans(som_data, 4, iter.max = 300, nstart = 10)

# add cluster number to som data and replot

k_means <-as.data.frame(kmeans$cluster)

som_data[,"kmeans"] <- k_means[,1]


# remake som plot with individual y axis scaling

# count cluster DMPs

som_data %>% 
        filter(kmeans == "4") %>% nrow()


# quick plot

som_data %>% 
        pivot_longer(names_to = "condition", values_to = "mean_M", cols = 1:4) %>% 
        mutate(condition = factor(condition, levels = c("avg_BH", "avg_BM", "avg_PH", "avg_PM")),
               cluster = factor(kmeans)) %>% 
        dplyr::group_by(cluster,condition) %>% 
        dplyr::summarize(m = mean(mean_M),
                         s = sd(mean_M)) %>% 
        ggplot(aes(x = condition, y = m))+
        geom_line(aes(group = cluster))+
        geom_errorbar(aes(ymin = m-s,
                          ymax = m+s), width = 0.2)+
        facet_grid(~cluster)



#################################################################################

### total DMPs

###################################################################################

# count all DMPs and plot

# count DMPs and plot

hypo_col = "#453781FF"
hyper_col = "#DCE319FF"


DMPs_PH_vs_BH %>%  
        dplyr::select(delta_M) %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0)) %>% 
        mutate(comp = "Homo") -> x


DMPs_PM_vs_BM %>%  
        dplyr::select(delta_M) %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0)) %>% 
        mutate(comp = "Myo") %>% 
        base::rbind(x,.)->x





p1 <- x %>% 
        mutate(percent_hyper = hyper/(hypo+hyper)*100,
               percent_hypo = hypo/(hypo+hyper)*100) %>% 
        ggplot(aes(x = comp))+
        geom_bar(aes(y = -hypo), stat = "identity", position = "dodge", width = 0.8, fill = hypo_col)+
        geom_bar(aes(y = hyper), stat = "identity", position = "dodge", width = 0.8, fill = hyper_col)+
        geom_label(aes(x = comp, y = hyper, label = paste(hyper, "/", round(percent_hyper, 1), "%")), alpha = 0.8)+
        geom_label(aes(x = comp, y = -hypo, label = paste(hypo, "/", round(percent_hypo, 1), "%")), alpha = 0.8)+
        theme_classic(base_size = 20) +
        scale_x_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())+
        labs(y = "DMPs after 7 weeks RT (un adj.p < 0.05" )


# plot DMPs split into regulatory features

hDMP <- DMPs_PH_vs_BH%>% 
        merge(.,anno, by = "cpg") %>% 
        group_by(Relation_to_Island) %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) 

mDMP <- DMPs_PM_vs_BM%>% 
        merge(.,anno, by = "cpg") %>% 
        group_by(Relation_to_Island) %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100)


# plot homogenate post vs. homogenate baseline DMPs

p2 <- ggplot(data = hDMP, aes(x = Relation_to_Island))+
        geom_bar(aes(y = hyper), stat = "identity", width = 0.9, fill = hyper_col)+
        geom_bar(aes(y = -hypo), stat = "identity", width = 0.9, fill = hypo_col)+
        theme_classic(base_size = 20)+
        labs(y = "MYO + INT DMPs after 7 weeks RT")+
        geom_label(aes(x = Relation_to_Island, y = -6000, label = total))+
        geom_label(aes(x = Relation_to_Island, y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 0), "%")), vjust = c(0.3,-0.4,0.3,1.5,-0.4,0.3), alpha = 0.6)+
        geom_label(aes(x = Relation_to_Island, y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 0), "%")), vjust = c(0.7,1.3,0.7,-0.6,1.3,0.7), alpha = 0.6)+
        theme(axis.title.x = element_blank())

# plot myonuclei post vs. myonuclei baseline DMPs

p3 <- ggplot(data = mDMP, aes(x = Relation_to_Island))+
        geom_bar(aes(y = hyper), stat = "identity", fill = hyper_col, width = 0.9)+
        geom_bar(aes(y = -hypo), stat = "identity", fill = hypo_col, width = 0.9)+
        theme_classic(base_size = 20)+
        labs(y = "MYO DMPs after 7 weeks RT")+
        geom_label(aes(x = Relation_to_Island, y = -11000, label = total))+
        geom_label(aes(x = Relation_to_Island, y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 0), "%")), vjust = c(0.3,-0.4,0.3,1.5,-0.4,0.3), alpha = 0.6)+
        geom_label(aes(x = Relation_to_Island, y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 0), "%")), vjust = c(0.7,1.3,0.7,-0.6,1.3,0.7), alpha = 0.6)+
        theme(axis.title.x = element_blank())+
        scale_fill_identity(name = 'the fill', guide = "legend",labels = c("hypo", "hyper"),aes(y = 0, x = 7))



plot_grid(p1,p2,nrow = 1, labels = c("A","B"))


# isolate homogenate and myonuclei DMPs within promoter assosiated islands and plot together


hDMP <- DMPs_PH_vs_BH %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island") %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) %>% 
        mutate(condition = "Homo")

mDMP <- DMPs_PM_vs_BM%>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island") %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) %>% 
        mutate(condition = "Myo")

p4 <- rbind(hDMP, mDMP) %>% 
        ggplot(aes(x = condition))+
        geom_bar(aes(y = hyper), stat = "identity", fill = hyper_col, width = 0.8)+
        geom_bar(aes(y = -hypo), stat = "identity", fill = hypo_col, width = 0.8)+
        geom_label(aes(y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), alpha = 0.8)+
        geom_label(aes(y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), alpha = 0.8)+
        theme_classic(base_size = 20) +
        scale_x_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())+
        labs(y = "DMPs in Islands and Promoters")

DMPs <- list(
        Homogenate_hypo = DMPs_PH_vs_BH %>% 
                filter(delta_M < 0) %>% 
                pull(cpg),
        "MYO+INT" = DMPs_PH_vs_BH %>% 
                filter(delta_M > 0) %>% 
                pull(cpg),
        "MYO" = DMPs_PM_vs_BM %>% 
                filter(delta_M < 0) %>% 
                pull(cpg),
        Myonuclei_hyper = DMPs_PM_vs_BM %>% 
                filter(delta_M > 0) %>% 
                pull(cpg)
)

p5 <- ggvenn(DMPs, set_name_size = 10, stroke_size = 1, set_name_color = c("White", "Black", "Black","White"),
       fill_color = c("#453781FF", "#DCE319FF","#453781FF", "#DCE319FF"),text_size = 8,stroke_alpha = 0.8)

p6 <- plot_grid(p5, labels = "B", label_x = 0.1, label_size = 25)


# figure 8 in manuscript

part1 <- plot_grid(p1,p5, nrow = 1, labels = c("A", "B"), rel_widths = c(0.65,1), label_size = 25)

part2 <- plot_grid(p3,p2,p4, nrow = 1, labels = c("C","D", "E"), rel_widths = c(1,1,0.5), label_size = 25)

DMP_RT_plot <- plot_grid(part1, part2, ncol = 1)

ggsave2(DMP_RT_plot, filename = "DMP_RT_plot.pdf",units = "cm", width = 19, height = 21, bg = "white")


# count promoter DMPs at baseline

DMPs_BM_vs_BH %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated") %>%  nrow()

####################################################################################################

### volcano plots - most significant probes 

####################################################################################################

# DMPs in homogenate after RT  

# create functions to transform axis

reverselog_trans <- function(base = exp(1)) {
        trans <- function(x) -log(x, base)
        inv <- function(x) base^(-x)
        trans_new(paste0("reverselog-", format(base)), trans, inv, 
                  log_breaks(base = base), 
                  domain = c(1e-100, Inf))
}

# calculate log2 fold change on beta values og significant DMPs

beta %>% 
        as.data.frame() %>% 
        dplyr::select("1BH","2BH","4BH","5BH","6BH","7BH","8BH","12BH",
                      "1BM","2BM","4BM","5BM","6BM","7BM","8BM","12BM",
                      "1PH","2PH","4PH","5PH","6PH","7PH","8PH","12PH",
                      "1PM","2PM","4PM","5PM","6PM","7PM","8PM","12PM") %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% DMPs_PH_vs_BH$cpg) %>%
        mutate(logFC = ((((`1PH`+`2PH`+`4PH`+`5PH`+`6PH`+`7PH`+`8PH`+`12PH`)/8)/
                                    ((`1BH`+`2BH`+`4BH`+`5BH`+`6BH`+`7BH`+`8BH`+`12BH`)/8))-1)) %>% 
        dplyr::select(cpg, logFC) -> logFC_homo

set.seed(1)
DMPs_PH_vs_BH %>% 
        merge(., anno, by = "cpg") %>%
        mutate(Gene = ifelse(delta_M < -0.8 | delta_M > 0.8 | p.value < 0.0001, as.character(UCSC_RefGene_Name), ""),
               Gene = ifelse(Gene == "NA", "", Gene)) %>%
        ggplot(aes(x = delta_M, y = p.value))+
        geom_point(aes(color = abs(delta_M), shape = Relation_to_Island), size = 3)+
        scale_color_viridis()+
        scale_y_continuous(trans = reverselog_trans(10),
                           labels = label_comma())+
        scale_x_continuous(n.breaks = 6)+
        geom_text_repel(aes(label = Gene), force = 2, max.overlaps = Inf, box.padding = 0.6, point.size = 2)+
        theme_classic(base_size = 20)+
        theme(panel.grid.major.y = element_line(color = "red", 
                                                size = 0.5,
                                                linetype = 2),
              axis.title.x = element_blank())+
        labs(color = "delta M", shape = "Relation to CpG Island",
             y = "un-adj. P value")


# get the dataframe

DMPs_PH_vs_BH %>% 
        merge(., anno, by = "cpg") %>%
        mutate(Gene = ifelse(delta_M < -0.8 | delta_M > 0.8 | p.value < 0.0001, as.character(UCSC_RefGene_Name), ""),
               Gene = ifelse(Gene == "NA", "", Gene)) %>%
        filter(Gene != "") %>% 
        dplyr::select(!7) %>% 
        arrange(desc(-p.value)) %>% pull(cpg) -> y
        write.csv(file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/Homogenate_DMPs_after_RT_volcano-plot.csv")

# DMPs in myonuclei after RT

set.seed(1)
DMPs_PM_vs_BM %>% 
        merge(., anno, by = "cpg") %>%
        mutate(Gene = ifelse(delta_M < -0.8 | delta_M > 0.8 | p.value < 0.0001, as.character(UCSC_RefGene_Name), ""),
               Gene = ifelse(Gene == "NA", "", Gene)) %>%
        ggplot(aes(x = delta_M, y = p.value))+
        geom_point(aes(color = abs(delta_M), shape = Relation_to_Island), size = 3)+
        scale_color_viridis()+
        scale_y_continuous(trans = reverselog_trans(10),
                           labels = label_comma())+
        scale_x_continuous(n.breaks = 6)+
        geom_text_repel(aes(label = Gene), force = 2, max.overlaps = Inf, box.padding = 0.6, point.size = 2)+
        theme_classic(base_size = 20)+
        theme(panel.grid.major.y = element_line(color = "red", 
                                                size = 0.5,
                                                linetype = 2),
              axis.title.x = element_blank())+
        labs(color = "delta M", shape = "Relation to CpG Island",
             y = "un-adj. P value")

DMPs_PM_vs_BM %>% 
        merge(., anno, by = "cpg") %>%
        mutate(Gene = ifelse(delta_M < -0.8 | delta_M > 0.8 | p.value < 0.0001, as.character(UCSC_RefGene_Name), ""),
               Gene = ifelse(Gene == "NA", "", Gene)) %>%
        filter(Gene != "") %>% 
        dplyr::select(!7) %>% 
        arrange(desc(-p.value)) %>% pull(cpg) -> x
        write.csv(file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/Myonuclear_DMPs_after_RT_volcano-plot.csv")



# get probe/gene data on homogenate and myonuclei DMPs in same dataframe
        
        myo_int <- readRDS("DMPs_PH_vs_BH.RDATA")
        myo <- readRDS("DMPs_PM_vs_BM.RDATA")
        
myo %>% 
        filter(cpg %in% x) %>% 
        mutate(sample = "MYO") %>% 
        rbind(., myo_int %>% filter(cpg %in% x) %>%  mutate(sample = "MYO_INT"))



# MYO DMPS: hypo MYO - Hyper MYO+INT
myo %>% 
        filter(cpg %in% x) %>% 
        mutate(sample = "MYO") %>% 
        merge(., myo_int %>% filter(cpg %in% x) %>%  mutate(sample = "MYO+INT"), by = "cpg") %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(delta_M.x <0 & delta_M.y > 0) %>% pull(UCSC_RefGene_Name)

# MYO DMPS: hyper MYO - Hypo MYO+INT

myo %>% 
        filter(cpg %in% x) %>% 
        mutate(sample = "MYO") %>% 
        merge(., myo_int %>% filter(cpg %in% x) %>%  mutate(sample = "MYO+INT"), by = "cpg") %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(delta_M.x >0 & delta_M.y < 0) %>% pull(UCSC_RefGene_Name)


anno %>% 
        filter(cpg %in% x | cpg %in% y)

myo_int%>% 
        filter(cpg %in% y) %>% 
        mutate(sample = "MYO+INT") %>% 
        merge(.,myo, by = "cpg") %>% 
        mutate("myo_vs_myo_int" = delta_M.y-delta_M.x) %>% 
        filter(delta_M.x < 0 & delta_M.y > 0) %>% 
        merge(., anno, by = "cpg") %>% 
        arrange(desc(myo_vs_myo_int))

myo_int %>% 
        filter(cpg %in% y)

merge(DMPs_PH_vs_BH, DMPs_PM_vs_BM, by = "cpg") %>% nrow()
nrow(DMPs_PM_vs_BM)

#####################################################################################################

### horizontal mean +- confidence interval for the top 20 highest delta M

##################################################################################################

# load and merge homogenate
"#453781FF", "#DCE319FF"

# homogenate

# plot the 10 highest and lowest delta_M values

DMPs_PH_vs_BH %>% 
        merge(.,anno,by = "cpg") %>% 
        arrange(delta_M) %>% 
        head(10)-> x
        
        

p1 <- DMPs_PH_vs_BH %>% 
        merge(.,anno,by = "cpg") %>% 
        arrange(delta_M) %>% 
        tail(10) %>% 
        rbind(x,.) %>% 
        mutate(name = fct_reorder(paste(cpg,";", UCSC_RefGene_Name), desc(delta_M)))%>%
        ggplot(aes(x = delta_M, y = name)) +
        geom_rect(aes(xmin = -Inf, xmax = 0, 
                      ymin = -Inf, ymax = Inf), fill = "#453781FF", alpha = 0.05)+
        geom_rect(aes(xmin = 0, xmax = Inf, 
                      ymin = -Inf, ymax = Inf), fill = "#DCE319FF", alpha = 0.05)+
        geom_point(size = 3)+
        geom_errorbar(aes(xmin = conf.int_0.025, xmax = conf.int_0.975), width = 0.2, size = 1)+
        geom_vline(xintercept = 0, color = "Red", size = 1.1) + 
        theme_classic()+
        theme(axis.title.y = element_blank())
        
        


" plot the lowest 20 p.values"

p2 <- DMPs_PH_vs_BH %>% 
        merge(.,anno,by = "cpg") %>% 
        arrange(p.value) %>% 
        head(20) %>% 
        mutate(name = fct_reorder(paste(cpg,";", UCSC_RefGene_Name), desc(delta_M)))%>%
        ggplot(aes(x = delta_M, y = name)) +
        geom_rect(aes(xmin = -Inf, xmax = 0, 
                      ymin = -Inf, ymax = Inf), fill = "#453781FF", alpha = 0.05)+
        geom_rect(aes(xmin = 0, xmax = Inf, 
                      ymin = -Inf, ymax = Inf), fill = "#DCE319FF", alpha = 0.05)+
        geom_point(size = 3)+
        geom_errorbar(aes(xmin = conf.int_0.025, xmax = conf.int_0.975), width = 0.2, size = 1)+
        geom_vline(xintercept = 0, color = "Red", size = 1.1) + 
        theme_classic()+
        theme(axis.title.y = element_blank())


# myonuclei

# plot the 10 highest and lowest delta_M values

DMPs_PM_vs_BM %>% 
        merge(.,anno,by = "cpg") %>% 
        arrange(delta_M) %>% 
        head(10)-> x



p3 <- DMPs_PM_vs_BM %>% 
        merge(.,anno,by = "cpg") %>% 
        arrange(delta_M) %>% 
        tail(10) %>% 
        rbind(x,.) %>% 
        mutate(name = fct_reorder(paste(cpg,";", UCSC_RefGene_Name), desc(delta_M)))%>%
        ggplot(aes(x = delta_M, y = name)) +
        geom_rect(aes(xmin = -Inf, xmax = 0, 
                      ymin = -Inf, ymax = Inf), fill = "#453781FF", alpha = 0.05)+
        geom_rect(aes(xmin = 0, xmax = Inf, 
                      ymin = -Inf, ymax = Inf), fill = "#DCE319FF", alpha = 0.05)+
        geom_point(size = 3)+
        geom_errorbar(aes(xmin = conf.int_0.025, xmax = conf.int_0.975), width = 0.2, size = 1)+
        geom_vline(xintercept = 0, color = "Red", size = 1.1) + 
        theme_classic()+
        theme(axis.title.y = element_blank())




" plot the lowest 20 p.values"

p4 <- DMPs_PM_vs_BM %>% 
        merge(.,anno,by = "cpg") %>% 
        arrange(p.value) %>% 
        head(20) %>% 
        mutate(name = fct_reorder(paste(cpg,";", UCSC_RefGene_Name), desc(delta_M)))%>%
        ggplot(aes(x = delta_M, y = name)) +
        geom_rect(aes(xmin = -Inf, xmax = 0, 
                      ymin = -Inf, ymax = Inf), fill = "#453781FF", alpha = 0.05)+
        geom_rect(aes(xmin = 0, xmax = Inf, 
                      ymin = -Inf, ymax = Inf), fill = "#DCE319FF", alpha = 0.05)+
        geom_point(size = 3)+
        geom_errorbar(aes(xmin = conf.int_0.025, xmax = conf.int_0.975), width = 0.2, size = 1)+
        geom_vline(xintercept = 0, color = "Red", size = 1.1) + 
        theme_classic()+
        theme(axis.title.y = element_blank())

plot_grid(p1,p2,p3,p4, nrow = 1, labels = c("A","B","C","D"))



















##############################################################################################

### DMRs

########################################################################################

# I will use the DMRcate method from new oshlack workflow



# create design matrix

dfh %>% 
        filter(condition == "BH" | condition == "PH") -> dfh3


m_vals_homo <- B2M(beta) %>% 
        as.data.frame() %>% 
        dplyr::select(rownames(dfh3)) 


# Sample names
sample_names <- c("1BH", "2BH", "4BH", "5BH", "6BH", "7BH", "8BH", "12BH",
                  "1PH", "2PH", "4PH", "5PH", "6PH", "7PH", "8PH", "12PH")

# Create a vector of group labels
group_labels <- ifelse(grepl("B", sample_names), "Baseline", "Post")

# Create a factor variable with the group labels
group_factor <- factor(group_labels, levels = c("Baseline", "Post"))

# Create the design matrix
design_matrix <- model.matrix(~ group_factor)


# Print the design matrix
design_matrix

annotated_data_homo <- cpg.annotate(datatype = "array", 
                               object = as.matrix(m_vals_homo), 
                               what = "M", 
                               arraytype = "EPIC", 
                               analysis.type = "differential", 
                               design_matrix,
                               fdr = "none", 
                               pval = 0.05,
                               coef = "group_factorPost")

dmrcate_results <- dmrcate(annotated_data_homo, 
                           lambda = 1000, 
                           C = 2, 
                           pcutoff = "fdr")


results.ranges <- extractRanges(dmrcate_results)
results.ranges %>% 
        as.data.frame() %>% 
        arrange(desc(-min_smoothed_fdr))-> DMRs_homo

write.csv(DMRs_homo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/DMRs_Homogenate.csv")


# plot DMR with lowest min smoothed FDR (gene MLNR)

DMR.plot(ranges = results.ranges, dmr = 67146, CpGs = beta, phen.col = sample_names, 
         what = "Beta", arraytype = "EPIC", genome = "hg19")


# plot DMR with 2 lowest min smoothed FDR (gene GRAMD1C)

DMR.plot(ranges = results.ranges, dmr = 69662, CpGs = beta, phen.col = sample_names, 
         what = "Beta", arraytype = "EPIC", genome = "hg19")

# plot DMR IGF2 

DMR.plot(ranges = results.ranges, dmr = 58302, CpGs = beta, phen.col = sample_names, 
         what = "Beta", arraytype = "EPIC", genome = "hg19")


# Myonuclear DMRs

# create design matrix

dfh %>% 
        filter(condition == "BM" | condition == "PM") -> dfh4


m_vals_myo <- B2M(beta) %>% 
        as.data.frame() %>% 
        dplyr::select(rownames(dfh4)) 


# Sample names
sample_names <- c("1BM", "2BM", "4BM", "5BM", "6BM", "7BM", "8BM", "12BM",
                  "1PM", "2PM", "4PM", "5PM", "6PM", "7PM", "8PM", "12PM")

# Create a vector of group labels
group_labels <- ifelse(grepl("B", sample_names), "Baseline", "Post")

# Create a factor variable with the group labels
group_factor <- factor(group_labels, levels = c("Baseline", "Post"))

# Create the design matrix
design_matrix <- model.matrix(~ group_factor)


# Print the design matrix
design_matrix

annotated_data_myo <- cpg.annotate(datatype = "array", 
                                    object = as.matrix(m_vals_myo), 
                                    what = "M", 
                                    arraytype = "EPIC", 
                                    analysis.type = "differential", 
                                    design_matrix,
                                    fdr = "none", 
                                    pval = 0.05,
                                    coef = "group_factorPost")

dmrcate_results <- dmrcate(annotated_data_myo, 
                           lambda = 1000, 
                           C = 2, 
                           pcutoff = "fdr")


results.ranges <- extractRanges(dmrcate_results)
results.ranges %>% 
        as.data.frame() %>% 
        arrange(desc(-min_smoothed_fdr))-> DMRs_myo

write.csv(DMRs_myo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/DMRs_myonuclei.csv")

# plot DMR with lowest min smoothed FDR (gene MLNR)

DMR.plot(ranges = results.ranges, dmr = 67146, CpGs = beta, phen.col = sample_names, 
         what = "Beta", arraytype = "EPIC", genome = "hg19")





# create plot like in oshlack workflow


# indicate which genome is being used

gen <- "hg19"

# the index of the DMR that we will plot 

dmrIndex <- 58302

# extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))
# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))

# create plot

iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")

# Order annotation by chromosome and position

colnames(Illumina_anno)

annEPIC <- Illumina_anno[order(Illumina_anno$CHR, Illumina_anno$MAPINFO),] %>% 
        as.data.frame() %>% 
        drop_na(MAPINFO)


bValsOrd <- beta[match(annEPIC$Name,rownames(beta)),] 

head(bValsOrd)
# create genomic ranges object from methylation data

cpgData <- GRanges(seqnames=Rle(annEPIC$CHR),
                   ranges=IRanges(start=annEPIC$MAPINFO, end=annEPIC$MAPINFO),
                   strand=Rle(rep("*",nrow(annEPIC))),
                   betas=bValsOrd)
as.data.frame(cpgData)-> x
# extract data on CpGs in DMR
cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])



# methylation data track
methTrack <- DataTrack(range=cpgData, groups=c("Baseline", "Post"),genome = gen,
                       chromosome=chrom, ylim=c(-0.05,1.05),
                       type=c("a","p"), name="DNA Meth.\n(beta value)",
                       background.panel="white", legend=TRUE, cex.title=0.8,
                       cex.axis=0.8, cex.legend=0.8)

# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom,fill="darkred")


tracks <- list(iTrack, gTrack,dmrTrack)

sizes <- c(2,2,2) # set up the relative sizes of the tracks
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))









dmr <- GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end))

overlapping_cpgs <- findOverlaps(cpgData, dmr)


# plot DMR manually


Illumina_anno %>% 
        filter(CHR == "11",
               MAPINFO %in% 2158438:2165961) %>% 
        dplyr::select(cpg = Name, MAPINFO) -> IGF2_dmps

beta[rownames(beta) %in% IGF2_dmps$cpg,] %>% 
        as.data.frame() %>% 
        dplyr::select("1BH", "2BH", "4BH", "5BH", "6BH", "7BH", "8BH", "12BH",
                      "1PH", "2PH", "4PH", "5PH", "6PH", "7PH", "8PH", "12PH",
                      "1BM", "2BM", "4BM", "5BM", "6BM", "7BM", "8BM", "12BM",
                      "1PM", "2PM", "4PM", "5PM", "6PM", "7PM", "8PM", "12PM") %>% 
        mutate(BH = ((`1BH`+ `2BH`+ `4BH`+ `5BH`+ `6BH`+ `7BH`+ `8BH`+ `12BH`)/8),
               PH = ((`1PH`+ `2PH`+ `4PH`+ `5PH`+ `6PH`+ `7PH`+ `8PH`+ `12PH`)/8),
               BM = ((`1BM`+ `2BM`+ `4BM`+ `5BM`+ `6BM`+ `7BM`+ `8BM`+ `12BM`)/8),
               PM = ((`1PM`+ `2PM`+ `4PM`+ `5PM`+ `6PM`+ `7PM`+ `8PM`+ `12PM`)/8)) %>% 
        dplyr::select(BH,PH,BM,PM) %>% 
        rownames_to_column(var = "cpg") %>% 
        merge(.,IGF2_dmps, by = "cpg") %>% 
        arrange(desc(-MAPINFO)) %>% 
        pivot_longer(names_to = "Group", values_to = "beta", cols = 2:5) %>% 
        ggplot(aes(x = MAPINFO, y = beta, color = Group, Group = Group))+
        geom_point()+
        geom_line()





#################################################################################################

### GSEA - KEGG

#################################################################################################

# create annotation file with only one gene name





x <- Illumina_anno %>% 
        mutate(UCSC_RefGene_Name = str_split(UCSC_RefGene_Name, ";", simplify = TRUE)[, 1]) 

unique_gene <- unique(anno$UCSC_RefGene_Name) 
unique_gene <- Filter(function(x) x != "", unique_gene)

gene_cpgs <- list()

# create subset of probes annotated to individual genes

for (i in 1:length(unique_gene)) {
        y <- anno %>% 
                filter(UCSC_RefGene_Name == unique_gene[i]) %>% 
                pull(cpg)
        gene_cpgs[[unique_gene[i]]] <- y
        print(i)
}



# calculate mean delta M for each gene

mean_change_df <- M_change[,39:40]


# Assuming you have a data frame named 'mean_change_df' with columns 'PH_vs_BH' and 'PM_vs_BM', and rownames as probe names

# Initialize an empty data frame to store the results
results_df <- data.frame(
        gene = character(length(unique_gene)),
        mean_change_PH_vs_BH = numeric(length(unique_gene)),
        mean_change_PM_vs_BM = numeric(length(unique_gene)),
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
        results_df[i, "mean_change_PH_vs_BH"] <- sum(gene_mean_change$PH_vs_BH) / length(gene_probes)
        results_df[i, "mean_change_PM_vs_BM"] <- sum(gene_mean_change$PM_vs_BM) / length(gene_probes)
        
        print(i)
}

# View the results data frame
results_df %>% unique(results_df$gene)


# Load the package
library(org.Hs.eg.db)

# Convert gene symbols to Entrez Gene IDs
entrezIDs <- mapIds(org.Hs.eg.db, keys = results_df$gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

as.data.frame(entrezIDs) %>% 
        rownames_to_column(var = "gene") %>% 
        merge(results_df,., by = "gene") %>% 
        filter(entrezIDs != "<NA>") -> change_df

rownames(change_df) <- change_df$entrezIDs
change_df <- change_df[,-4]


# run gage on pathway list

# kegg

data(kegg.sets.hs)


### homogenate

exprsMat <- as.matrix(change_df[2]) # homogenate

kegg_res_homo <- gage(exprs = exprsMat, gsets = kegg.sets.hs, same.dir = TRUE, ref = NULL, samp = NULL)

view(kegg_res_homo$less)                              ### view less methylated KEGG pathways post vs. Baseline
view(kegg_res_homo$greater)                           ### view more methylated KEGG pathways post vs. Baseline

### myonuclei

exprsMat <- as.matrix(change_df[3]) # myonuclei

kegg_res_myo <- gage(exprs = exprsMat, gsets = kegg.sets.hs, same.dir = TRUE, ref = NULL, samp = NULL)

view(kegg_res_myo$less)                              ### view less methylated KEGG pathways post vs. Baseline
view(kegg_res_myo$greater)                           ### view more methylated KEGG pathways post vs. Baseline


###______________________________________________________

### rerun with individual probes and gsameth function

x$cpg = x$Name

m_change_df <- M_change[39:40] %>% 
        rownames_to_column(var = "cpg") %>% 
        merge(., x, by = "cpg") %>% 
        dplyr::select(1:3, UCSC_RefGene_Name)

entrezIDs <- mapIds(org.Hs.eg.db, keys = m_change_df$UCSC_RefGene_Name, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

exprsMat <- as.matrix(m_change_df[2])
rownames(exprsMat) <- entrezIDs

library(missMethyl)

# check for updated KEGG sets
library(gage)
KEGG_new <- kegg.gsets(species = "hsa", id.type = "kegg", check.new=TRUE)



kegg_res_homo <- gsameth(sig.cpg = DMPs_PH_vs_BH$cpg,
                         all.cpg = rownames(beta), 
                         collection = KEGG_new$kg.sets, 
                         array.type = "EPIC")

kegg_res_myo <- gsameth(sig.cpg = DMPs_PM_vs_BM$cpg,
                         all.cpg = rownames(beta), 
                         collection = KEGG_new$kg.sets, 
                         array.type = "EPIC")

# rerun gsameth with only signalling and metabolic kegg pathways

subset <- KEGG_new[["sigmet.idx"]]

kegg.subset = KEGG_new$kg.sets[subset]

kegg.subset_res_homo <- gsameth(sig.cpg = DMPs_PH_vs_BH$cpg,
                         all.cpg = rownames(beta), 
                         collection = kegg.subset, 
                         array.type = "EPIC")

kegg.subset_res_myo <- gsameth(sig.cpg = DMPs_PM_vs_BM$cpg,
                        all.cpg = rownames(beta), 
                        collection = kegg.subset, 
                        array.type = "EPIC")

# add direction of methylation for the entire pathway

rownames(kegg.subset_res_myo)

kegg.subset

entrezIDs <- mapIds(org.Hs.eg.db, keys = results_df$gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
results_df$entrezIDs <- entrezIDs

pathway_dir = data.frame()

for (i in 1:length(KEGG_new$kg.sets)) {
        pathway = names(KEGG_new$kg.sets[i])
        pathway_genes = KEGG_new$kg.sets[[pathway]]
        
        gene_mean_change <- results_df[results_df$entrezIDs %in% pathway_genes,]
        
        pathway_dir[i, "pathway"] <- pathway
        pathway_dir[i, "mean_change_PH_vs_BH"] <- sum(gene_mean_change$mean_change_PH_vs_BH) / length(pathway_genes)
        pathway_dir[i, "mean_change_PM_vs_BM"] <- sum(gene_mean_change$mean_change_PM_vs_BM) / length(pathway_genes)
        
        print(i)
}


kegg.subsets_res_homo <- kegg.subset_res_homo %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir, by = "pathway") %>% 
        mutate(dir = ifelse(mean_change_PH_vs_BH < 0, "hypo", "hyper")) %>% 
        dplyr::select(1,2,3,4,5,6,8) %>% 
        arrange(desc(-P.DE))

write.csv(kegg.subset_res_homo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/kegg_signaling_metabolic_pathways_homogenate.csv")

kegg.subsets_res_myo <- kegg.subset_res_myo %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir, by = "pathway") %>% 
        mutate(dir = ifelse(mean_change_PM_vs_BM < 0, "hypo", "hyper")) %>% 
        dplyr::select(1,2,3,4,5,7,8) %>% 
        arrange(desc(-P.DE))

write.csv(kegg.subset_res_myo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/kegg_signaling_metabolic_pathways_myonuclei.csv")


# add direction to all kegg pathways

kegg_res_homo <- kegg_res_homo %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir, by = "pathway") %>% 
        mutate(dir = ifelse(mean_change_PH_vs_BH < 0, "hypo", "hyper")) %>% 
        dplyr::select(1,2,3,4,5,6,8) %>% 
        arrange(desc(-P.DE))

write.csv(kegg_res_homo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/kegg_all_pathways_homogenate.csv")


kegg_res_myo <- kegg_res_myo %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir, by = "pathway") %>% 
        mutate(dir = ifelse(mean_change_PM_vs_BM < 0, "hypo", "hyper")) %>% 
        dplyr::select(1,2,3,4,5,7,8) %>% 
        arrange(desc(-P.DE))

write.csv(kegg_res_myo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/kegg_all_pathways_myonuclei.csv")


# plot the significant pathways for all 4 groups

### Map KEGG pathways to png



keggresids <- kegg_res_homo %>% 
        filter(P.DE < 0.05) %>%
        mutate(id = substr(pathway, start = 1, stop = 8)) %>% 
        pull(id)

res_df <- matrix(results_df[,2])
rownames(res_df) <- results_df$entrezIDs

scaled_values <- ifelse(res_df < 0,
                        rescale(res_df, to = c(-1, 0)),
                        rescale(res_df, to = c(0, 1)))

rownames(scaled_values)<- rownames(res_df)

setwd("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/GSEA figures/KEGG_all_pathways_homogenate/")
tmp = sapply(keggresids, function(pid) pathview(gene.data = scaled_values, pathway.id = pid, species = "hsa", low = "blue", mid = "grey", high = "yellow"))


keggresids <- kegg_res_myo %>% 
        filter(P.DE < 0.05) %>%
        mutate(id = substr(pathway, start = 1, stop = 8)) %>% 
        pull(id)

res_df <- matrix(results_df[,3])
rownames(res_df) <- results_df$entrezIDs

scaled_values <- ifelse(res_df < 0,
                        rescale(res_df, to = c(-1, 0)),
                        rescale(res_df, to = c(0, 1)))
setwd("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/GSEA figures/KEGG_all_pathways_myonuclei/")
tmp = sapply(keggresids, function(pid) pathview(gene.data = scaled_values, pathway.id = pid, species = "hsa", low = "blue", mid = "grey", high = "yellow"))


keggresids <- kegg.subsets_res_homo %>% 
        filter(P.DE < 0.05) %>%
        mutate(id = substr(pathway, start = 1, stop = 8)) %>% 
        pull(id)

res_df <- matrix(results_df[,2])
rownames(res_df) <- results_df$entrezIDs

scaled_values <- ifelse(res_df < 0,
                        rescale(res_df, to = c(-1, 0)),
                        rescale(res_df, to = c(0, 1)))

rownames(scaled_values)<- rownames(res_df)

setwd("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/GSEA figures/KEGG_sigmet_pathways_homogenate//")
tmp = sapply(keggresids, function(pid) pathview(gene.data = scaled_values, pathway.id = pid, species = "hsa", low = "blue", mid = "grey", high = "yellow"))


keggresids <- kegg.subsets_res_myo %>% 
        filter(P.DE < 0.05) %>%
        mutate(id = substr(pathway, start = 1, stop = 8)) %>% 
        pull(id)

res_df <- matrix(results_df[,3])
rownames(res_df) <- results_df$entrezIDs

scaled_values <- ifelse(res_df < 0,
                        rescale(res_df, to = c(-1, 0)),
                        rescale(res_df, to = c(0, 1)))
setwd("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/GSEA figures/KEGG_sigmet_pathways_myonuclei//")
tmp = sapply(keggresids, function(pid) pathview(gene.data = scaled_values, pathway.id = pid, species = "hsa", low = "blue", mid = "grey", high = "yellow"))



# plot WNT pathway in homo and myo

hsa04310

# get the min and max M value within the pathway
pathway = "hsa04310 Wnt signaling pathway"     # write name of pathway

# find average methylation of the individual genes

results_df %>% 
        filter(entrezIDs %in% kegg.subset[[pathway]]) %>%  pull(entrezIDs) -> x


min(results_df[entrezIDs %in% x,2:3])
max(results_df[entrezIDs %in% x,2:3])




res_df <- matrix(results_df[,2])
rownames(res_df) <- results_df$entrezIDs

scaled_values <- ifelse(res_df < 0,
                        rescale(res_df, to = c(-1, 0), from = c(-0.09638417, 0)),
                        rescale(res_df, to = c(0, 1), from = c(0,0.09510493)))

rownames(scaled_values)<- rownames(res_df)

setwd("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/GSEA figures/")
tmp = sapply("hsa04310", function(pid) pathview(gene.data = scaled_values, pathway.id = "hsa04310", species = "hsa", low = "blue", mid = "grey", high = "yellow"))

res_df <- matrix(results_df[,3])
rownames(res_df) <- results_df$entrezIDs



scaled_values <- ifelse(res_df < 0,
                        rescale(res_df, to = c(-1, 0), from = c(-0.09638417, 0)),
                        rescale(res_df, to = c(0, 1), from = c(0,0.09510493)))

rownames(scaled_values)<- rownames(res_df)

setwd("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/GSEA figures/")
tmp = sapply("hsa04310", function(pid) pathview(gene.data = scaled_values, pathway.id = pid, species = "hsa", low = "blue", mid = "grey", high = "yellow"))


# plot leukocyte transedothelial migration pathway


# get the min and max M value within the pathway
pathway = "hsa04670 Leukocyte transendothelial migration"     # write name of pathway

# find average methylation of the individual genes

results_df %>% 
        filter(entrezIDs %in% kegg.subset[[pathway]]) %>%  pull(entrezIDs) -> x


min(results_df[entrezIDs %in% x,2:3])
max(results_df[entrezIDs %in% x,2:3])




res_df <- matrix(results_df[,2])
rownames(res_df) <- results_df$entrezIDs

scaled_values <- ifelse(res_df < 0,
                        rescale(res_df, to = c(-1, 0), from = c(-0.1313915, 0)),
                        rescale(res_df, to = c(0, 1), from = c(0,0.05768334)))

rownames(scaled_values)<- rownames(res_df)

setwd("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/GSEA figures/")
tmp = sapply("hsa04670", function(pid) pathview(gene.data = scaled_values, pathway.id = "hsa04670", species = "hsa", low = "blue", mid = "grey", high = "yellow"))

res_df <- matrix(results_df[,3])
rownames(res_df) <- results_df$entrezIDs



scaled_values <- ifelse(res_df < 0,
                        rescale(res_df, to = c(-1, 0), from = c(-0.1313915, 0)),
                        rescale(res_df, to = c(0, 1), from = c(0,0.05768334)))

rownames(scaled_values)<- rownames(res_df)

setwd("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/GSEA figures/")
tmp = sapply("hsa04670", function(pid) pathview(gene.data = scaled_values, pathway.id = "hsa04670", species = "hsa", low = "blue", mid = "grey", high = "yellow"))


# check number of overlap between MYO+INT and MYO pathways

merge(kegg_res_homo %>% filter(P.DE < 0.05), 
      kegg_res_myo %>% filter(P.DE < 0.05), by = "pathway") %>% 
        mutate(same = ifelse(dir.x == dir.y, "same", "diff")) %>% 
        filter(same == "diff")



##################################################################################

### KEGG plot

############################################################################################

# make horizontal plot of mean pathway change for top 20 most significant pathways in homogenate

kegg_res_homo %>% 
        arrange(desc(-P.DE)) %>% 
        head(20) -> top_kegg

pathway_levels = top_kegg %>% pull(pathway)

kegg_res_myo %>% 
        filter(pathway %in% top_kegg$pathway) %>% 
        merge(top_kegg, ., by = "pathway") %>% 
        dplyr::select(pathway, PH_vs_BH = mean_change_PH_vs_BH, PM_vs_BM = mean_change_PM_vs_BM) %>% 
        pivot_longer(names_to = "contrast", values_to = "mean_change", cols = 2:3) %>% 
        mutate(pathway = factor(pathway, levels = rev(pathway_levels))) %>% 
        ggplot(aes(y = pathway, x = mean_change, group = contrast, color = contrast))+
        geom_point()+
        geom_vline(xintercept = 0)+
        theme_classic()+
        theme(axis.title.y = element_blank())+
        labs(title = "Topp 20 KEGG pathways in Homogenate")


# repeat for KEGG subsets

# add mean change

kegg.subset_res_homo %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,kegg_res_homo, by = "pathway") %>%
        dplyr::select(1,8,10) %>% 
        arrange(desc(-P.DE.y)) %>% 
        head(20) -> top_kegg.sub

pathway_levels.sub = top_kegg.sub %>% pull(pathway)

kegg_res_myo %>% 
        filter(pathway %in% top_kegg.sub$pathway) %>% 
        merge(top_kegg.sub, ., by = "pathway") %>% 
        dplyr::select(pathway, "MYO+INT" = mean_change_PH_vs_BH, MYO = mean_change_PM_vs_BM) %>% 
        pivot_longer(names_to = "Sample", values_to = "mean_change", cols = 2:3) %>% 
        mutate(pathway = factor(pathway, levels = rev(pathway_levels.sub))) %>% 
        ggplot(aes(y = pathway, x = mean_change, group = Sample, color = Sample))+
        geom_point(size = 1.5)+
        geom_vline(xintercept = 0)+
        theme_classic()+
        theme(axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.position = "right")

#########################################################################################

### GSEA - GO (gene ontology)

#########################################################################################


data("go.sets.hs")

GO_res_homo <- gsameth(sig.cpg = DMPs_PH_vs_BH$cpg,
                      all.cpg = rownames(beta), 
                      collection = go.sets.hs, 
                      array.type = "EPIC")

GO_res_myo <- gsameth(sig.cpg = DMPs_PM_vs_BM$cpg,
        all.cpg = rownames(beta), 
        collection = go.sets.hs, 
        array.type = "EPIC")


# add direction of methylation for GO pathways

pathway_dir_go = data.frame()

for (i in 1:length(go.sets.hs)) {
        pathway = names(go.sets.hs[i])
        pathway_genes = go.sets.hs[[pathway]]
        
        gene_mean_change <- results_df[results_df$entrezIDs %in% pathway_genes,]
        
        pathway_dir_go[i, "pathway"] <- pathway
        pathway_dir_go[i, "mean_change_PH_vs_BH"] <- sum(gene_mean_change$mean_change_PH_vs_BH) / length(pathway_genes)
        pathway_dir_go[i, "mean_change_PM_vs_BM"] <- sum(gene_mean_change$mean_change_PM_vs_BM) / length(pathway_genes)
        
        print(i)
}

data("go.subs.hs")

BP_GO <- names(go.sets$go.sets[go.sets$go.subs$BP])
CC_GO <- names(go.sets.hs[go.subs.hs$CC])
MF_GO <- names(go.sets.hs[go.subs.hs$MF])

GO_res_homo <- GO_res_homo %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(dir = ifelse(mean_change_PH_vs_BH < 0, "hypo", "hyper")) %>% 
        dplyr::select(1,2,3,4,5,6,8) %>% 
        arrange(desc(-P.DE)) %>% 
        mutate(subset = ifelse(pathway %in% BP_GO, "Biological processes", 
                               ifelse(pathway %in% CC_GO, "Cellular components", "Molecular function"))) 

write.csv(GO_res_homo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/GO_all_pathways_homogenate.csv")

GO_res_myo <- GO_res_myo %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(dir = ifelse(mean_change_PM_vs_BM < 0, "hypo", "hyper")) %>% 
        dplyr::select(1,2,3,4,5,7,8) %>% 
        arrange(desc(-P.DE)) %>% 
        mutate(subset = ifelse(pathway %in% BP_GO, "Biological processes", 
                                ifelse(pathway %in% CC_GO, "Cellular components", "Molecular function"))) 


write.csv(GO_res_myo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/GO_all_pathways_myonuclei.csv")



##############################################################################

### GO plot

############################################################################

GO_res_homo %>% 
        arrange(desc(-P.DE)) %>% 
        head(20) -> top_go


pathway_levels = top_go %>% pull(pathway)

GO_res_myo %>% 
        filter(pathway %in% top_go$pathway) %>% 
        merge(top_go, ., by = "pathway") %>% 
        dplyr::select(pathway, PH_vs_BH = mean_change_PH_vs_BH, PM_vs_BM = mean_change_PM_vs_BM, subset = subset.x) %>% 
        pivot_longer(names_to = "contrast", values_to = "mean_change", cols = 2:3) %>% 
        mutate(pathway = factor(pathway, levels = pathway_levels)) %>% 
        ggplot(aes(y = pathway, x = mean_change, group = contrast, color = contrast))+
        geom_point(aes(shape = subset))+
        geom_vline(xintercept = 0)+
        theme_classic()+
        theme(axis.title.y = element_blank())+
        labs(title = "Topp 20 GO pathways in Homogenate")


############################################################################

### GSEA - MsigDB

#################################################################################################

# load MsigDB dataset from Champ package
data("PathwayList")

# convert symbol gene names to entrez ids

MsigDB_champ <- list()

for (i in 1:length(PathwayList)) {
        x <- PathwayList[[i]]
        
        y <- mapIds(org.Hs.eg.db, keys = x, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first") %>% 
                as.data.frame() %>% 
                na.omit() %>% 
                pull(.)
        
        MsigDB_champ[[names(PathwayList[i])]] <- y
        print(i)
        
}

# rewritten with skipping of errors

for (i in 1:length(PathwayList)) {
        tryCatch({
                x <- PathwayList[[i]]
                
                y <- mapIds(org.Hs.eg.db, keys = x, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first") %>% 
                        as.data.frame() %>% 
                        na.omit() %>% 
                        pull(.)
                
                MsigDB_champ[[names(PathwayList[i])]] <- y
                print(i)
        }, error = function(e) {
                cat("Error in iteration", i, ":", conditionMessage(e), "\n")
        })
}



library(missMethyl)




MsigDB_res_homo <- gsameth(sig.cpg = DMPs_PH_vs_BH$cpg,
                       all.cpg = rownames(beta), 
                       collection = MsigDB_champ, 
                       array.type = "EPIC")

MsigDB_res_myo <- gsameth(sig.cpg = DMPs_PM_vs_BM$cpg,
                      all.cpg = rownames(beta), 
                      collection = MsigDB_champ, 
                      array.type = "EPIC")

# arrange MsigDB lists ascending

MsigDB_res_homo <-  MsigDB_res_homo %>% 
        as.data.frame() %>% 
        arrange(desc(-P.DE))

MsigDB_res_myo <-  MsigDB_res_myo %>% 
        as.data.frame() %>% 
        arrange(desc(-P.DE))

write.csv(MsigDB_res_homo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/MsigDB_pathways_homogenate.csv")

write.csv(MsigDB_res_myo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/MsigDB_pathways_myonuclei.csv")

# add mean methylation direction

# add direction of methylation for GO pathways

pathway_dir_MsigDB = data.frame()

for (i in 1:length(MsigDB_champ)) {
        pathway = names(MsigDB_champ[i])
        pathway_genes = MsigDB_champ[[pathway]]
        
        gene_mean_change <- results_df[results_df$entrezIDs %in% pathway_genes,]
        
        pathway_dir_MsigDB[i, "pathway"] <- pathway
        pathway_dir_MsigDB[i, "mean_change_PH_vs_BH"] <- sum(gene_mean_change$mean_change_PH_vs_BH) / length(pathway_genes)
        pathway_dir_MsigDB[i, "mean_change_PM_vs_BM"] <- sum(gene_mean_change$mean_change_PM_vs_BM) / length(pathway_genes)
        
        print(i)
}

######################################################################################################

### MsigDB plot

###########################################################################

MsigDB_res_homo %>% 
        rownames_to_column(var = "pathway") %>% 
        arrange(desc(-P.DE)) %>% 
        head(20) -> top_MsigDB



pathway_levels = top_MsigDB %>% pull(pathway)

top_MsigDB %>% 
        merge(.,pathway_dir_MsigDB, by = "pathway") %>% 
        dplyr::select(pathway, PH_vs_BH = mean_change_PH_vs_BH, PM_vs_BM = mean_change_PM_vs_BM) %>% 
        pivot_longer(names_to = "contrast", values_to = "mean_change", cols = 2:3) %>% 
        mutate(pathway = factor(pathway, levels = rev(pathway_levels))) %>% 
        ggplot(aes(y = pathway, x = mean_change, group = contrast, color = contrast))+
        geom_point()+
        geom_vline(xintercept = 0)+
        theme_classic()+
        theme(axis.title.y = element_blank())+
        labs(title = "Topp 20 MsigDB pathways in Homogenate")

###########################################################################

### check for cell population shift

##########################################################################################

# Use the same method as Voisin et al. 2023 https://onlinelibrary.wiley.com/doi/10.1111/acel.13859

# download the datasets from Rubenstein_skeletal_muscle gene sets in MSigDB



library(xml2)


B_CELLS <- read_xml("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/MsigDB/RUBENSTEIN_SKELETAL_MUSCLE_B_CELLS.v2023.1.Hs.xml")

ENDOTHELIAL_CELLS <- read_xml("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/MsigDB/RUBENSTEIN_SKELETAL_MUSCLE_ENDOTHELIAL_CELLS.v2023.1.Hs.xml")

FAP_CELLS <- read_xml("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/MsigDB/RUBENSTEIN_SKELETAL_MUSCLE_FAP_CELLS.v2023.1.Hs.xml")

FBN1_FAP_CELLS <- read_xml("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/MsigDB/RUBENSTEIN_SKELETAL_MUSCLE_FBN1_FAP_CELLS.v2023.1.Hs.xml")

MYELOID_CELLS <- read_xml("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/MsigDB/RUBENSTEIN_SKELETAL_MUSCLE_MYELOID_CELLS.v2023.1.Hs.xml")

NK_CELLS <- read_xml("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/MsigDB/RUBENSTEIN_SKELETAL_MUSCLE_NK_CELLS.v2023.1.Hs.xml")

PCV_ENDOTHELIAL_CELLS <- read_xml("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/MsigDB/RUBENSTEIN_SKELETAL_MUSCLE_PCV_ENDOTHELIAL_CELLS.v2023.1.Hs.xml")

PERICYTES <- read_xml("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/MsigDB/RUBENSTEIN_SKELETAL_MUSCLE_PERICYTES.v2023.1.Hs.xml")

SATELLITE_CELLS <- read_xml("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/MsigDB/RUBENSTEIN_SKELETAL_MUSCLE_SATELLITE_CELLS.v2023.1.Hs.xml")

SMOOTH_MUSCLE_CELLS <- read_xml("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/MsigDB/RUBENSTEIN_SKELETAL_MUSCLE_SMOOTH_MUSCLE_CELLS.v2023.1.Hs.xml")

T_CELLS <- read_xml("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/MsigDB/RUBENSTEIN_SKELETAL_MUSCLE_T_CELLS.v2023.1.Hs.xml")



cell_populations <- list()

cell_populations["B_CELLS"] <- xml_attrs(xml_child(B_CELLS, 1))[["MEMBERS_EZID"]]
cell_populations["ENDOTHELIAL_CELLS"] <- xml_attrs(xml_child(ENDOTHELIAL_CELLS, 1))[["MEMBERS_EZID"]]
cell_populations["FAP_CELLS"] <- xml_attrs(xml_child(FAP_CELLS, 1))[["MEMBERS_EZID"]]
cell_populations["FBN1_FAP_CELLS"] <- xml_attrs(xml_child(FBN1_FAP_CELLS, 1))[["MEMBERS_EZID"]]
cell_populations["MYELOID_CELLS"] <- xml_attrs(xml_child(MYELOID_CELLS, 1))[["MEMBERS_EZID"]]
cell_populations["NK_CELLS"] <- xml_attrs(xml_child(NK_CELLS, 1))[["MEMBERS_EZID"]]
cell_populations["PCV_ENDOTHELIAL_CELLS"] <- xml_attrs(xml_child(PCV_ENDOTHELIAL_CELLS, 1))[["MEMBERS_EZID"]]
cell_populations["PERICYTES"] <- xml_attrs(xml_child(PERICYTES, 1))[["MEMBERS_EZID"]]
cell_populations["SATELLITE_CELLS"] <- xml_attrs(xml_child(SATELLITE_CELLS, 1))[["MEMBERS_EZID"]]
cell_populations["SMOOTH_MUSCLE_CELLS"] <- xml_attrs(xml_child(SMOOTH_MUSCLE_CELLS, 1))[["MEMBERS_EZID"]]
cell_populations["T_CELLS"] <- xml_attrs(xml_child(T_CELLS, 1))[["MEMBERS_EZID"]]


split_string <- function(x) {
        unlist(strsplit(x, split = ",", fixed = TRUE))
}

cell_populations <- lapply(cell_populations, split_string)


# run over-representation analysis for homogenate 

cell_population_genes_homo <- gsameth(sig.cpg = DMPs_PH_vs_BH$cpg,
                       all.cpg = rownames(beta), 
                       collection = cell_populations, 
                       array.type = "EPIC")


# run over-representation analysis for myonuclei 

cell_population_genes_myo <- gsameth(sig.cpg = DMPs_PM_vs_BM$cpg,
                                      all.cpg = rownames(beta), 
                                      collection = cell_populations, 
                                      array.type = "EPIC")


write.csv(cell_population_genes_homo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/cell_population_genes_homogenate.csv")

write.csv(cell_population_genes_myo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/cell_population_genes_myonuclei.csv")

# cell population direction

pathway_dir_cells = data.frame()

for (i in 1:length(cell_populations)) {
        pathway = names(cell_populations[i])
        pathway_genes = cell_populations[[pathway]]
        
        gene_mean_change <- results_df[results_df$entrezIDs %in% pathway_genes,]
        
        pathway_dir_cells[i, "cell_population"] <- pathway
        pathway_dir_cells[i, "mean_change_PH_vs_BH"] <- sum(gene_mean_change$mean_change_PH_vs_BH) / length(pathway_genes)
        pathway_dir_cells[i, "mean_change_PM_vs_BM"] <- sum(gene_mean_change$mean_change_PM_vs_BM) / length(pathway_genes)
        
        print(i)
}

# plot cell populations

set.seed(1)
cell_population_genes_homo %>% 
        rownames_to_column(var = "cell_population") %>% 
        merge(.,pathway_dir_cells, by = "cell_population") %>% 
        arrange(desc(-P.DE)) %>%     #pull(cell_population) -> cell_levels
        dplyr::select(cell_population, PH_vs_BH = mean_change_PH_vs_BH, PM_vs_BM = mean_change_PM_vs_BM, pval_homo = P.DE) %>% 
        merge(.,rownames_to_column(cell_population_genes_myo, var = "cell_population"), by = "cell_population" ) %>% 
        dplyr::select(1:4,pval_myo = P.DE) %>% 
        pivot_longer(names_to = "contrast", values_to = "mean_change", cols = 2:3) %>% 
        mutate(pval = ifelse(contrast == "PH_vs_BH", pval_homo, pval_myo)) %>% 
        mutate(pval = ifelse(pval <0.06, paste("p =",round(pval, 2)), "")) %>% 
        mutate(cell_population = factor(cell_population, levels = rev(cell_levels))) %>% 
        ggplot(aes(y = cell_population, x = mean_change, group = contrast, color = contrast))+
        geom_point(size = 2)+
        geom_vline(xintercept = 0)+
        theme_classic(base_size = 20)+
        scale_y_discrete(labels = gsub("_"," ", rev(cell_levels)))+
        theme(axis.title.y = element_blank(), 
              axis.title.x = element_blank())+
        geom_text_repel(aes(label = pval))+
        labs(color = "Sample")+
        scale_color_discrete(labels = c("MYO+INT","MYO"))
        
        
        
        
        
#############################################################################################

### search code

#############################################################################################


# find genes in pathway

_________________________________________________
# MsigDB

pathway = "KEGG_FOCAL_ADHESION"      # write name of pathway

# find average methylation of the individual genes

results_df %>% 
        filter(gene %in% PathwayList[[pathway]])


_________________________________________________
# KEGG

pathway = "hsa04150 mTOR signaling pathway"     # write name of pathway

# find average methylation of the individual genes

results_df %>% 
        filter(entrezIDs %in% kegg.subset[[pathway]]) %>% 
        mutate(difference = mean_change_PM_vs_BM-mean_change_PH_vs_BH) %>% 
        arrange(desc(abs(mean_change_PM_vs_BM))) %>% 
        summarise(mean_homo = mean(mean_change_PH_vs_BH),
                  mean_myo = mean(mean_change_PM_vs_BM))



pathway = "hsa04020 Calcium signaling pathway"     # write name of pathway

# find average methylation of the individual genes

results_df %>% 
        filter(entrezIDs %in% kegg.subset[[pathway]]) %>% 
        mutate(difference = mean_change_PM_vs_BM-mean_change_PH_vs_BH) %>% 
        arrange(desc(abs(difference))) %>% 
        summarise(mean_homo = mean(mean_change_PH_vs_BH),
                  mean_myo = mean(mean_change_PM_vs_BM))

pathway = "hsa04310 Wnt signaling pathway"     # write name of pathway

# find average methylation of the individual genes

results_df %>% 
        filter(entrezIDs %in% kegg.subset[[pathway]]) %>% 
        mutate(difference = mean_change_PM_vs_BM-mean_change_PH_vs_BH) %>% 
        arrange(desc(abs(difference))) %>% 
        summarise(mean_homo = mean(mean_change_PH_vs_BH),
                  mean_myo = mean(mean_change_PM_vs_BM))

_________________________________________________
# GO

pathway = "GO:0048468 cell development"       # write name of pathway

# find average methylation of the individual genes

results_df %>% 
        filter(entrezIDs %in% go.sets.hs[[pathway]]) 


_________________________________________________

# find probes in genes

_________________________________________________

Illumina_anno <- Illumina_anno %>% 
        mutate("cpg" = IlmnID) 

gene_name <- "GRK5"       # write gene name

Probes <- anno %>% 
        filter(UCSC_RefGene_Name == gene_name) 

x = c("cg25431166", "cg11866473") # extra probes if misidentified

Chromosome <- merge(Probes, Illumina_anno, by = "cpg") %>% 
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
        geom_vline(xintercept = tss_data %>% 
                           as.data.frame() %>% 
                           filter(external_gene_name == gene_name) %>% 
                           pull(transcription_start_site) , color = "red")

# add TSS to figure

BiocManager::install("biomaRt")

library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

tss_data <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "transcription_start_site", "chromosome_name", "strand"),
                  mart = ensembl)

tss_data %>% 
        as.data.frame() %>% 
        filter(external_gene_name == gene_name) %>% 
        pull(transcription_start_site) %>% mean()


##################################################################

### linear regression of individual genes in Wnt pathway

###############################################################################

library(lme4)
library(lmerTest)
library(emmeans)

# identify genes in pathway and probes in genes


pathway = "hsa04310 Wnt signaling pathway"    

# find average methylation of the individual genes

results_df %>% 
        filter(entrezIDs %in% kegg.subset[[pathway]]) %>% 
        mutate(difference = mean_change_PM_vs_BM-mean_change_PH_vs_BH) %>% 
        arrange(desc(abs(difference))) -> Wnt



for (i in 1:length(Wnt$gene)) {
        
        tryCatch({
        anno %>% 
                filter(UCSC_RefGene_Name == Wnt[i,1]) %>% 
                pull(cpg) -> x
        
        M_change[,1:32] %>% 
                rownames_to_column(var = "cpg") %>% 
                filter(cpg %in% x) %>% 
                pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
                mutate(Condition = substr_right(FP, 2),
                       ID = as.factor(str_sub(FP, end = -3)),
                       Timepoint = as.factor(str_sub(Condition, end = 1)),
                       Sample = as.factor(substr_right(Condition, 1))) %>% 
                dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
                as.data.frame() %>% 
                lmer(M_value ~ cpg +  Timepoint*Sample + (1|ID), data = .) -> rmaModel
        
        model <- coef(summary(rmaModel)) %>%  tail(3) %>% head(1) %>% as.data.frame()
        
        Wnt[i,"p.val(timepoint)"] <- model$`Pr(>|t|)`
        print(i)
        }, error= function(e) {
                cat("Error in iteration", i, ":", conditionMessage(e), "\n")
        })
        
        
}

# get list of genes that are significnat for time and sample

Wnt %>% 
        dplyr::select(1:5, p.val = 6) %>% 
        filter(p.val < 0.05) %>% pull(gene)

# plot emmeans for these genes

Wnt_sig = c("FRAT2",   "DVL3" ,   "WNT1"  ,  "AXIN1"   ,"CCDC88C", "CSNK2B" , "LGR6"  ,  "CHD8"  ,  "CTBP1")

# supplementary file
sig = data.frame()

for (i in 1:length(Wnt_sig)) {
        
        anno %>% 
                filter(UCSC_RefGene_Name == Wnt_sig[i]) %>% 
                pull(cpg) -> x
        
        M_change[,1:32] %>% 
                rownames_to_column(var = "cpg") %>% 
                filter(cpg %in% x) %>% 
                pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
                mutate(Condition = substr_right(FP, 2),
                       ID = as.factor(str_sub(FP, end = -3)),
                       Timepoint = as.factor(str_sub(Condition, end = 1)),
                       Sample = as.factor(substr_right(Condition, 1))) %>% 
                dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
                as.data.frame() %>% 
                lmer(M_value ~ cpg +  Timepoint*Sample + (1|ID), data = .) -> rmaModel
        
        est <- emmeans(rmaModel, specs = ~Timepoint|Sample, pbkrtest.limit = 4992, lmerTest.limit = 4992)
        
        sig <- summary(est) %>% 
                dplyr::select(Sample, Timepoint, emmean, SE, df, lower.CL, upper.CL) %>% 
                mutate(gene = Wnt_sig[i]) %>% 
                rbind(sig,.)
        
        print(i)
}

# fix sig dataframe and save to csv

sig %>% 
        mutate(Sample = ifelse(Sample == "H", "MYO+INT", "MYO"),
               Timepoint = ifelse(Timepoint == "B", "Baseline", "Post")) %>% 
        dplyr::select(gene, 1:7) %>% 
        write.csv(file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/Wnt_emmeans.csv")

# FRAT2

anno %>% 
        filter(UCSC_RefGene_Name == Wnt_sig[1]) %>% 
        pull(cpg) -> x
        

M_change[,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% x) %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Condition = substr_right(FP, 2),
               ID = as.factor(str_sub(FP, end = -3)),
               Timepoint = as.factor(str_sub(Condition, end = 1)),
               Sample = as.factor(substr_right(Condition, 1))) %>% 
        dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
        as.data.frame() %>% 
        lmer(M_value ~ cpg +  Timepoint*Sample + (1|ID), data = .) -> rmaModel


qqnorm(resid(rmaModel));qqline(resid(rmaModel))

plot(rmaModel)

est <- emmeans(rmaModel, specs = ~Timepoint|Sample)



# plot the estimated marginal means 


p1 <- est %>%
        data.frame() %>%
        mutate(Timepoint = factor(Timepoint, levels = c("B", "P"))) %>%
        ggplot(aes(Timepoint, emmean, group = Sample, color = Sample)) + 
        geom_line(position = position_dodge(width = 0.2), size = 1.2) +
        geom_point(position = position_dodge(width = 0.2), size = 4)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      position = position_dodge(width = 0.2), 
                      width = 0.2, size = 1.2)+
        theme_classic(base_size = 20)+
        scale_x_discrete(labels = c("Baseline", "Post"))+
        scale_color_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())
       


# DVL3
anno %>% 
        filter(UCSC_RefGene_Name == Wnt_sig[2]) %>% 
        pull(cpg) -> x


M_change[,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% x) %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Condition = substr_right(FP, 2),
               ID = as.factor(str_sub(FP, end = -3)),
               Timepoint = as.factor(str_sub(Condition, end = 1)),
               Sample = as.factor(substr_right(Condition, 1))) %>% 
        dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
        as.data.frame() %>% 
        lmer(M_value ~ cpg +  Timepoint*Sample + (1|ID), data = .) -> rmaModel


qqnorm(resid(rmaModel));qqline(resid(rmaModel))

plot(rmaModel)

est <- emmeans(rmaModel, specs = ~Timepoint|Sample)

# plot the estimated marginal means 


p2 <- est %>%
        data.frame() %>%
        mutate(Timepoint = factor(Timepoint, levels = c("B", "P"))) %>%
        ggplot(aes(Timepoint, emmean, group = Sample, color = Sample)) + 
        geom_line(position = position_dodge(width = 0.2), size = 1.2) +
        geom_point(position = position_dodge(width = 0.2), size = 4)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      position = position_dodge(width = 0.2), 
                      width = 0.2, size = 1.2)+
        theme_classic(base_size = 20)+
        scale_x_discrete(labels = c("Baseline", "Post"))+
        scale_color_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())
# WNT1
anno %>% 
        filter(UCSC_RefGene_Name == Wnt_sig[3]) %>% 
        pull(cpg) -> x


M_change[,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% x) %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Condition = substr_right(FP, 2),
               ID = as.factor(str_sub(FP, end = -3)),
               Timepoint = as.factor(str_sub(Condition, end = 1)),
               Sample = as.factor(substr_right(Condition, 1))) %>% 
        dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
        as.data.frame() %>% 
        lmer(M_value ~ cpg +  Timepoint*Sample + (1|ID), data = .) -> rmaModel


qqnorm(resid(rmaModel));qqline(resid(rmaModel))

plot(rmaModel)

est <- emmeans(rmaModel, specs = ~Timepoint|Sample)

# plot the estimated marginal means 


p3 <- est %>%
        data.frame() %>%
        mutate(Timepoint = factor(Timepoint, levels = c("B", "P"))) %>%
        ggplot(aes(Timepoint, emmean, group = Sample, color = Sample)) + 
        geom_line(position = position_dodge(width = 0.2), size = 1.2) +
        geom_point(position = position_dodge(width = 0.2), size = 4)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      position = position_dodge(width = 0.2), 
                      width = 0.2, size = 1.2)+
        theme_classic(base_size = 20)+
        scale_x_discrete(labels = c("Baseline", "Post"))+
        scale_color_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())
# AXIN1
anno %>% 
        filter(UCSC_RefGene_Name == Wnt_sig[4]) %>% 
        pull(cpg) -> x


M_change[,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% x) %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Condition = substr_right(FP, 2),
               ID = as.factor(str_sub(FP, end = -3)),
               Timepoint = as.factor(str_sub(Condition, end = 1)),
               Sample = as.factor(substr_right(Condition, 1))) %>% 
        dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
        as.data.frame() %>% 
        lmer(M_value ~ cpg +  Timepoint*Sample + (1|ID), data = .) -> rmaModel


qqnorm(resid(rmaModel));qqline(resid(rmaModel))

plot(rmaModel)

est <- emmeans(rmaModel, specs = ~Timepoint|Sample)

# plot the estimated marginal means 


p4 <- est %>%
        data.frame() %>%
        mutate(Timepoint = factor(Timepoint, levels = c("B", "P"))) %>%
        ggplot(aes(Timepoint, emmean, group = Sample, color = Sample)) + 
        geom_line(position = position_dodge(width = 0.2), size = 1.2) +
        geom_point(position = position_dodge(width = 0.2), size = 4)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      position = position_dodge(width = 0.2), 
                      width = 0.2, size = 1.2)+
        theme_classic(base_size = 20)+
        scale_x_discrete(labels = c("Baseline", "Post"))+
        scale_color_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())
# CCDC88C
anno %>% 
        filter(UCSC_RefGene_Name == Wnt_sig[5]) %>% 
        pull(cpg) -> x


M_change[,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% x) %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Condition = substr_right(FP, 2),
               ID = as.factor(str_sub(FP, end = -3)),
               Timepoint = as.factor(str_sub(Condition, end = 1)),
               Sample = as.factor(substr_right(Condition, 1))) %>% 
        dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
        as.data.frame() %>% 
        lmer(M_value ~ cpg +  Timepoint*Sample + (1|ID), data = .) -> rmaModel


qqnorm(resid(rmaModel));qqline(resid(rmaModel))

plot(rmaModel)

est <- emmeans(rmaModel, specs = ~Timepoint|Sample, pbkrtest.limit = 4992, lmerTest.limit = 4992)

# plot the estimated marginal means 


p5 <- est %>%
        data.frame() %>%
        mutate(Timepoint = factor(Timepoint, levels = c("B", "P"))) %>%
        ggplot(aes(Timepoint, emmean, group = Sample, color = Sample)) + 
        geom_line(position = position_dodge(width = 0.2), size = 1.2) +
        geom_point(position = position_dodge(width = 0.2), size = 4)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      position = position_dodge(width = 0.2), 
                      width = 0.2, size = 1.2)+
        theme_classic(base_size = 20)+
        scale_x_discrete(labels = c("Baseline", "Post"))+
        scale_color_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())
# CSNK2B
anno %>% 
        filter(UCSC_RefGene_Name == Wnt_sig[6]) %>% 
        pull(cpg) -> x


M_change[,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% x) %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Condition = substr_right(FP, 2),
               ID = as.factor(str_sub(FP, end = -3)),
               Timepoint = as.factor(str_sub(Condition, end = 1)),
               Sample = as.factor(substr_right(Condition, 1))) %>% 
        dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
        as.data.frame() %>% 
        lmer(M_value ~ cpg +  Timepoint*Sample + (1|ID), data = .) -> rmaModel


qqnorm(resid(rmaModel));qqline(resid(rmaModel))

plot(rmaModel)

est <- emmeans(rmaModel, specs = ~Timepoint|Sample)

# plot the estimated marginal means 


p6 <- est %>%
        data.frame() %>%
        mutate(Timepoint = factor(Timepoint, levels = c("B", "P"))) %>%
        ggplot(aes(Timepoint, emmean, group = Sample, color = Sample)) + 
        geom_line(position = position_dodge(width = 0.2), size = 1.2) +
        geom_point(position = position_dodge(width = 0.2), size = 4)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      position = position_dodge(width = 0.2), 
                      width = 0.2, size = 1.2)+
        theme_classic(base_size = 20)+
        scale_x_discrete(labels = c("Baseline", "Post"))+
        scale_color_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())
# LGR6
anno %>% 
        filter(UCSC_RefGene_Name == Wnt_sig[7]) %>% 
        pull(cpg) -> x


M_change[,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% x) %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Condition = substr_right(FP, 2),
               ID = as.factor(str_sub(FP, end = -3)),
               Timepoint = as.factor(str_sub(Condition, end = 1)),
               Sample = as.factor(substr_right(Condition, 1))) %>% 
        dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
        as.data.frame() %>% 
        lmer(M_value ~ cpg +  Timepoint*Sample + (1|ID), data = .) -> rmaModel


qqnorm(resid(rmaModel));qqline(resid(rmaModel))

plot(rmaModel)

est <- emmeans(rmaModel, specs = ~Timepoint|Sample)

# plot the estimated marginal means 


p7 <- est %>%
        data.frame() %>%
        mutate(Timepoint = factor(Timepoint, levels = c("B", "P"))) %>%
        ggplot(aes(Timepoint, emmean, group = Sample, color = Sample)) + 
        geom_line(position = position_dodge(width = 0.2), size = 1.2) +
        geom_point(position = position_dodge(width = 0.2), size = 4)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      position = position_dodge(width = 0.2), 
                      width = 0.2, size = 1.2)+
        theme_classic(base_size = 20)+
        scale_x_discrete(labels = c("Baseline", "Post"))+
        scale_color_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())
# CHD8
anno %>% 
        filter(UCSC_RefGene_Name == Wnt_sig[8]) %>% 
        pull(cpg) -> x


M_change[,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% x) %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Condition = substr_right(FP, 2),
               ID = as.factor(str_sub(FP, end = -3)),
               Timepoint = as.factor(str_sub(Condition, end = 1)),
               Sample = as.factor(substr_right(Condition, 1))) %>% 
        dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
        as.data.frame() %>% 
        lmer(M_value ~ cpg +  Timepoint*Sample + (1|ID), data = .) -> rmaModel


qqnorm(resid(rmaModel));qqline(resid(rmaModel))

plot(rmaModel)

est <- emmeans(rmaModel, specs = ~Timepoint|Sample)

# plot the estimated marginal means 


p8 <- est %>%
        data.frame() %>%
        mutate(Timepoint = factor(Timepoint, levels = c("B", "P"))) %>%
        ggplot(aes(Timepoint, emmean, group = Sample, color = Sample)) + 
        geom_line(position = position_dodge(width = 0.2), size = 1.2) +
        geom_point(position = position_dodge(width = 0.2), size = 4)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      position = position_dodge(width = 0.2), 
                      width = 0.2, size = 1.2)+
        theme_classic(base_size = 20)+
        scale_x_discrete(labels = c("Baseline", "Post"))+
        scale_color_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())
# CTBP1
anno %>% 
        filter(UCSC_RefGene_Name == Wnt_sig[9]) %>% 
        pull(cpg) -> x


M_change[,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% x) %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Condition = substr_right(FP, 2),
               ID = as.factor(str_sub(FP, end = -3)),
               Timepoint = as.factor(str_sub(Condition, end = 1)),
               Sample = as.factor(substr_right(Condition, 1))) %>% 
        dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
        as.data.frame() %>% 
        lmer(M_value ~ cpg +  Timepoint*Sample + (1|ID), data = .) -> rmaModel


qqnorm(resid(rmaModel));qqline(resid(rmaModel))

plot(rmaModel)

est <- emmeans(rmaModel, specs = ~Timepoint|Sample)

# plot the estimated marginal means 


p9 <- est %>%
        data.frame() %>%
        mutate(Timepoint = factor(Timepoint, levels = c("B", "P"))) %>%
        ggplot(aes(Timepoint, emmean, group = Sample, color = Sample)) + 
        geom_line(position = position_dodge(width = 0.2), size = 1.2) +
        geom_point(position = position_dodge(width = 0.2), size = 4)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      position = position_dodge(width = 0.2), 
                      width = 0.2, size = 1.2)+
        theme_classic(base_size = 20)+
        scale_x_discrete(labels = c("Baseline", "Post"))+
        scale_color_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())


library(ggpubr)
fig <- ggarrange(p1+rremove("ylab"),
          p2+rremove("ylab"),
          p3+rremove("ylab"),
          p4+rremove("ylab"),
          p5+rremove("ylab"),
          p6+rremove("ylab"),
          p7+rremove("ylab"),
          p8+rremove("ylab"),
          p9+rremove("ylab"), 
          nrow = 1, labels = c("C","D","E","F","G","H","I","J","K"), common.legend = TRUE, 
          hjust = 0.09, font.label = list(size = 25))

annotate_figure(fig, left = "emmenas M-value")



# run WNT9a

anno %>% 
        filter(UCSC_RefGene_Name == "TNF") %>% 
        pull(cpg) -> x

M_change[,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% x) %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Condition = substr_right(FP, 2),
               ID = as.factor(str_sub(FP, end = -3)),
               Timepoint = as.factor(str_sub(Condition, end = 1)),
               Sample = as.factor(substr_right(Condition, 1))) %>% 
        dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
        as.data.frame() %>% 
        lmer(M_value ~ cpg +  Timepoint*Sample + (1|ID), data = .) -> rmaModel

model <- coef(summary(rmaModel)) %>%  tail(3) %>% head(1) %>% as.data.frame()

qqnorm(resid(rmaModel));qqline(resid(rmaModel))

plot(rmaModel)

est <- emmeans(rmaModel, specs = ~Timepoint|Sample)

# plot the estimated marginal means 


est %>%
        data.frame() %>%
        mutate(Timepoint = factor(Timepoint, levels = c("B", "P"))) %>%
        ggplot(aes(Timepoint, emmean, group = Sample, color = Sample)) + 
        geom_line(position = position_dodge(width = 0.2), size = 1.2) +
        geom_point(position = position_dodge(width = 0.2), size = 4)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      position = position_dodge(width = 0.2), 
                      width = 0.2, size = 1.2)+
        theme_classic(base_size = 20)+
        scale_x_discrete(labels = c("Baseline", "Post"))+
        scale_color_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())

# rerun wth only island probes

anno %>% 
        filter(UCSC_RefGene_Name == "WNT9A" & Relation_to_Island == "Island") %>% 
        pull(cpg) -> x

x = "cg21467614"

M_change[,1:32] %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% x) %>% 
        pivot_longer(names_to = "FP", values_to = "M_value", cols = 2:33) %>% 
        mutate(Condition = substr_right(FP, 2),
               ID = as.factor(str_sub(FP, end = -3)),
               Timepoint = as.factor(str_sub(Condition, end = 1)),
               Sample = as.factor(substr_right(Condition, 1))) %>% 
        dplyr::select(cpg, ID, Timepoint, Sample, M_value) %>% 
        as.data.frame() %>% 
        lmer(M_value ~ #cpg +  
                     Timepoint*Sample + (1|ID), data = .) -> rmaModel

model <- coef(summary(rmaModel)) %>%  tail(3) %>% head(1) %>% as.data.frame()

qqnorm(resid(rmaModel));qqline(resid(rmaModel))

plot(rmaModel)

est <- emmeans(rmaModel, specs = ~Timepoint|Sample)

# plot the estimated marginal means 


est %>%
        data.frame() %>%
        mutate(Timepoint = factor(Timepoint, levels = c("B", "P"))) %>%
        ggplot(aes(Timepoint, emmean, group = Sample, color = Sample)) + 
        geom_line(position = position_dodge(width = 0.2), size = 1.2) +
        geom_point(position = position_dodge(width = 0.2), size = 4)+
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      position = position_dodge(width = 0.2), 
                      width = 0.2, size = 1.2)+
        theme_classic(base_size = 20)+
        scale_x_discrete(labels = c("Baseline", "Post"))+
        scale_color_discrete(labels = c("MYO+INT","MYO"))+
        theme(axis.title.x = element_blank())


Illumina_anno  %>% 
        select(Name, 16, 21,22) %>% 
        filter(UCSC_RefGene_Name == "WNT9A") %>% 
        pull(cpg) -> x

DMPs_PH_vs_BH %>% 
        filter(cpg %in% x)

DMPs_PM_vs_BM %>% 
        filter(cpg %in% x) %>% 
        merge(.,anno, by = "cpg")

anno %>% 
        filter(UCSC_RefGene_Name %in% c("TNF", "DIF", "TNF-alpha", "TNFA", "TNFSF2", "TNLG1F"))
##############################################################################

### epigenetic age

################################################################

### MsetExProbes dataset from "Oshlack workflow - filtering and QC"





epiage_horvath <- agep(beta, coeff=NULL, method="horvath")

epiage_horvath %>% 
        rownames_to_column(var = "participant") -> x

epiage <- dfh %>% 
        rownames_to_column(var = "participant") %>% 
        merge(.,x, by = "participant")


### boxplot of ages for the groups

epiage %>% 
        ggplot(aes(x = condition, y = horvath.age, fill = condition))+
        geom_boxplot()+
        theme_classic()


###     MEAT 2.0 muscle tissue clock



library(SummarizedExperiment); library(MEAT)




### create Summarized experiment element for the epiage_estimation 
### as seen below it can be run with and without a pheno dataframe

cancer <- SummarizedExperiment(assays = list(beta = samp),
                               colData = pheno)                 

beta <- SummarizedExperiment(assays = list(beta = beta))


### then the CpGs need to be "cleaned" so the dataset only contains the relevant 18747 CpG sites

beta_clean <- clean_beta(SE = beta,
                           version = "MEAT2.0")

### calibrate beta values 



beta_clean_calibrated <- BMIQcalibration(beta_clean)


### then you estimate the epigenetic age with this function, where the "age_col_name" is optional, 
###     but reqired if you wish to get difference and residuals between chronological age and predicted age


epiage_meat <- epiage_estimation(SE = beta_clean,
                                 age_col_name = NULL,
                                 version = "MEAT2.0")

epiage_meat_calibrated <- epiage_estimation(SE = beta_clean_calibrated,
                                 age_col_name = NULL,
                                 version = "MEAT2.0")

### get only the age estimations


DNAmage <- as.data.frame(epiage_meat$DNAmage)

DNAmage_calibrated <- as.data.frame(epiage_meat_calibrated$DNAmage)

### add to earlyer created participant data with horvath clock estimations



epiage <- cbind(epiage, DNAmage)

epiage <- cbind(epiage, DNAmage_calibrated) 

# plot epiage

epiage %>% 
        ggplot(aes(x = condition, y = `epiage_meat$DNAmage`+20, fill = condition))+
        geom_boxplot()+
        theme_classic()

epiage %>% 
        ggplot(aes(x = condition, y = `epiage_meat_calibrated$DNAmage`, fill = condition))+
        geom_boxplot()+
        theme_classic()

BH = epiage %>% filter(condition == "BH") %>% pull(`epiage_meat_calibrated$DNAmage`)
PH = epiage %>% filter(condition == "PH") %>% pull(`epiage_meat_calibrated$DNAmage`)
BM = epiage %>% filter(condition == "BM") %>% pull(`epiage_meat_calibrated$DNAmage`)
PM = epiage %>% filter(condition == "PM") %>% pull(`epiage_meat_calibrated$DNAmage`)

        
t.test(x = PH, y = BH, paired = TRUE)
t.test(x = PM, y = BM, paired = TRUE)







x = c("cg25431166", "cg11866473", "cg08114542", "cg05600174", "cg05386815", "cg04036064",
      "cg10656332", "cg06262380")

NMRK1 <- anno %>% 
        filter(UCSC_RefGene_Name == "NMRK1" | cpg %in% x)

DMPs_PH_vs_BH %>% 
        filter(cpg %in% x)

DMPs_PM_vs_BM %>% 
        filter(cpg %in% x)

y = c("cg14128911","cg19433767", "cg22782812", "cg01074840", "cg13811969", "cg21690947")

NMRK2 <- anno %>% 
        filter(UCSC_RefGene_Name == "NMRK2" | cpg %in% y)

DMPs_PH_vs_BH %>% 
        filter(cpg %in% y)

DMPs_PM_vs_BM %>% 
        filter(cpg %in% y)

# not changed with training

DMPs_BM_vs_BH %>% 
        filter(cpg %in% x)

DMPs_PM_vs_PH %>% 
        filter(cpg %in% x)

DMPs_BM_vs_BH %>% 
        filter(cpg %in% y)

DMPs_PM_vs_PH %>% 
        filter(cpg %in% y)


# plot NMRK1 probes MYO+int vs. MYO

M_change %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% x) %>% 
        dplyr::select(cpg, BM_vs_BH, PM_vs_PH) %>% 
        merge(., Illumina_anno, by = "cpg") %>% 
        dplyr::select(cpg, BM_vs_BH, PM_vs_PH, MAPINFO,  Relation_to_UCSC_CpG_Island, CHR) %>% 
        pivot_longer(names_to = "contrast", values_to = "mean_m", cols = 2:3) %>% 
        arrange(desc(MAPINFO)) %>% 
        ggplot(aes(x = MAPINFO, y = mean_m, color = contrast, group = contrast))+
        geom_hline(yintercept = 0)+
        geom_point(aes(shape = Relation_to_UCSC_CpG_Island), size = 3)+
        geom_smooth(method = "loess", alpha = 0.5, se = FALSE)+
        labs(title = paste("NMRK1",";", "N Probes =", 8, ";","CHR", 9),
             y = "mean M-value difference",
             x = "nucleotide")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90))+
        scale_x_continuous(n.breaks = 20)


M_change %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% y) %>% 
        dplyr::select(cpg, BM_vs_BH, PM_vs_PH) %>% 
        merge(., Illumina_anno, by = "cpg") %>% 
        dplyr::select(cpg, BM_vs_BH, PM_vs_PH, MAPINFO,  Relation_to_UCSC_CpG_Island, CHR) %>% 
        pivot_longer(names_to = "contrast", values_to = "mean_m", cols = 2:3) %>% 
        arrange(desc(MAPINFO)) %>% 
        ggplot(aes(x = MAPINFO, y = mean_m, color = contrast, group = contrast))+
        geom_hline(yintercept = 0)+
        geom_point(aes(shape = Relation_to_UCSC_CpG_Island), size = 3)+
        geom_smooth(method = "loess", alpha = 0.5, se = FALSE)+
        labs(title = paste("NMRK2",";", "N Probes =", 6, ";","CHR", 19),
             y = "mean M-value difference",
             x = "nucleotide")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90))+
        scale_x_continuous(n.breaks = 20)








##############################################################################

### PCA and gene correlation

################################################################


# PCA 1 correlates with cell population

pca.out <- prcomp(t(M_change[,1:32]), scale. = FALSE)

loadings <- pca.out$rotation[,1]

sorted_loadings <- sort((loadings), decreasing = TRUE) 
top_loadings <- names(sorted_loadings)[1:100]

# load annotation df

anno <- readRDS("anno.RDATA")

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

# isolate the abs arranged list og PC1 drivers




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
        
        
        
        
        
        
        
        
        



data("go.subs.hs")

BP_GO <- names(go.sets$go.sets[go.sets$go.subs$BP])
CC_GO <- names(go.sets.hs[go.subs.hs$CC])
MF_GO <- names(go.sets.hs[go.subs.hs$MF])

GO_res_homo <- GO_res_homo %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(dir = ifelse(mean_change_PH_vs_BH < 0, "hypo", "hyper")) %>% 
        dplyr::select(1,2,3,4,5,6,8) %>% 
        arrange(desc(-P.DE)) %>% 
        mutate(subset = ifelse(pathway %in% BP_GO, "Biological processes", 
                               ifelse(pathway %in% CC_GO, "Cellular components", "Molecular function"))) 

write.csv(GO_res_homo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/GO_all_pathways_homogenate.csv")

GO_res_myo <- GO_res_myo %>% 
        rownames_to_column(var = "pathway") %>% 
        merge(.,pathway_dir_go, by = "pathway") %>% 
        mutate(dir = ifelse(mean_change_PM_vs_BM < 0, "hypo", "hyper")) %>% 
        dplyr::select(1,2,3,4,5,7,8) %>% 
        arrange(desc(-P.DE)) %>% 
        mutate(subset = ifelse(pathway %in% BP_GO, "Biological processes", 
                               ifelse(pathway %in% CC_GO, "Cellular components", "Molecular function"))) 


write.csv(GO_res_myo, file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/GO_all_pathways_myonuclei.csv")



##############################################################################

### GO plot

############################################################################

GO_res_homo %>% 
        arrange(desc(-P.DE)) %>% 
        head(20) -> top_go


pathway_levels = top_go %>% pull(pathway)

GO_res_myo %>% 
        filter(pathway %in% top_go$pathway) %>% 
        merge(top_go, ., by = "pathway") %>% 
        dplyr::select(pathway, PH_vs_BH = mean_change_PH_vs_BH, PM_vs_BM = mean_change_PM_vs_BM, subset = subset.x) %>% 
        pivot_longer(names_to = "contrast", values_to = "mean_change", cols = 2:3) %>% 
        mutate(pathway = factor(pathway, levels = pathway_levels)) %>% 
        ggplot(aes(y = pathway, x = mean_change, group = contrast, color = contrast))+
        geom_point(aes(shape = subset))+
        geom_vline(xintercept = 0)+
        theme_classic()+
        theme(axis.title.y = element_blank())+
        labs(title = "Topp 20 GO pathways in Homogenate")


























# check genes

anno %>% 
        filter(cpg %in% top_loadings)




plot(pca.out$x[,1:2])

plot(pca.out$x[,2:3])
text(pca.out$x[,1:2], labels = row.names(pca.out$x), pos = 4)


# run correlation between pca.out$x[,1] with M_value from M_change dataset

cor_results <- apply(M_change[1:32], 1, function(x) {
        cor(x, pca.out$x[,1])
})

# took 2 min, way faster than expected

as.data.frame(cor_results) %>% 
        mutate(r_squared = cor_results^2) %>% 
        arrange(-abs(r_squared)) %>% 
        filter(r_squared > 0.95) %>% 
        rownames_to_column(var = "cpg") %>% 
        merge(., anno, by = "cpg") %>% 
        filter(UCSC_RefGene_Name != "NA") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1),
               entrezID = mapIds(org.Hs.eg.db, keys = UCSC_RefGene_Name, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")) -> cor_data

        
# run ORA on most correlated probes

cor_kegg <- gsameth(sig.cpg = cor_data$cpg,
                                all.cpg = rownames(beta), 
                                collection = KEGG_new$kg.sets, 
                                array.type = "EPIC")

cor_kegg %>% 
        arrange(P.DE)

write.csv(cor_kegg %>% 
                  arrange(P.DE), "/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/kegg_corr_PCA1.csv")


# check which PCA1 correlated probes change significantly with time in MYO

DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        filter(cpg %in% cor_data$cpg) %>% 
        merge(., anno, by = "cpg") %>% 
        arrange(-abs(delta_M)) 

DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        filter(cpg %in% cor_data$cpg) %>% 
        merge(., anno, by = "cpg") %>% 
        arrange(-abs(delta_M))
        


# plot CPGs that were correlated with PCA1, and significant in MYO following RT


pca.out$x[,1] %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "FP") %>% 
        dplyr::select(FP, "PCA1" = 2) -> pca_df

# color based on cell population, x = PCA1, y = M_value

DMPs_PM_vs_BM %>% filter(p.value < 0.05) %>% pull(cpg) -> z

# filter cpgs that were NOT significant in myo+int

DMPs_PH_vs_BH %>% filter(p.value > 0.05) %>% pull(cpg) -> z2


# plot CPGs that are correlated with cell population, significant in MYO, but not in MYO+INT

M_change %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% cor_data$cpg) %>% 
        filter(cpg %in% z) %>% 
        filter(cpg %in% z2) %>% 
        dplyr::select(1:33) %>% 
        pivot_longer(names_to = "FP", values_to = "M_val", cols = 2:33) %>% 
        merge(., pca_df, by = "FP") %>% 
        mutate(cell_population = substr(FP, nchar(FP), nchar(FP))) %>% 
        ggplot(aes(x = PCA1, y = M_val, fill = cell_population, color = cell_population))+
        geom_point()


# extract df

M_change %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% cor_data$cpg) %>% 
        filter(cpg %in% z) %>% 
        filter(cpg %in% z2)


# get gene name of genes hypo and hyper methylation


# hypo

M_change %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% cor_data$cpg) %>% 
        filter(cpg %in% z) %>% 
        filter(cpg %in% z2) %>% 
        filter(PM_vs_BM < 0) %>% 
        merge(., anno, by ="cpg") %>% 
        distinct(UCSC_RefGene_Name, .keep_all = TRUE) %>% 
        #filter(Relation_to_Island == "Island") %>% 
        filter(UCSC_RefGene_Name %in% rownames(pooled_res3)) %>% pull(UCSC_RefGene_Name)

# get genes with multiple cpgs

M_change %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% cor_data$cpg) %>% 
        filter(cpg %in% z) %>% 
        filter(cpg %in% z2) %>% 
        filter(PM_vs_BM < 0) %>% 
        merge(., anno, by ="cpg") %>% 
        group_by(UCSC_RefGene_Name) %>% 
        filter(n() > 1) %>% pull(UCSC_RefGene_Name) %>% unique()
    
        


# hyper

M_change %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% cor_data$cpg) %>% 
        filter(cpg %in% z) %>% 
        filter(cpg %in% z2) %>% 
        filter(PM_vs_BM > 0) %>% 
        merge(., anno, by ="cpg") %>% 
        distinct(UCSC_RefGene_Name, .keep_all = TRUE)%>% 
        filter(UCSC_RefGene_Name %in% rownames(pooled_res3)) %>% pull(UCSC_RefGene_Name)



# get genes with multiple cpgs

M_change %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% cor_data$cpg) %>% 
        filter(cpg %in% z) %>% 
        filter(cpg %in% z2) %>% 
        filter(PM_vs_BM > 0) %>% 
        merge(., anno, by ="cpg") %>% 
        group_by(UCSC_RefGene_Name) %>% 
        filter(n() > 1) %>% pull(UCSC_RefGene_Name) %>% unique()        # zero


#### add KEGG pathway to genes

M_change %>% 
        rownames_to_column(var = "cpg") %>% 
        filter(cpg %in% cor_data$cpg) %>% 
        filter(cpg %in% z) %>% 
        filter(cpg %in% z2) %>% 
        filter(PM_vs_BM < 0) %>% 
        merge(., anno, by ="cpg")


DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        filter(cpg %in% cor_data$cpg) %>% 
        filter(cpg %in% z) %>% 
        filter(cpg %in% z2) %>% 
        merge(., anno, by ="cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1),
               entrezID = mapIds(org.Hs.eg.db, keys = UCSC_RefGene_Name, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")) %>% 
        mutate(KEGG_pathway = imap(KEGG_new$kg.sets, ~if(entrezID %in% .x) .y else NULL))

DMPs_PM_vs_BM %>%
        filter(p.value < 0.05) %>%
        filter(cpg %in% cor_data$cpg) %>%
        filter(cpg %in% z) %>%
        filter(cpg %in% z2) %>%
        merge(., anno, by ="cpg") %>%
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1),
               entrezID = mapIds(org.Hs.eg.db, keys = UCSC_RefGene_Name, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")) %>%
        rowwise() %>%  # Important: This ensures that each row is treated separately
        mutate(KEGG_pathway = paste(
                names(KEGG_new$kg.sets)[sapply(KEGG_new$kg.sets, function(x) entrezID %in% x)],
                collapse = "; "
        )) %>%
        ungroup() -> PCA_MYO_GENES


write.csv(PCA_MYO_GENES, "/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/MYO_only__genes.csv")


##############################################################################

### Cell specific gene lists - differentially methylated

################################################################

# load cell specific gene lists from https://cells.ucsc.edu/?ds=muscle-cell-atlas

files <- list.files("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Cell_specific_genes/", full.names = TRUE)


for (i in 1:length(files)) {
        gunzip(filename = files[i], remove = TRUE)
        print(i)
}


tsv_files <- list.files("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Cell_specific_genes/", full.names = TRUE)


for (i in 1:length(tsv_files)) {
        z <- sapply(strsplit(tsv_files[i], split =  "/"), "[", 10)
        
        z <- str_remove(z, ".tsv")
         
        x <- read_tsv(file = tsv_files[i])
        
        assign(z, x)
        print(i)
}


# check if myonuclei genes significant in MYO or MYO+INT

DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% Myonuclei$symbol) %>% 
        distinct(cpg, .keep_all = TRUE) %>% 
        mutate(symbol = UCSC_RefGene_Name) %>% 
        merge(., Myonuclei, by = "symbol") %>% 
        dplyr::select(1:10, 13) -> MYO_myonuclei

write.csv(MYO_myonuclei, "./MYO_myonuclei.csv")

DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% Myonuclei$symbol) %>% 
        distinct(cpg, .keep_all = TRUE) %>% 
        mutate(symbol = UCSC_RefGene_Name) %>% 
        merge(., Myonuclei, by = "symbol") %>% 
        dplyr::select(1:10, 13) -> MYOINT_myonuclei

write.csv(MYOINT_myonuclei, "./MYOINT_myonuclei.csv")


# check satellite cell genes in MYO and MYO+INT


DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% MuSCsandprogenitors1$symbol) %>% 
        distinct(cpg, .keep_all = TRUE) %>% 
        mutate(symbol = UCSC_RefGene_Name) %>% 
        merge(., MuSCsandprogenitors1, by = "symbol") %>% 
        dplyr::select(1:10, 13)-> MYO_MuSC_Progenitor1

write.csv(MYO_MuSC_Progenitor1, "./MYO_MuSC_Progenitor1.csv")

DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% MuSCsandprogenitors1$symbol) %>% 
        distinct(cpg, .keep_all = TRUE) %>% 
        mutate(symbol = UCSC_RefGene_Name) %>% 
        merge(., MuSCsandprogenitors1, by = "symbol") %>% 
        dplyr::select(1:10, 13)-> MYOINT_MuSC_Progenitor1

write.csv(MYOINT_MuSC_Progenitor1, "./MYOINT_MuSC_Progenitor1.csv")


DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% MuSCsandprogenitors2$symbol) %>% 
        distinct(cpg, .keep_all = TRUE)%>% 
        mutate(symbol = UCSC_RefGene_Name) %>% 
        merge(., MuSCsandprogenitors2, by = "symbol") %>% 
        dplyr::select(1:10, 13) -> MYO_MuSC_Progenitor2

write.csv(MYO_MuSC_Progenitor2, "./MYO_MuSC_Progenitor2.csv")

DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% MuSCsandprogenitors2$symbol) %>% 
        distinct(cpg, .keep_all = TRUE)%>% 
        mutate(symbol = UCSC_RefGene_Name) %>% 
        merge(., MuSCsandprogenitors2, by = "symbol") %>% 
        dplyr::select(1:10, 13) -> MYOINT_MuSC_Progenitor2

write.csv(MYOINT_MuSC_Progenitor2, "./MYOINT_MuSC_Progenitor2.csv")



# check innflamatory machrofages in MYO and MYO+INT

DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% Inflammatorymacrophagesandmonocytes$symbol) %>% 
        distinct(cpg, .keep_all = TRUE) %>% 
        mutate(symbol = UCSC_RefGene_Name) %>% 
        merge(., Inflammatorymacrophagesandmonocytes, by = "symbol") %>% 
        dplyr::select(1:10, 13)-> MYO_Inf_macrophages

write.csv(MYO_Inf_macrophages, "./MYO_Inf_macrophages.csv")

DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% Inflammatorymacrophagesandmonocytes$symbol) %>% 
        distinct(cpg, .keep_all = TRUE) %>% 
        mutate(symbol = UCSC_RefGene_Name) %>% 
        merge(., Inflammatorymacrophagesandmonocytes, by = "symbol") %>% 
        dplyr::select(1:10, 13)-> MYOINT_Inf_macrophages

write.csv(MYOINT_Inf_macrophages, "./MYOINT_Inf_macrophages.csv")


# check fibroblasts in MYO and MYO + INT

# pool Fibroblast genes

Fibroblasts1 <- Fibroblasts1 %>% mutate(cell_population = "Fibroblast1")
Fibroblasts2 <- Fibroblasts2 %>% mutate(cell_population = "Fibroblast2")
Fibroblasts3 <- Fibroblasts3 %>% mutate(cell_population = "Fibroblast3")

rbind(Fibroblasts1, Fibroblasts2) %>% 
        rbind(., Fibroblasts3) %>% 
        distinct(symbol, .keep_all = TRUE)-> Fibroblasts



DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% Fibroblasts$symbol) %>% 
        distinct(cpg, .keep_all = TRUE) %>% 
        mutate(symbol = UCSC_RefGene_Name) %>% 
        merge(., Fibroblasts, by = "symbol") %>% 
        dplyr::select(1:10, 13, 16)-> MYO_fibroblasts

write.csv(MYO_fibroblasts, "./MYO_fibroblasts.csv")

DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% Fibroblasts$symbol) %>% 
        distinct(cpg, .keep_all = TRUE) %>% 
        mutate(symbol = UCSC_RefGene_Name) %>% 
        merge(., Fibroblasts, by = "symbol") %>% 
        dplyr::select(1:10, 13, 16)-> MYOINT_fibroblasts

write.csv(MYOINT_fibroblasts, "./MYOINT_fibroblasts.csv")


      
# check overlap with pooled transcriptome

MYO_myonuclei %>% 
        filter(symbol %in% rownames(pooled_res3))

MYOINT_myonuclei %>% 
        filter(symbol %in% rownames(pooled_res3))

MYO_MuSC_Progenitor1 %>% 
        filter(symbol %in% rownames(pooled_res3))

MYOINT_MuSC_Progenitor1 %>% 
        filter(symbol %in% rownames(pooled_res3))

MYO_MuSC_Progenitor2 %>% 
        filter(symbol %in% rownames(pooled_res3))

MYOINT_MuSC_Progenitor2 %>% 
        filter(symbol %in% rownames(pooled_res3))



##############################################################################

### Muscle fiber type specific genes, supplementary file 4 in Rubenstein et al. 2020

################################################################

# https://doi.org/10.1038/s41598-019-57110-6


Rubenstein_type1 <- c("ANKRD2",
                      "ATP2A2",
                      "CA3",
                      "CASQ2",
                      "CD36",
                      "CYB5R1",
                      "FABP3",
                      "LDHB",
                      "MYH7",
                      "MYL12A",
                      "MYL2",
                      "MYL3",
                      "MYL6B",
                      "MYOZ2",
                      "PDLIM1",
                      "PLN",
                      "TNNC1",
                      "TNNI1",
                      "TNNT1",
                      "TPM3"
)

Rubenstein_type2a <- c("ALDOA",
                       "ATP2A1",
                       "DDIT4L",
                       "ENO3",
                       "G0S2",
                       "GAPDH",
                       "LDHA",
                       "MYBPC2",
                       "MYH1",
                       "MYH2",
                       "MYL1",
                       "MYLPF",
                       "PFKM",
                       "PGM1",
                       "PKM",
                       "SLN",
                       "TNNC2",
                       "TNNI2",
                       "TNNT3",
                       "TPM1"
)


# check if type1 or 2 genes are significant in MYO or MYO+INT




DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% Rubenstein_type1) 

DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% Rubenstein_type2a) 


DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% Rubenstein_type1) 

DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        filter(UCSC_RefGene_Name %in% Rubenstein_type2a) 


# other myonuclei genes

Myonuclei %>% 
        mutate(overlap = ifelse(symbol %in% Rubenstein_type1, "Type1", ifelse(symbol %in% Rubenstein_type2a, "Type2a", "NA"))) %>% 
        filter(overlap == "NA") %>% 
        dplyr::select("UCSC_RefGene_Name" = symbol, "hprd_class" = 4) -> x


DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        merge(., x, by = "UCSC_RefGene_Name")


DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% head()
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1)) %>% 
        merge(., x, by = "UCSC_RefGene_Name")





##############################################################################
        
### get DMPs with lmFit and eBayes
        
################################################################

library(limma)

# make the design matrix


FP = rep(c("FP1","FP2","FP4","FP5","FP6","FP7","FP8","FP12"), times = 4)

condition = dfh$condition

design <- model.matrix(~FP + condition)

colnames(design) <- c("Intercept","FPFP12","FPFP2","FPFP4","FPFP5","FPFP6","FPFP7","FPFP8","conditionBM","conditionPH","conditionPM")


# fit the model


M_change %>%
        dplyr::select(1:32) %>% 
        lmFit(., design) -> fit


colnames(fit$coefficients)

# design contrast matrix



cont.matrix <- makeContrasts(
        Baseline_DIFF = conditionBM - Intercept,
        Post_DIFF = conditionPM - conditionPH,
        Change_MYO = conditionPM - conditionBM,
        "Change_MYO+INT" = conditionPH - Intercept,
        Interaction = (conditionPM - conditionBM) - (conditionPH - Intercept),
        levels=design
)


# fit the contrasts

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

topTable(fit2, coef="Baseline_DIFF", number = Inf) %>% filter(adj.P.Val < 0.05) %>%  nrow()
topTable(fit2, coef="Post_DIFF", number = Inf) %>% filter(adj.P.Val < 0.05) %>%  nrow()
topTable(fit2, coef="Change_MYO", number = Inf) %>% filter(P.Value < 0.05) %>%  nrow()
topTable(fit2, coef="Change_MYO+INT", number = Inf) %>% filter(adj.P.Val < 0.05) %>%  nrow()
topTable(fit2, coef="Interaction", number = Inf) %>% filter(adj.P.Val < 0.05) %>%  nrow()



# run the ebayes on the MYO samples alone


FP = rep(c("FP1","FP2","FP4","FP5","FP6","FP7","FP8","FP12"), times = 2)

condition = dfh %>% 
        filter(condition == "BM" | condition == "PM") %>% 
        pull(condition) %>% 
        factor(., levels = c("BM", "PM"))

design <- model.matrix(~FP + condition)

colnames(design) <- c("Intercept","FPFP12","FPFP2","FPFP4","FPFP5","FPFP6","FPFP7","FPFP8","conditionPM")


# fit the model


M_change %>%
        dplyr::select(9:16,25:32) %>% 
        lmFit(., design) -> fit


colnames(fit$coefficients)

# design contrast matrix



cont.matrix <- makeContrasts(
        Change_MYO = conditionPM - Intercept,
         levels=design
)


# fit the contrasts

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

topTable(fit2, coef="Change_MYO", number = Inf) %>% filter(adj.P.Val < 0.05) %>%  nrow()


plotSA(fit)


# check var in all probes

row_variance <- apply(M_change[,1:32], 1, var)

row_variance %>% 
        as.data.frame() %>% 
        dplyr::select("var" = 1) %>% 
        ggplot(aes(x = var))+
        geom_density()+
        scale_x_continuous(trans = "sqrt", n.breaks = 20)

row_variance %>% 
        as.data.frame() %>% 
        dplyr::select("var" = 1) %>%
        filter(var > 0.1) %>% nrow()



##################################################################################

### Intra-individual within timepoint correlation between methylation M value

###################################################################################


# extract each particiant into their own dataframe

# FP 1

M_change %>% 
        dplyr::select("1BH", "1BM") %>% 
        head(10000) %>% 
        ggplot(aes(x = `1BH`, y = `1BM`))+
        geom_point(size = 0.1)+
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        labs(title = "Correlation between 1BH and 1BM",
             x = "FP 1 MYO+INT",
             y = "FP 1 MYO")  

M_change %>% 
        dplyr::select("1PH", "1PM") %>% 
        head(10000) %>% 
        ggplot(aes(x = `1PH`, y = `1PM`))+
        geom_point(size = 0.1)+
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        labs(title = "Correlation between 1PH and 1PM",
             x = "FP 1 MYO+INT",
             y = "FP 1 MYO") 


# FP 2

M_change %>% 
        dplyr::select("2BH", "2BM") %>% 
        head(10000) %>% 
        ggplot(aes(x = `2BH`, y = `2BM`))+
        geom_point(size = 0.1)+
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        labs(title = "Correlation between 2BH and 2BM",
             x = "FP 2 MYO+INT",
             y = "FP 2 MYO")  

M_change %>% 
        dplyr::select("2PH", "2PM") %>% 
        head(10000) %>% 
        ggplot(aes(x = `2PH`, y = `2PM`))+
        geom_point(size = 0.1)+
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        labs(title = "Correlation between 2PH and 2PM",
             x = "FP 2 MYO+INT",
             y = "FP 2 MYO") 


# plot DMPs changed with exercise

# list of DMPs in MYO and MYO+INT

MYO <- DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        pull(cpg)

MYOINT <- DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        pull(cpg)

# combine lists, and keep only unique rownames
 
all_sig <- unique(c(MYO, MYOINT))

shared <- c(MYO, MYOINT)

shared <- shared[duplicated(shared) ]

DMPs_PH_vs_BH %>% 
        filter(cpg %in% all_sig) %>% 
        mutate(Group = ifelse(cpg %in% MYO, "MYO", ifelse(cpg %in% MYOINT, "MYO+INT", ifelse(cpg %in% shared, "Shared", "error")))) %>% 
        dplyr::select(cpg, delta_M_MYOINT = delta_M, Group) %>% 
        merge(.,DMPs_PM_vs_BM %>% 
        filter(cpg %in% all_sig) %>% 
        mutate(Group = ifelse(cpg %in% MYO, "MYO", ifelse(cpg %in% MYOINT, "MYO+INT", ifelse(cpg %in% shared, "Shared", "error")))) %>% 
        dplyr::select(cpg, delta_M_MYO = delta_M, Group)) %>% 
        mutate(Group = ifelse(cpg %in% shared, "Shared", Group),
               Size = abs(delta_M_MYO)+abs(delta_M_MYOINT)) %>% 
        ggplot(aes(x = delta_M_MYOINT, y = delta_M_MYO, color = Group, size = Size))+
        geom_point()+
        scale_color_manual(values = c("MYO" = "green", "MYO+INT" = "purple", "Shared" = "black"))+
        scale_size(range = c(0.01, 3))

# combined figure is not easy t interpret.
# split into MYO significant, MYO+INT significant and shared on both plots

DMPs_PH_vs_BH %>% 
        filter(cpg %in% MYOINT) %>% 
        mutate(Group = ifelse(cpg %in% MYO, "MYO", ifelse(cpg %in% MYOINT, "MYO+INT", ifelse(cpg %in% shared, "Shared", "error")))) %>% 
        dplyr::select(cpg, delta_M_MYOINT = delta_M, Group) %>% 
        merge(.,DMPs_PM_vs_BM %>% 
                      filter(cpg %in% MYOINT) %>% 
                      mutate(Group = ifelse(cpg %in% MYO, "MYO", ifelse(cpg %in% MYOINT, "MYO+INT", ifelse(cpg %in% shared, "Shared", "error")))) %>% 
                      dplyr::select(cpg, delta_M_MYO = delta_M, Group)) %>% 
        mutate(Group = ifelse(cpg %in% shared, "Shared", Group),
               Size = abs(delta_M_MYO)+abs(delta_M_MYOINT)) %>% 
        ggplot(aes(x = abs(delta_M_MYOINT), y = abs(delta_M_MYO), color = Group, size = Size))+
        geom_point()+
        scale_color_manual(values = c("MYO" = "green", "MYO+INT" = "purple", "Shared" = "black"))+
        scale_size(range = c(0.01, 3)) +
        geom_smooth(method = "lm")


DMPs_PM_vs_BM %>% 
        filter(cpg %in% MYO) %>% 
        mutate(Group = ifelse(cpg %in% MYO, "MYO", ifelse(cpg %in% MYOINT, "MYO+INT", ifelse(cpg %in% shared, "Shared", "error")))) %>% 
        dplyr::select(cpg, delta_M_MYOINT = delta_M, Group) %>% 
        merge(.,DMPs_PH_vs_BH %>% 
                      filter(cpg %in% MYO) %>% 
                      mutate(Group = ifelse(cpg %in% MYO, "MYO", ifelse(cpg %in% MYOINT, "MYO+INT", ifelse(cpg %in% shared, "Shared", "error")))) %>% 
                      dplyr::select(cpg, delta_M_MYO = delta_M, Group)) %>% 
        mutate(Group = ifelse(cpg %in% shared, "Shared", Group),
               Size = abs(delta_M_MYO)+abs(delta_M_MYOINT)) %>% 
        ggplot(aes(x = abs(delta_M_MYOINT), y = abs(delta_M_MYO), color = Group, size = Size))+
        geom_point()+
        scale_color_manual(values = c("MYO" = "green", "MYO+INT" = "purple", "Shared" = "black"))+
        scale_size(range = c(0.01, 3)) +
        geom_smooth(method = "lm")+
        theme_classic()+
        scale_y_continuous()



