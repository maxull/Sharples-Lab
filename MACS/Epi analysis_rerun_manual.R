############################################################################################
###                                                                                      ###
###  Redoing analysis with manual DMPs                                                   ###
###                                                                                      ###
############################################################################################

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

setwd("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Epigenetics/")

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
        summarise(count(.,Relation_to_Island)) %>% 
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
        theme_classic()+
        theme(legend.title = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_text(hjust = 0.58))+
        scale_y_continuous(n.breaks = 6)+
        labs(y = "Number of DMPs",
             title = "Baseline DMPs in Myonuclei vs. Homogenate")


# find number of significant cpgs that were altered more than <0.01 B < 0.1 B < 0.2 B <


change_b_vals_baseline %>% 
        mutate(change = abs(change)) %>% # change all numbers to positive
        merge(baseline_dmps,.,by.x = "cpg") %>% 
        # find values in ranges
        mutate(skew = ifelse(change <0.01, "<0.01", 
                             ifelse(change >0.01 & change < 0.1, "0.01<0.1", 
                                    ifelse(change > 0.1 & change < 0.2, "0.1<0.2", 
                                           ifelse(change > 0.2 & change < 0.3, "0.2<0.3",">0.3"))))) %>% 
        summarise(count(.,skew)) -> df1

df2 <- data.frame(category = factor(df1[,1], levels = c("<0.01", "0.01<0.1", "0.1<0.2", "0.2<0.3", ">0.3")),
                  count = df1[,2])

library(scales)
p3 <- ggplot(data = df2, aes(x = category))+
        geom_bar(aes(y = count), stat = "identity")+
        geom_text(aes(label = count, y = count), vjust = -0.5)+
        theme_classic()+
        labs(title = "number of DMPs by skewed B-values at baseline",
             x = "B-value skew")+
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
        summarise(count(.,Relation_to_Island)) %>% 
        as.data.frame()

p_count2 <- merge(post_dmps, change_b_vals_post, by.x = "cpg") %>% 
        filter(dir == "Positive_skew") %>% 
        mutate(Relation_to_Island = factor(Relation_to_Island, levels = c("S_Shelf","N_Shelf","Island", "S_Shore", "N_Shore", "OpenSea"))) %>% 
        summarise(count(.,Relation_to_Island)) %>% 
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
        theme_classic()+
        theme(legend.title = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_text(hjust = 0.58))+
        scale_y_continuous(labels = label_comma(),
                           n.breaks = 7)+
        labs(y = "Number of DMPs",
             title = "Post DMPs in Myonuclei vs. Homogenate")


# find number of significant cpgs that were altered more than <0.01 B < 0.1 B < 0.2 B <


change_b_vals_post %>% 
        mutate(change = abs(change)) %>% # change all numbers to positive
        merge(post_dmps,.,by.x = "cpg") %>% 
        # find values in ranges
        mutate(skew = ifelse(change <0.01, "<0.01", 
                             ifelse(change >0.01 & change < 0.1, "0.01<0.1", 
                                    ifelse(change > 0.1 & change < 0.2, "0.1<0.2", 
                                           ifelse(change > 0.2 & change < 0.3, "0.2<0.3",">0.3"))))) %>% 
        summarise(count(.,skew)) -> df3

df4 <- data.frame(category = factor(df3[,1], levels = c("<0.01", "0.01<0.1", "0.1<0.2", "0.2<0.3", ">0.3")),
                  count = df3[,2])

p4 <- ggplot(data = df4, aes(x = category))+
        geom_bar(aes(y = count), stat = "identity")+
        geom_text(aes(label = count, y = count), vjust = -0.5)+
        theme_classic()+
        labs(title = "number of DMPs by skewed B-values at post",
             x = "B-value skew")+
        theme(axis.title.y = element_blank())+
        scale_y_continuous(trans = "log2", expand  = c(0, 1), n.breaks = 9)


plot_grid(p1,p3,p2,p4, rel_widths = c(1.5,1), labels = c("A", "B", "C", "D"))



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
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island") %>% 
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


#################################################################################

### total DMPs

###################################################################################

# count all DMPs and plot

# count DMPs and plot

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
        theme_classic() +
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
        theme_classic()+
        labs(y = "MYO + INT DMPs after 7 weeks RT")+
        geom_label(aes(x = Relation_to_Island, y = -6000, label = total))+
        geom_label(aes(x = Relation_to_Island, y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), vjust = c(0.3,-0.4,0.3,1.5,-0.4,0.3), alpha = 0.6)+
        geom_label(aes(x = Relation_to_Island, y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), vjust = c(0.7,1.3,0.7,-0.6,1.3,0.7), alpha = 0.6)+
        theme(axis.title.x = element_blank())

# plot myonuclei post vs. myonuclei baseline DMPs

p3 <- ggplot(data = mDMP, aes(x = Relation_to_Island))+
        geom_bar(aes(y = hyper), stat = "identity", fill = hyper_col, width = 0.9)+
        geom_bar(aes(y = -hypo), stat = "identity", fill = hypo_col, width = 0.9)+
        theme_classic()+
        labs(y = "MYO DMPs after 7 weeks RT")+
        geom_label(aes(x = Relation_to_Island, y = -11000, label = total))+
        geom_label(aes(x = Relation_to_Island, y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), vjust = c(0.3,-0.4,0.3,1.5,-0.4,0.3), alpha = 0.6)+
        geom_label(aes(x = Relation_to_Island, y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), vjust = c(0.7,1.3,0.7,-0.6,1.3,0.7), alpha = 0.6)+
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
        theme_classic()+
        theme(axis.title.x = element_blank())+
        labs(y = "DMPs in Islands and Promoters")

DMPs <- list(
        Homogenate_hypo = DMPs_PH_vs_BH %>% 
                filter(delta_M < 0) %>% 
                pull(cpg),
        Homogenate = DMPs_PH_vs_BH %>% 
                filter(delta_M > 0) %>% 
                pull(cpg),
        Myonuclei = DMPs_PM_vs_BM %>% 
                filter(delta_M < 0) %>% 
                pull(cpg),
        Myonuclei_hyper = DMPs_PM_vs_BM %>% 
                filter(delta_M > 0) %>% 
                pull(cpg)
)

p5 <- ggvenn(DMPs, set_name_size = 10, stroke_size = 1, set_name_color = c("White", "Black", "Black","White"),
       fill_color = c("#453781FF", "#DCE319FF","#453781FF", "#DCE319FF"),text_size = 8,stroke_alpha = 0.8)

p6 <- plot_grid(p5, labels = "C", label_x = 0.1, label_size = 14)


# figure 8 in manuscript

part1 <- plot_grid(p1,p4,p6, nrow = 1, labels = c("A", "B"), rel_widths = c(0.5,0.5,2))

part2 <- plot_grid(p2,p3, nrow = 1, labels = c("D", "E"))

DMP_RT_plot <- plot_grid(part1, part2, ncol = 1)

ggsave2(DMP_RT_plot, filename = "DMP_RT_plot.pdf",units = "cm", width = 19, height = 21, bg = "white")


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
        theme_classic()+
        theme(panel.grid.major.y = element_line(color = "red", 
                                                size = 0.5,
                                                linetype = 2),
              axis.title.x = element_blank())+
        labs(title = "Homogenate DMPs after 7 weeks of RT",
             color = "delta_M",
             y = "un-adj. P value")


# get the dataframe

DMPs_PH_vs_BH %>% 
        merge(., anno, by = "cpg") %>%
        mutate(Gene = ifelse(delta_M < -0.8 | delta_M > 0.8 | p.value < 0.0001, as.character(UCSC_RefGene_Name), ""),
               Gene = ifelse(Gene == "NA", "", Gene)) %>%
        filter(Gene != "") %>% 
        dplyr::select(!7) %>% 
        arrange(desc(-p.value)) %>% 
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
        theme_classic()+
        theme(panel.grid.major.y = element_line(color = "red", 
                                                size = 0.5,
                                                linetype = 2),
              axis.title.x = element_blank())+
        labs(title = "Myonuclei DMPs after 7 weeks of RT",
             color = "delta_M",
             y = "un-adj. P value")

DMPs_PM_vs_BM %>% 
        merge(., anno, by = "cpg") %>%
        mutate(Gene = ifelse(delta_M < -0.8 | delta_M > 0.8 | p.value < 0.0001, as.character(UCSC_RefGene_Name), ""),
               Gene = ifelse(Gene == "NA", "", Gene)) %>%
        filter(Gene != "") %>% 
        dplyr::select(!7) %>% 
        arrange(desc(-p.value)) %>% 
        write.csv(file = "C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/Myonuclear_DMPs_after_RT_volcano-plot.csv")






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

unique_gene <- unique(x$UCSC_RefGene_Name) 
unique_gene <- Filter(function(x) x != "", unique_gene)

gene_cpgs <- list()

# create subset of probes annotated to individual genes

for (i in 1:length(unique_gene)) {
        y <- x %>% 
                filter(UCSC_RefGene_Name == unique_gene[i]) %>% 
                pull(Name)
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
results_df


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



kegg_res_homo <- gsameth(sig.cpg = DMPs_PH_vs_BH$cpg,
                         all.cpg = rownames(beta), 
                         collection = kegg.sets.hs, 
                         array.type = "EPIC")

kegg_res_myo <- gsameth(sig.cpg = DMPs_PM_vs_BM$cpg,
                         all.cpg = rownames(beta), 
                         collection = kegg.sets.hs, 
                         array.type = "EPIC")

# rerun gsameth with only signalling and metabolic kegg pathways

data("sigmet.idx.hs")
kegg.subset = kegg.sets.hs[sigmet.idx.hs]

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

for (i in 1:length(kegg.sets.hs)) {
        pathway = names(kegg.sets.hs[i])
        pathway_genes = kegg.sets.hs[[pathway]]
        
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

BP_GO <- names(go.sets.hs[go.subs.hs$BP])
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


#################################################################################################

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


######################################################################################################

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

pathway = "hsa04510 Focal adhesion"     # write name of pathway

# find average methylation of the individual genes

results_df %>% 
        filter(entrezIDs %in% kegg.sets.hs[[pathway]])


_________________________________________________
# GO

pathway = "GO:0048468 cell development"       # write name of pathway

# find average methylation of the individual genes

results_df %>% 
        filter(entrezIDs %in% go.sets.hs[[pathway]])


