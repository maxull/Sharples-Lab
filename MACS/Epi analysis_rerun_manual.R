############################################################################################
###                                                                                      ###
###  Redoing analysis with manual DMPs                                                   ###
###                                                                                      ###
############################################################################################

library(cowplot)

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

library(dendextend)

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



library(ggpubr)

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


library(ggrepel)
options()
# DMPs in homogenate after RT    
library(scales)

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


















