############################################################################################
###                                                                                      ###
###  Redoing analysis with manual DMPs                                                   ###
###                                                                                      ###
############################################################################################

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
        filter(.<0.05) %>% nrow()

# filter all dmp lists for p.value <=0.05


DMPs_BM_vs_BH <- DMPs_BM_vs_BH %>% 
        filter(adj.p.val < 0.05)

DMPs_PM_vs_PH <- DMPs_PM_vs_PH %>% 
        filter(adj.p.val < 0.05) 

DMPs_PM_vs_BM <- DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05)

DMPs_PH_vs_BH <- DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05)%>% nrow()


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

tally(as.factor(som_data$kmeans))


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
              axis.title.x = element_blank())+
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

# plot with heatmap (doesnt wirk)
plot_grid(plot_heat, plot)

