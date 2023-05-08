################################################################
###                                                          ###
### MACS - effect of 7 weeks RT                              ###
###                                                          ###
################################################################



# files


saveRDS(myNorm_fun, file = "fun_norm.RDATA")    # functional normalized and filtered B values

saveRDS(myDMP_BH_PH, file = "myDMP_BH_PH.RDATA")   # homogenate unadjusterd p value 0.05, DMPs calculated from M values

saveRDS(myDMP_BM_PM, file = "myDMP_BM_PM.RDATA")   # myonuclear unadjusterd p value 0.05, DMPs calculated from M values

saveRDS(beta, file = "beta.RDATA")

saveRDS(M_change, file = "M_change.RDATA")




# count DMPs and plot

myDMP_BH_PH$BH_to_PH %>% 
        as.data.frame() %>% 
        dplyr::select(deltaBeta) %>% 
        dplyr::summarise(hypo = sum(deltaBeta < 0),
                         hyper = sum(deltaBeta > 0)) %>% 
        mutate(comp = "BH_to_PH") -> x




myDMP_BM_PM$BM_to_PM %>% 
        as.data.frame() %>% 
        dplyr::select(deltaBeta) %>% 
        dplyr::summarise(hypo = sum(deltaBeta < 0),
                         hyper = sum(deltaBeta > 0)) %>% 
        mutate(comp = "BM_to_PM") %>% 
        base::rbind(x,.)->x





x %>% 
        mutate(percent_hyper = hyper/(hypo+hyper)*100,
               percent_hypo = hypo/(hypo+hyper)*100) %>% 
        ggplot(aes(x = comp))+
        geom_bar(aes(y = -hypo), stat = "identity", position = "dodge", width = 0.6, fill = hypo_col)+
        geom_bar(aes(y = hyper), stat = "identity", position = "dodge", width = 0.6, fill = hyper_col)+
        geom_label(aes(x = comp, y = hyper, label = paste(hyper, "/", round(percent_hyper, 1), "%")), alpha = 0.8)+
        geom_label(aes(x = comp, y = -hypo, label = paste(hypo, "/", round(percent_hypo, 1), "%")), alpha = 0.8)+
        theme_classic() +
        theme(axis.title.x = element_blank())+
        labs(y = "DMPs after 7 weeks RT (un adj.p < 0.05" )



# DMPs in fun-normalized vary considerably from BMIQ normalized
# Will therefore continue the analysis for time effect in fun-normalized M-values
library(ggrepel)
# color codes

hypo_col <- as.character("#453781FF")
hyper_col <- as.character("#DCE319FF")


hDMP <- myDMP_BH_PH$BH_to_PH %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, BH_AVG, PH_AVG, deltaM = deltaBeta, gene) %>% 
        merge(.,anno, by = "cpg") %>% 
        group_by(Relation_to_Island) %>% 
        dplyr::summarise(hypo = sum(deltaM < 0),
                         hyper = sum(deltaM > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) 

mDMP <- myDMP_BM_PM$BM_to_PM %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, BM_AVG, PM_AVG, deltaM = deltaBeta, gene) %>% 
        merge(.,anno, by = "cpg") %>% 
        group_by(Relation_to_Island) %>% 
        dplyr::summarise(hypo = sum(deltaM < 0),
                         hyper = sum(deltaM > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100)

# plot homogenate post vs. homogenate baseline DMPs

p1 <- ggplot(data = hDMP, aes(x = Relation_to_Island))+
        geom_bar(aes(y = hyper), stat = "identity", width = 0.9, fill = hyper_col)+
        geom_bar(aes(y = -hypo), stat = "identity", width = 0.9, fill = hypo_col)+
        theme_classic()+
        labs(y = "MYO + INT DMPs after 7 weeks RT")+
        geom_label(aes(x = Relation_to_Island, y = -3000, label = total))+
        geom_label(aes(x = Relation_to_Island, y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), vjust = c(0.3,-0.4,0.3,1.5,-0.4,0.3), alpha = 0.6)+
        geom_label(aes(x = Relation_to_Island, y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), vjust = c(0.7,1.3,0.7,-0.6,1.3,0.7), alpha = 0.6)+
        theme(axis.title.x = element_blank())

# plot myonuclei post vs. myonuclei baseline DMPs

p2 <- ggplot(data = mDMP, aes(x = Relation_to_Island))+
        geom_bar(aes(y = hyper), stat = "identity", fill = hyper_col, width = 0.9)+
        geom_bar(aes(y = -hypo), stat = "identity", fill = hypo_col, width = 0.9)+
        theme_classic()+
        labs(y = "MYO DMPs after 7 weeks RT")+
        geom_label(aes(x = Relation_to_Island, y = -4400, label = total))+
        geom_label(aes(x = Relation_to_Island, y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), vjust = c(0.3,-0.4,0.3,1.5,-0.4,0.3), alpha = 0.6)+
        geom_label(aes(x = Relation_to_Island, y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), vjust = c(0.7,1.3,0.7,-0.6,1.3,0.7), alpha = 0.6)+
        theme(axis.title.x = element_blank())+
        scale_fill_identity(name = 'the fill', guide = "legend",labels = c("hypo", "hyper"),aes(y = 0, x = 7))

library(cowplot)

plot_grid(p1,p2,nrow = 1, labels = c("A","B"))


# isolate homogenate and myonuclei DMPs within promoter assosiated islands and plot together


hDMP <- myDMP_BH_PH$BH_to_PH %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, BH_AVG, PH_AVG, deltaM = deltaBeta, gene) %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island") %>% 
        dplyr::summarise(hypo = sum(deltaM < 0),
                         hyper = sum(deltaM > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) %>% 
        mutate(condition = "homo")

mDMP <- myDMP_BM_PM$BM_to_PM %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, BM_AVG, PM_AVG, deltaM = deltaBeta, gene) %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island") %>% 
        dplyr::summarise(hypo = sum(deltaM < 0),
                         hyper = sum(deltaM > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) %>% 
        mutate(condition = "myo")

rbind(hDMP, mDMP) %>% 
        ggplot(aes(x = condition))+
        geom_bar(aes(y = hyper), stat = "identity", fill = hyper_col)+
        geom_bar(aes(y = -hypo), stat = "identity", fill = hypo_col)+
        geom_label(aes(y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), alpha = 0.8)+
        geom_label(aes(y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), alpha = 0.8)
        


###

# isolate cpgs in islands and promoters and check directional overlap

# isolate data



hDMP <- myDMP_BH_PH$BH_to_PH %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, BH_AVG, PH_AVG, deltaM = deltaBeta, gene) %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island")
        
mDMP <- myDMP_BM_PM$BM_to_PM %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, BM_AVG, PM_AVG, deltaM = deltaBeta, gene) %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island")

DMPs_pro_isl <- merge(hDMP, mDMP, by = "cpg")

# hypo homo and hypo myo

DMPs_pro_isl %>% 
        filter(deltaM.x < 0 & deltaM.y < 0) %>% 
        dplyr::select(cpg, gene.x)

# hyper homo and myo

DMPs_pro_isl %>% 
        filter(deltaM.x > 0 & deltaM.y > 0) %>% 
        dplyr::select(cpg, gene.x)

# hyper homo and hypo myo

DMPs_pro_isl %>% 
        filter(deltaM.x > 0 & deltaM.y < 0) %>% 
        dplyr::select(cpg, gene.x)

# hypo homo and hyper myo

DMPs_pro_isl %>% 
        filter(deltaM.x < 0 & deltaM.y > 0) %>% 
        dplyr::select(cpg, gene.x)


# all CPGs, but with higher delta M

# isolate data

hDMP <- myDMP_BH_PH$BH_to_PH %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, BH_AVG, PH_AVG, deltaM = deltaBeta, gene) %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(deltaM >= 0.3 | deltaM <= -0.3)

mDMP <- myDMP_BM_PM$BM_to_PM %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, BM_AVG, PM_AVG, deltaM = deltaBeta, gene) %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(deltaM >= 0.3 | deltaM <= -0.3)

DMPs_0.3 <- merge(hDMP, mDMP, by = "cpg")

# hypo homo and hypo myo

DMPs_0.3 %>% 
        filter(deltaM.x < 0 & deltaM.y < 0) %>% 
        dplyr::select(cpg, gene.x)

# hyper homo and myo

DMPs_0.3 %>% 
        filter(deltaM.x > 0 & deltaM.y > 0) %>% 
        dplyr::select(cpg, gene.x)

# hyper homo and hypo myo

DMPs_0.3 %>% 
        filter(deltaM.x > 0 & deltaM.y < 0) %>% 
        dplyr::select(cpg, gene.x)

# hypo homo and hyper myo

DMPs_0.3 %>% 
        filter(deltaM.x < 0 & deltaM.y > 0) %>% 
        dplyr::select(cpg, gene.x)



# find DMPs that are only in homo and only in myo

hDMP <- myDMP_BH_PH$BH_to_PH %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, BH_AVG, PH_AVG, deltaM = deltaBeta, gene) %>% 
        merge(.,anno, by = "cpg")


mDMP <- myDMP_BM_PM$BM_to_PM %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, BM_AVG, PM_AVG, deltaM = deltaBeta, gene) %>% 
        merge(.,anno, by = "cpg")

hDMP[hDMP$cpg %in%  mDMP$cpg == FALSE,] -> x

mDMP[mDMP$cpg %in%  hDMP$cpg == FALSE,] -> y

merge(x,y, by = "cpg")


# plot venn diagram


DMPs <- list(hDMP = hDMP$cpg,
             mDMP = mDMP$cpg)


library(ggvenn)

# plot venn-diagram 

venn <- ggvenn(DMPs, set_name_size = 10, stroke_size = 1, text_size = 10)

# re-do venn with island and promoter DMPs

hDMP <- myDMP_BH_PH$BH_to_PH %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, BH_AVG, PH_AVG, deltaM = deltaBeta, gene) %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island")

mDMP <- myDMP_BM_PM$BM_to_PM %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "cpg") %>% 
        dplyr::select(cpg, P.Value, BM_AVG, PM_AVG, deltaM = deltaBeta, gene) %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island")

DMPs <- list(hDMP = hDMP$cpg,
             mDMP = mDMP$cpg)

venn <- ggvenn(DMPs, set_name_size = 10, stroke_size = 1, text_size = 10)


###

# add fantom 4 and 5 enhancers to annotation
Illumina_anno %>% 
        dplyr::select(cpg = 1, 21:22) %>% 
        merge(anno, ., by = "cpg") -> anno


# Homogenate DMPs



# count hyper and hypo methylated probes in islands and promoters, to compare with adams percentages

hDMP %>% 
        mutate(dir = ifelse(deltaM < 0, "hypo", "hyper")) %>% 
        dplyr::summarise(count(.,dir))

mDMP %>% nrow()
        mutate(dir = ifelse(deltaM < 0, "hypo", "hyper")) %>% 
        dplyr::summarise(count(.,dir))


# load adams DMP lists to overlap the cpgs
        
Adam_dmps <- read.csv("C:/Users/maxul/Downloads/DMPs.xlsx", skip = 2)





##########################################################################################

### re-do analysis on manual DMP calculations

##########################################################################################

# manual DMPs, total count

DMPs_PH_vs_BH%>% 
        as.data.frame() %>% 
        dplyr::select(delta_M) %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0)) %>% 
        mutate(comp = "BH_to_PH") -> x




DMPs_PM_vs_BM%>% 
        as.data.frame() %>% 
        dplyr::select(delta_M) %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0)) %>% 
        mutate(comp = "BM_to_PM") %>% 
        base::rbind(x,.)->x





x %>% 
        mutate(percent_hyper = hyper/(hypo+hyper)*100,
               percent_hypo = hypo/(hypo+hyper)*100) %>% 
        ggplot(aes(x = comp))+
        geom_bar(aes(y = -hypo), stat = "identity", position = "dodge", width = 0.6, fill = hypo_col)+
        geom_bar(aes(y = hyper), stat = "identity", position = "dodge", width = 0.6, fill = hyper_col)+
        geom_label(aes(x = comp, y = hyper, label = paste(hyper, "/", round(percent_hyper, 1), "%")), alpha = 0.8)+
        geom_label(aes(x = comp, y = -hypo, label = paste(hypo, "/", round(percent_hypo, 1), "%")), alpha = 0.8)+
        theme_classic() +
        theme(axis.title.x = element_blank())+
        labs(y = "DMPs after 7 weeks RT (un adj.p < 0.05" )


###__________________________________________-

# plot DMPs in different regulatory features

hypo_col <- as.character("#453781FF")
hyper_col <- as.character("#DCE319FF")


hDMP <- DMPs_PH_vs_BH %>% 
        as.data.frame() %>% 
        merge(.,anno, by = "cpg") %>% 
        group_by(Relation_to_Island) %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) 

mDMP <- DMPs_PM_vs_BM %>% 
        as.data.frame() %>% 
        merge(.,anno, by = "cpg") %>% 
        group_by(Relation_to_Island) %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) 

# plot homogenate post vs. homogenate baseline DMPs

p1 <- ggplot(data = hDMP, aes(x = Relation_to_Island))+
        geom_bar(aes(y = hyper), stat = "identity", width = 0.9, fill = hyper_col)+
        geom_bar(aes(y = -hypo), stat = "identity", width = 0.9, fill = hypo_col)+
        theme_classic()+
        labs(y = "MYO + INT DMPs after 7 weeks RT")+
        geom_label(aes(x = Relation_to_Island, y = -6000, label = total))+
        geom_label(aes(x = Relation_to_Island, y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), vjust = c(0.3,-0.4,0.3,1.5,-0.4,0.3), alpha = 0.6)+
        geom_label(aes(x = Relation_to_Island, y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), vjust = c(0.7,1.3,0.7,-0.6,1.3,0.7), alpha = 0.6)+
        theme(axis.title.x = element_blank())


# plot myonuclei post vs. myonuclei baseline DMPs

p2 <- ggplot(data = mDMP, aes(x = Relation_to_Island))+
        geom_bar(aes(y = hyper), stat = "identity", fill = hyper_col, width = 0.9)+
        geom_bar(aes(y = -hypo), stat = "identity", fill = hypo_col, width = 0.9)+
        theme_classic()+
        labs(y = "MYO DMPs after 7 weeks RT")+
        geom_label(aes(x = Relation_to_Island, y = -10500, label = total))+
        geom_label(aes(x = Relation_to_Island, y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), vjust = c(0.3,-0.4,0.3,1.5,-0.4,0.3), alpha = 0.6)+
        geom_label(aes(x = Relation_to_Island, y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), vjust = c(0.7,1.3,0.7,-0.6,1.3,0.7), alpha = 0.6)+
        theme(axis.title.x = element_blank())+
        scale_fill_identity(name = 'the fill', guide = "legend",labels = c("hypo", "hyper"),aes(y = 0, x = 7))

library(cowplot)

plot_grid(p1,p2,nrow = 1, labels = c("A","B"))




# isolate homogenate and myonuclei DMPs within promoter assosiated islands and plot together


hDMP <- DMPs_PH_vs_BH %>% 
        as.data.frame() %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island") %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) %>% 
        mutate(condition = "homo")

mDMP <- DMPs_PM_vs_BM %>% 
        as.data.frame() %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island") %>%
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) %>% 
        mutate(condition = "myo")


rbind(hDMP, mDMP) %>% 
        ggplot(aes(x = condition))+
        geom_bar(aes(y = hyper), stat = "identity", fill = hyper_col)+
        geom_bar(aes(y = -hypo), stat = "identity", fill = hypo_col)+
        geom_label(aes(y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), alpha = 0.8)+
        geom_label(aes(y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), alpha = 0.8)+
        labs(y = "DMPs", title = "DMPs in islands and promoters")



#######################################################33

### check direction of methylation

DMPs_PM_vs_BM %>% 
        as.data.frame() %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island") -> x

x = x$cpg

M_change[rownames(M_change) %in%  x,] %>% 
        dplyr::summarise(hypo = sum(PM_vs_BM <0),
                         hyper = sum(PM_vs_BM >0))


### check overlap with previous DMPs

# get adams DMPs

library(readxl)

adam_dmps_baseline <- read_excel("C:/Users/maxul/Downloads/DMPs.xlsx", sheet = "MyoN_Base vs. HomoG_Base_170K", skip = 1)


# baseline

DMPs <- list(
        champ.dmp = rownames(myDMP_BH_BM$BH_to_BM),
        manual.dmp = DMPs_BM_vs_BH$cpg,
        adam.dmp = adam_dmps_baseline$`Probeset ID`
)

ggvenn(DMPs, set_name_size = 10, stroke_size = 1, text_size = 8)+
        labs(title = "baseline DMPs")

# post

adam_dmps_post <- read_excel("C:/Users/maxul/Downloads/DMPs.xlsx", sheet = "MyoN_Post vs. HomoG_Post_165K", skip = 1)


DMPs <- list(
        champ.dmp = rownames(myDMP_PH_PM$PH_to_PM),
        manual.dmp = DMPs_PM_vs_PH$cpg,
        adam.dmp = adam_dmps_post$`Probeset ID`
)

ggvenn(DMPs, set_name_size = 10, stroke_size = 1, text_size = 8)+
        labs(title = "post DMPs")

# change homogenate

adam_dmps_homogenate <- read_excel("C:/Users/maxul/Downloads/DMPs.xlsx", sheet = "HomoG_Post v HomoG_Base_11037", skip = 1)


DMPs <- list(
        champ.dmp = rownames(myDMP_BH_PH$BH_to_PH),
        manual.dmp = DMPs_PH_vs_BH$cpg,
        adam.dmp = adam_dmps_homogenate$`Probeset ID`
)

ggvenn(DMPs, set_name_size = 10, stroke_size = 1, text_size = 8)+
        labs(title = "homogenate post vs baseline DMPs")

# islands and promoters

DMPs <- list(
        champ.dmp = myDMP_BH_PH$BH_to_PH %>% 
                rownames_to_column(var = "cpg") %>% 
                as.data.frame %>% 
                merge(.,anno, by = "cpg") %>% 
                filter(Relation_to_Island == "Island" & Regulatory_Feature_Group == "Promoter_Associated") %>% 
                pull(cpg),
        manual.dmp = DMPs_PH_vs_BH %>% 
                merge(.,anno, by = "cpg") %>% 
                filter(Relation_to_Island == "Island" & Regulatory_Feature_Group == "Promoter_Associated") %>% 
                pull(cpg),
        adam.dmp = adam_dmps_homogenate %>% 
                as.data.frame() %>% 
                filter(Relation_to_UCSC_CpG_Island == "Island" & 
                               Regulatory_Feature_Group == "Promoter_Associated") %>% 
                pull(`Probeset ID`)
)

ggvenn(DMPs, set_name_size = 10, stroke_size = 1, text_size = 8)+
        labs(title = "homogenate post vs baseline DMPs, islands in promoters")


# change myonuclei

adam_dmps_myonuclei <- read_excel("C:/Users/maxul/Downloads/DMPs.xlsx", sheet = "MyoN_Post v MyoN_Base_11324", skip = 1)


DMPs <- list(
        champ.dmp = rownames(myDMP_BM_PM$BM_to_PM),
        manual.dmp = DMPs_PM_vs_BM$cpg,
        adam.dmp = adam_dmps_myonuclei$`Probeset ID`
)

ggvenn(DMPs, set_name_size = 10, stroke_size = 1, text_size = 8)+
        labs(title = "myonuclei post vs baseline DMPs")

# islands and promoters

DMPs <- list(
        champ.dmp = myDMP_BM_PM$BM_to_PM %>% 
                rownames_to_column(var = "cpg") %>% 
                as.data.frame %>% 
                merge(.,anno, by = "cpg") %>% 
                filter(Relation_to_Island == "Island" & Regulatory_Feature_Group == "Promoter_Associated") %>% 
                pull(cpg),
        manual.dmp = DMPs_PM_vs_BM %>% 
                merge(.,anno, by = "cpg") %>% 
                filter(Relation_to_Island == "Island" & Regulatory_Feature_Group == "Promoter_Associated") %>% 
                pull(cpg),
        adam.dmp = adam_dmps_myonuclei %>% 
                as.data.frame() %>% 
                filter(Relation_to_UCSC_CpG_Island == "Island" & 
                               Regulatory_Feature_Group == "Promoter_Associated") %>% 
                pull(`Probeset ID`)
)

ggvenn(DMPs, set_name_size = 10, stroke_size = 1, text_size = 8)+
        labs(title = "myonuclei post vs baseline DMPs, islands in promoters")



# check adams DMPs % and % in islands and promoters



adam_dmps_homogenate %>% 
        as.data.frame() %>% 
        dplyr::select(2,9,11,diff_homo = 37) %>% 
        dplyr::summarise(hypo = sum(diff_homo < 0),
                         hyper = sum(diff_homo > 0)) %>% 
        mutate(comp = "BH_to_PH") -> x
        
adam_dmps_myonuclei %>% 
        as.data.frame() %>% 
        dplyr::select(2,9,11,diff_homo = 37) %>% 
        dplyr::summarise(hypo = sum(diff_homo < 0),
                         hyper = sum(diff_homo > 0)) %>% 
        mutate(comp = "BM_to_PM") %>% 
        base::rbind(x,.)->x


x %>% 
        mutate(percent_hyper = hyper/(hypo+hyper)*100,
               percent_hypo = hypo/(hypo+hyper)*100) %>% 
        ggplot(aes(x = comp))+
        geom_bar(aes(y = -hypo), stat = "identity", position = "dodge", width = 0.6, fill = hypo_col)+
        geom_bar(aes(y = hyper), stat = "identity", position = "dodge", width = 0.6, fill = hyper_col)+
        geom_label(aes(x = comp, y = hyper, label = paste(hyper, "/", round(percent_hyper, 1), "%")), alpha = 0.8)+
        geom_label(aes(x = comp, y = -hypo, label = paste(hypo, "/", round(percent_hypo, 1), "%")), alpha = 0.8)+
        theme_classic() +
        theme(axis.title.x = element_blank())+
        labs(y = "DMPs after 7 weeks RT (un adj.p < 0.05" )


# repeat with only island and promoter DMPs

adam_dmps_homogenate %>% 
        as.data.frame() %>% 
        dplyr::select(2,9,11,diff_homo = 37) %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_UCSC_CpG_Island == "Island") %>%
        dplyr::summarise(hypo = sum(diff_homo < 0),
                         hyper = sum(diff_homo > 0)) %>% 
        mutate(comp = "BH_to_PH") -> x

adam_dmps_myonuclei %>% 
        as.data.frame() %>%
        dplyr::select(2,9,11,diff_myo = 45) %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_UCSC_CpG_Island == "Island") %>%
        dplyr::summarise(hypo = sum(diff_myo < 0),
                         hyper = sum(diff_myo > 0)) %>% 
        mutate(comp = "BM_to_PM") %>% 
        base::rbind(x,.)->x


x %>% 
        mutate(percent_hyper = hyper/(hypo+hyper)*100,
               percent_hypo = hypo/(hypo+hyper)*100) %>% 
        ggplot(aes(x = comp))+
        geom_bar(aes(y = -hypo), stat = "identity", position = "dodge", width = 0.6, fill = hypo_col)+
        geom_bar(aes(y = hyper), stat = "identity", position = "dodge", width = 0.6, fill = hyper_col)+
        geom_label(aes(x = comp, y = hyper, label = paste(hyper, "/", round(percent_hyper, 1), "%")), alpha = 0.8)+
        geom_label(aes(x = comp, y = -hypo, label = paste(hypo, "/", round(percent_hypo, 1), "%")), alpha = 0.8)+
        theme_classic() +
        theme(axis.title.x = element_blank())+
        labs(y = "DMPs after 7 weeks RT (un adj.p < 0.05" )


# check that adams island and promoter DMPs are same as mine in myonuclei


adam_isl_pro <- adam_dmps_myonuclei %>% 
        as.data.frame() %>%
        dplyr::select(cpg = 2,9,11,diff_myo = 45) %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_UCSC_CpG_Island == "Island") %>% 
        pull(cpg)-> adam


merge(DMPs_PM_vs_BM, anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated" & Relation_to_Island == "Island") %>% 
        pull(cpg)-> max
        merge(adam_isl_pro, ., by = "cpg")
dmps <- list(adam = adam, max = max)

library(ggvenn)
ggvenn(dmps, set_name_size = 10, stroke_size = 1, text_size = 8)+
        labs(title = "Islands in promoters")


# plot mine and adams DMPs filterd for promoter associated

adam_dmps_homogenate %>% 
        as.data.frame() %>% 
        dplyr::select(2,9,11,diff_homo = 37) %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated") %>%
        dplyr::summarise(hypo = sum(diff_homo < 0),
                         hyper = sum(diff_homo > 0)) %>% 
        mutate(comp = "BH_to_PH") -> x

adam_dmps_myonuclei %>% 
        as.data.frame() %>%
        dplyr::select(2,9,11,diff_myo = 45) %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated") %>%
        dplyr::summarise(hypo = sum(diff_myo < 0),
                         hyper = sum(diff_myo > 0)) %>% 
        mutate(comp = "BM_to_PM") %>% 
        base::rbind(x,.)->x


p1 <- x %>% 
        mutate(percent_hyper = hyper/(hypo+hyper)*100,
               percent_hypo = hypo/(hypo+hyper)*100) %>% 
        ggplot(aes(x = comp))+
        geom_bar(aes(y = -hypo), stat = "identity", position = "dodge", width = 0.6, fill = hypo_col)+
        geom_bar(aes(y = hyper), stat = "identity", position = "dodge", width = 0.6, fill = hyper_col)+
        geom_label(aes(x = comp, y = hyper, label = paste(hyper, "/", round(percent_hyper, 1), "%")), alpha = 0.8)+
        geom_label(aes(x = comp, y = -hypo, label = paste(hypo, "/", round(percent_hypo, 1), "%")), alpha = 0.8)+
        theme_classic() +
        theme(axis.title.x = element_blank())+
        labs(y = "DMPs after 7 weeks RT (un adj.p < 0.05",
             title = "ADAM: DMPs in promoters")

hDMP <- DMPs_PH_vs_BH %>% 
        as.data.frame() %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated") %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) %>% 
        mutate(condition = "BH_to_PH")

mDMP <- DMPs_PM_vs_BM %>% 
        as.data.frame() %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Regulatory_Feature_Group == "Promoter_Associated") %>%
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) %>% 
        mutate(condition = "BM_to_PM")


p2 <- rbind(hDMP, mDMP) %>% 
        ggplot(aes(x = condition))+
        geom_bar(aes(y = hyper), stat = "identity", fill = hyper_col, width = 0.6)+
        geom_bar(aes(y = -hypo), stat = "identity", fill = hypo_col, width = 0.6)+
        geom_label(aes(y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), alpha = 0.8)+
        geom_label(aes(y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), alpha = 0.8)+
        labs(title = "MAX: DMPs in promoters")+
        theme_classic()+
        theme(axis.title.x = element_blank())+
        labs(y = "DMPs after 7 weeks RT (un adj.p < 0.05" )

plot_grid(p1,p2, nrow = 1)


# repeat with only islands

adam_dmps_homogenate %>% 
        as.data.frame() %>% 
        dplyr::select(2,9,11,diff_homo = 37) %>% 
        filter(Relation_to_UCSC_CpG_Island == "Island") %>%
        dplyr::summarise(hypo = sum(diff_homo < 0),
                         hyper = sum(diff_homo > 0)) %>% 
        mutate(comp = "BH_to_PH") -> x

adam_dmps_myonuclei %>% 
        as.data.frame() %>%
        dplyr::select(2,9,11,diff_myo = 45) %>% 
        filter(Relation_to_UCSC_CpG_Island == "Island") %>%
        dplyr::summarise(hypo = sum(diff_myo < 0),
                         hyper = sum(diff_myo > 0)) %>% 
        mutate(comp = "BM_to_PM") %>% 
        base::rbind(x,.)->x


p1 <- x %>% 
        mutate(percent_hyper = hyper/(hypo+hyper)*100,
               percent_hypo = hypo/(hypo+hyper)*100) %>% 
        ggplot(aes(x = comp))+
        geom_bar(aes(y = -hypo), stat = "identity", position = "dodge", width = 0.6, fill = hypo_col)+
        geom_bar(aes(y = hyper), stat = "identity", position = "dodge", width = 0.6, fill = hyper_col)+
        geom_label(aes(x = comp, y = hyper, label = paste(hyper, "/", round(percent_hyper, 1), "%")), alpha = 0.8)+
        geom_label(aes(x = comp, y = -hypo, label = paste(hypo, "/", round(percent_hypo, 1), "%")), alpha = 0.8)+
        theme_classic() +
        theme(axis.title.x = element_blank())+
        labs(y = "DMPs after 7 weeks RT (un adj.p < 0.05",
             title = "ADAM: DMPs in Islands")

hDMP <- DMPs_PH_vs_BH %>% 
        as.data.frame() %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Relation_to_Island == "Island") %>% 
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) %>% 
        mutate(condition = "BH_to_PH")

mDMP <- DMPs_PM_vs_BM %>% 
        as.data.frame() %>% 
        merge(.,anno, by = "cpg") %>% 
        filter(Relation_to_Island == "Island") %>%
        dplyr::summarise(hypo = sum(delta_M < 0),
                         hyper = sum(delta_M > 0),
                         total = hyper+hypo,
                         percent_hyper = (hyper/total)*100,
                         percent_hypo = (hypo/total)*100) %>% 
        mutate(condition = "BM_to_PM")


p2 <- rbind(hDMP, mDMP) %>% 
        ggplot(aes(x = condition))+
        geom_bar(aes(y = hyper), stat = "identity", fill = hyper_col, width = 0.6)+
        geom_bar(aes(y = -hypo), stat = "identity", fill = hypo_col, width = 0.6)+
        geom_label(aes(y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), alpha = 0.8)+
        geom_label(aes(y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), alpha = 0.8)+
        labs(title = "MAX: DMPs in Islands")+
        theme_classic()+
        theme(axis.title.x = element_blank())+
        labs(y = "DMPs after 7 weeks RT (un adj.p < 0.05" )

plot_grid(p1,p2, nrow = 1)


# read adams new paired t-terst data

adams_dmps_t_test <- read_tsv(col_names = TRUE, file = "C:/Users/maxul/Downloads/Myonuclei Post vs. Baseline_Paired T test_v2 Baseline denominator+ unadj p 0.05.txt")

dmps <- list(
        adams_dmps = adams_dmps_t_test %>% filter(Relation_to_UCSC_CpG_Island == "Island" & Regulatory_Feature_Group == "Promoter_Associated") %>% pull('Probeset ID'),
        max_dmps = DMPs_PM_vs_BM %>% merge(., anno, by = "cpg") %>% filter(Relation_to_Island == "Island" & Regulatory_Feature_Group == "Promoter_Associated") %>% pull(cpg)
)
library(ggvenn)
ggvenn(dmps, set_name_size = 10, stroke_size = 1, text_size = 8)