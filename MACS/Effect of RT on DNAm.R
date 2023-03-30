################################################################
###                                                          ###
### MACS - effect of 7 weeks RT                              ###
###                                                          ###
################################################################



# files

saveRDS(bmiq_norm, file = "bmiq_norm.RDATA")    # bmiq normalized and filtered B values

saveRDS(myNorm_fun, file = "fun_norm.RDATA")    # functional normalized and filtered B values

saveRDS(myDMP_BH_PH_bmiq, file = "myDMP_BH_PH_bmiq.RDATA")   # homogenate unadjusterd p value 0.05, DMPs calculated from M values

saveRDS(myDMP_BM_PM_bmiq, file = "myDMP_BM_PM_bmiq.RDATA")   # myonuclear unadjusterd p value 0.05, DMPs calculated from M values








# count DMPs and plot

myDMP_BH_PH_bmiq$BH_to_PH %>% 
        as.data.frame() %>% 
        dplyr::select(deltaBeta) %>% 
        dplyr::summarise(hypo = sum(deltaBeta < 0),
                         hyper = sum(deltaBeta > 0)) %>% 
        mutate(comp = "BH_to_PH",
               norm = "BMIQ") -> x




myDMP_BM_PM_bmiq$BM_to_PM %>% 
        as.data.frame() %>% 
        dplyr::select(deltaBeta) %>% 
        dplyr::summarise(hypo = sum(deltaBeta < 0),
                         hyper = sum(deltaBeta > 0)) %>% 
        mutate(comp = "BM_to_PM",
               norm = "BMIQ") %>% 
        base::rbind(x,.)->x



myDMP_BH_PH$BH_to_PH %>% 
        as.data.frame() %>% 
        dplyr::select(deltaBeta) %>% 
        dplyr::summarise(hypo = sum(deltaBeta < 0),
                         hyper = sum(deltaBeta > 0)) %>% 
        mutate(comp = "BH_to_PH",
               norm = "Fun-norm")%>% 
        base::rbind(x,.)->x




myDMP_BM_PM$BM_to_PM %>% 
        as.data.frame() %>% 
        dplyr::select(deltaBeta) %>% 
        dplyr::summarise(hypo = sum(deltaBeta < 0),
                         hyper = sum(deltaBeta > 0)) %>% 
        mutate(comp = "BM_to_PM",
               norm = "Fun-norm")%>% 
        base::rbind(x,.)->x

x %>% 
        ggplot(aes(x = comp, fill = norm))+
        geom_bar(aes(y = -hypo), stat = "identity", position = "dodge", width = 0.6)+
        geom_bar(aes(y = hyper), stat = "identity", position = "dodge", width = 0.6)


# DMPs in fun-normalized vary considerably from BMIQ normalized
# Will therefore continue the analysis for time effect in fun-normalized M-values
library(ggrepel)
# color codes

hypo_col <- "#453781FF"
hyper_col <- "#DCE319FF"


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

p1 <- ggplot(data = hDMP, aes(x = Relation_to_Island))+
        geom_bar(aes(y = hyper), stat = "identity", position = "dodge", fill = hyper_col, width = 0.9)+
        geom_bar(aes(y = -hypo), stat = "identity", position = "dodge", fill = hypo_col, width = 0.9)+
        theme_classic()+
        labs(y = "MYO + INT DMPs after 7 weeks RT")+
        geom_label(aes(x = Relation_to_Island, y = -3000, label = total))+
        geom_label(aes(x = Relation_to_Island, y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), vjust = c(0.3,-0.4,0.3,1.5,-0.4,0.3), alpha = 0.6)+
        geom_label(aes(x = Relation_to_Island, y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), vjust = c(0.7,1.3,0.7,-0.6,1.3,0.7), alpha = 0.6)+
        theme(axis.title.x = element_blank())

p2 <- ggplot(data = mDMP, aes(x = Relation_to_Island))+
        geom_bar(aes(y = hyper), stat = "identity", position = "dodge", fill = hyper_col, width = 0.9)+
        geom_bar(aes(y = -hypo), stat = "identity", position = "dodge", fill = hypo_col, width = 0.9)+
        theme_classic()+
        labs(y = "MYO DMPs after 7 weeks RT")+
        geom_label(aes(x = Relation_to_Island, y = -4400, label = total))+
        geom_label(aes(x = Relation_to_Island, y = hyper, label = paste(hyper, "/" ,round(percent_hyper, 1), "%")), vjust = c(0.3,-0.4,0.3,1.5,-0.4,0.3), alpha = 0.6)+
        geom_label(aes(x = Relation_to_Island, y = -hypo, label = paste(hypo, "/" ,round(percent_hypo, 1), "%")), vjust = c(0.7,1.3,0.7,-0.6,1.3,0.7), alpha = 0.6)+
        theme(axis.title.x = element_blank())

library(cowplot)

plot_grid(p1,p2,nrow = 1)













        
