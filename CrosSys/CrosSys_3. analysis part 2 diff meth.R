#
#               CrosSys
#
#               calculate differential methylation
#



library(methylKit)
library(tidyverse)


# Set working directory max's mac
setwd("/Users/maxullrich/OneDrive - UGent/CrosSys")

# Set working directory NIH pc
setwd("D:/OneDrive - UGent/CrosSys/")


#################################################################################
##########      load data                               #########################
#################################################################################

# normalized and filtered data
norm.filt.dat <- readRDS(file =  "./methylation_results/norm.filt.dat.rds")




##########################################################################################################
#########               separate samples into QF and ASAT and perform dmr analysis separate                                                    ####################
##########################################################################################################

sample_names <- c()

for (i in 1:length(norm.filt.dat)) {
        sample_names[i] <- norm.filt.dat[[i]]@sample.id
}


# split data into ASAT and QF
norm.filt.dat_ASAT <- new("methylRawList", norm.filt.dat[1:35])
norm.filt.dat_QF   <- new("methylRawList", norm.filt.dat[36:68])

# fix treatment after subsetting

sample_names_asat <- c()

for (i in 1:length(norm.filt.dat_ASAT)) {
        sample_names_asat[i] <- norm.filt.dat_ASAT[[i]]@sample.id
}

sample_names_qf <- c()

for (i in 1:length(norm.filt.dat_QF)) {
        sample_names_qf[i] <- norm.filt.dat_QF[[i]]@sample.id
}

treatment_qf <- c()

treatment_qf[grep('PRE-QF', sample_names_qf)] <- 0
treatment_qf[grep('POST-QF', sample_names_qf)] <- 1

treatment_asat <- c()

treatment_asat[grep('PRE-ASAT', sample_names_asat)] <- 0
treatment_asat[grep('POST-ASAT', sample_names_asat)] <- 1


# add tretment to subsets

norm.filt.dat_ASAT@treatment <- treatment_asat

norm.filt.dat_QF@treatment <- treatment_qf

##########################################################################################################
#########               UNITE ASAT and QF                                               ####################
##########################################################################################################




# Merge everything together based on common bases
minPerGroup <- 15L
meth_ASAT <- methylKit::unite(norm.filt.dat_ASAT, destrand = F, min.per.group = minPerGroup)

meth_QF <- methylKit::unite(norm.filt.dat_QF, destrand = F, min.per.group = minPerGroup)

##########################################################################################################
#########               calculate percent methylation and run pca                                               ####################
##########################################################################################################



# extract percent methylation

# filter probes for already IDd probes

united_ID_asat <- paste0(meth_ASAT$chr,".",meth_ASAT$start, ".", meth_ASAT$end)

percent_meth_ASAT <- data.frame(ID = united_ID_asat)


for (i in 1:length(sample_names_asat)) {
        x <- data.frame(ID = paste0(norm.filt.dat_ASAT[[i]]$chr,".",norm.filt.dat_ASAT[[i]]$start,".",norm.filt.dat_ASAT[[i]]$end),
                        percent_meth = norm.filt.dat_ASAT[[i]]$numCs/norm.filt.dat_ASAT[[i]]$coverage) %>% 
                filter(ID %in% united_ID_asat) 
        
        names(x) <- c("ID",sample_names_asat[i])
                

                
        
        percent_meth_ASAT <- full_join(percent_meth_ASAT, x, by = "ID")
        
        print(i)
        
}

pca.ASAT <- prcomp(t(percent_meth_ASAT[,2:36] %>% na.omit()))

pca.ASAT$x[,1:2] %>% 
        as.data.frame() %>% 
        rownames_to_column(var = 'ID') %>% 
        mutate(Sample = factor(ifelse(grepl('PRE-ASAT', ID), 'PRE-ASAT', 
                                      ifelse(grepl('POST-ASAT', ID), 'POST-ASAT', 
                                             ifelse(grepl('PRE-QF', ID), 'PRE-QF', 'POST-QF'))),
                               levels = c('PRE-ASAT','POST-ASAT','PRE-QF','POST-QF'))) %>% 
        mutate(Tissue = ifelse(grepl('ASAT', ID), 'ASAT', 'QF')) %>% 
        #filter(Tissue == "ASAT") %>% 
        ggplot(aes(x = PC1, y = PC2))+
        #geom_point(aes(size = 2+pca.out$sdev),shape = 3)+
        geom_point(aes(fill = Sample), size = 3, stroke = 2, color = 'Black', shape = 21)+
        scale_fill_manual(values = c("#440154FF","#B396B9"))+
        theme_classic(base_size = 20)+
        labs(x = 'Principal Component 1',
             y = 'Principal Component 2')+
        theme(legend.position = 'top')+
        # geom_label_repel(aes(label = ID), max.overlaps = 36, size = 3)+
        guides(fill = guide_legend(override.aes = list(size = 6)))





# extract percent methylation

# filter probes for already IDd probes

united_ID_qf <- paste0(meth_QF$chr,".",meth_QF$start, ".", meth_QF$end)

percent_meth_QF <- data.frame(ID = united_ID_qf)


for (i in 1:length(sample_names_qf)) {
        x <- data.frame(ID = paste0(norm.filt.dat_QF[[i]]$chr,".",norm.filt.dat_QF[[i]]$start,".",norm.filt.dat_QF[[i]]$end),
                        percent_meth = norm.filt.dat_QF[[i]]$numCs/norm.filt.dat_QF[[i]]$coverage) %>% 
                filter(ID %in% united_ID_asat) 
        
        names(x) <- c("ID",sample_names_qf[i])
        
        
        
        
        percent_meth_QF <- full_join(percent_meth_QF, x, by = "ID")
        
        print(i)
        
}

pca.QF <- prcomp(t(percent_meth_QF[,2:34] %>% na.omit()))

pca.QF$x[,1:2] %>% 
        as.data.frame() %>% 
        rownames_to_column(var = 'ID') %>% 
        mutate(Sample = factor(ifelse(grepl('PRE-ASAT', ID), 'PRE-ASAT', 
                                      ifelse(grepl('POST-ASAT', ID), 'POST-ASAT', 
                                             ifelse(grepl('PRE-QF', ID), 'PRE-QF', 'POST-QF'))),
                               levels = c('PRE-ASAT','POST-ASAT','PRE-QF','POST-QF'))) %>% 
        mutate(Tissue = ifelse(grepl('ASAT', ID), 'ASAT', 'QF')) %>% 
        #filter(Tissue == "ASAT") %>% 
        ggplot(aes(x = PC1, y = PC2))+
        #geom_point(aes(size = 2+pca.out$sdev),shape = 3)+
        geom_point(aes(fill = Sample), size = 3, stroke = 2, color = 'Black', shape = 21)+
        scale_fill_manual(values = c("#5DC863FF","#B3D7B1"))+
        theme_classic(base_size = 20)+
        labs(x = 'Principal Component 1',
             y = 'Principal Component 2')+
        theme(legend.position = 'top')+
        geom_label_repel(aes(label = ID), max.overlaps = 34, size = 3)+
        guides(fill = guide_legend(override.aes = list(size = 6)))
