#
#               CrosSys
#
#               combine data, PCA, normalization etc. 
#



library(methylKit)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# Set working directory max's mac
setwd("/Users/maxullrich/OneDrive - UGent/CrosSys")

# Set working directory NIH pc
setwd("D:/OneDrive - UGent/CrosSys/")

# increase data loading time
options(timeout = 600)

##########################################################################################################
#########               get all files                          ###########################################
##########################################################################################################


# Locate all files
all_meth_files <- list.files(pattern = '*.RDATA', 
                                   path = "./methylation_results/merged_meth_files", full.names = TRUE)

# names

sample_names <- gsub(pattern = '.RDATA', '',list.files(pattern = '*.RDATA', 
                                     path = "./methylation_results/merged_meth_files", full.names = FALSE))


# load data into list

samples <- c()

for (i in 1:length(all_meth_files)) {
        samples[[i]] <- readRDS(all_meth_files[i])
        print(i)
}

names(samples) <- sample_names


##########################################################################################################
#########               Check that numCs + numTs = coverage       ########################################
##########################################################################################################

# load data and skip
samples <- readRDS('methylation_results/samples.RDATA')
# 
# for (i in 1:length(samples)) {
#         
#         samples[[i]] <- samples[[i]] %>% 
#                 mutate(check = (numCs+numTs == coverage)) %>% 
#                 filter(check != FALSE) %>% 
#                 mutate(percent_meth = numCs/coverage) %>% 
#                 dplyr::select(chr, start, end, strand, coverage, numCs, numTs, percent_meth)
#         print(i)
# }
# 
# # save sample data
# saveRDS(samples, "./methylation_results/samples.RDATA")


##########################################################################################################
#########               plot counts, positions and percent meth for each sample       ####################
##########################################################################################################


qc <- tibble()

for (i in 1:length(samples)) {
        
        qc <- rbind(qc,samples[[i]] %>% 
                summarise('total_counts' = sum(coverage),
                          'positions' = length(coverage),
                          'percent_meth' = mean(percent_meth)))
}

qc$ID <- sample_names

ggplot(data = qc, aes(x = ID, y = total_counts))+
        geom_histogram(stat = 'identity') +
        scale_y_continuous(expand = c(0,0))+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90))+
        labs(y = 'Total counts after filtering')

ggplot(data = qc, aes(x = ID, y = positions))+
        geom_histogram(stat = 'identity') +
        scale_y_continuous(expand = c(0,0))+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90))+
        labs(y = 'Unique genomic position after filtering')

ggplot(data = qc, aes(x = ID, y = percent_meth*100))+
        geom_histogram(stat = 'identity') +
        scale_y_continuous(expand = c(0,0))+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90))+
        labs(y = 'Percent Methylation')


##########################################################################################################
#########               plot methylation distribution                                 ####################
##########################################################################################################

# plot muscle and fat separately

density_df <- tibble()

for (i in 1:length(samples)) {
        
        sample <- samples[[i]]
        
        sample$ID = sample_names[i]
        sample$Timepoint <- ifelse(grepl('PRE', sample_names[i]), 'Pre', 'Post' )
        sample$Tissue <- ifelse(grepl('ASAT', sample_names[i]), 'ASAT', 'QF')
        
        
        density_df <- rbind(density_df,sample %>% dplyr::select(percent_meth, ID, Timepoint, Tissue))
        print(i)
}

head(density_df)

density_df %>% 
        ggplot(aes(x = percent_meth, group = ID, color = Tissue))+
        geom_density()



##########################################################################################################
#########               normalize                                                     ####################
##########################################################################################################

# load data and skip
norm.filt.dat <- readRDS('./methylation_results/norm.filt.dat.RDATA')

norm.filt.dat <- read_rds("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/CrosSys/methylation_results/norm.filt.dat.RDATA")
# Check if the file exists
file.exists('./methylation_results/norm.filt.dat.RDATA')

getwd()
# 
# # Combine into a methylRawList
# # Convert data frames to methylRaw objects
# methylRawList_obj <- lapply(seq_along(samples), function(i) {
#         new("methylRaw",
#                 samples[[i]],
#                 sample.id = sample_names[i],
#                 assembly = "hg38",
#                 context = "CpG",
#                 resolution = "base"
#         )
# })
# 
# 
# # Combine into a methylRawList
# methylRawList_obj <- new("methylRawList", methylRawList_obj)
# 
# # Normalise coverage
# norm.filt.dat <- normalizeCoverage(methylRawList_obj)
# 
# # save normalized filetered data
# saveRDS(norm.filt.dat, "./methylation_results/norm.filt.dat.RDATA")

# combine for density plot
norm_density_df <- tibble()

for (i in 1:length(norm.filt.dat)) {
        
        sample <- as.data.frame(norm.filt.dat[[i]]@.Data)
        colnames(sample) <- c('chr','start','end','strand','coverage','numCs','numTs','percent_meth')
        
        sample$ID = sample_names[i]
        sample$Timepoint <- ifelse(grepl('PRE', sample_names[i]), 'Pre', 'Post' )
        sample$Tissue <- ifelse(grepl('ASAT', sample_names[i]), 'ASAT', 'QF')
        
        
        norm_density_df <- rbind(norm_density_df,sample %>% dplyr::select(percent_meth, ID, Timepoint, Tissue))
        print(i)
}

head(density_df)

norm_density_df %>% 
        ggplot(aes(x = percent_meth, group = ID, color = Tissue))+
        geom_density() +
        theme_classic()


##########################################################################################################
#########               combine and cluster                                           ####################
##########################################################################################################

# add treatment/condition to the normalized data

treatment <- c()

treatment[grep('PRE-ASAT', sample_names)] <- 0
treatment[grep('POST-ASAT', sample_names)] <- 1
treatment[grep('PRE-QF', sample_names)] <- 2
treatment[grep('POST-QF', sample_names)] <- 3

table(treatment)

norm.filt.dat@treatment <- treatment

# Merge everything together based on common bases
# Only merge CpGs that are in at least 6 samples per group after discussion with Adam (11/03/24)
minPerGroup <- 15L
meth <- methylKit::unite(norm.filt.dat, destrand = F, min.per.group = minPerGroup)

# get percent methylation df
percent_meth <- percMethylation(meth, rowids = TRUE)

# save percent methylation data
saveRDS(percent_meth, "./methylation_results/percent_meth.RDATA")



# Cluster samples
dendrogram <- clusterSamples(percent_meth,
                             dist = "correlation",
                             method = "average",
                             plot = T)

# pca plot


pca.out <- prcomp(t(percent_meth %>% na.omit()), scale. = FALSE)

plot(pca.out$x[,1:2])

# plot nice PCA plot

pca.out$x[,1:2] %>% 
        as.data.frame() %>% 
        rownames_to_column(var = 'ID') %>% 
        mutate(Sample = factor(ifelse(grepl('PRE-ASAT', ID), 'PRE-ASAT', 
                                      ifelse(grepl('POST-ASAT', ID), 'POST-ASAT', 
                                             ifelse(grepl('PRE-QF', ID), 'PRE-QF', 'POST-QF'))),
                               levels = c('PRE-ASAT','POST-ASAT','PRE-QF','POST-QF'))) %>% 
        mutate(Tissue = ifelse(grepl('ASAT', ID), 'ASAT', 'QF')) %>% 
        filter(Tissue == "ASAT") %>% 
        ggplot(aes(x = PC1, y = PC2))+
        #geom_point(aes(size = 2+pca.out$sdev),shape = 3)+
        geom_point(aes(fill = Sample), size = 3, stroke = 2, color = 'Black', shape = 21)+
        scale_fill_manual(values = c("#440154FF","#B396B9","#5DC863FF","#B3D7B1"))+
        theme_classic(base_size = 20)+
        labs(x = 'Principal Component 1',
             y = 'Principal Component 2')+
        theme(legend.position = 'top')+
        geom_label_repel(aes(label = ID), max.overlaps = 68)+
        guides(fill = guide_legend(override.aes = list(size = 6)))


##########################################################################################################
#########               save files                                              ####################
##########################################################################################################







