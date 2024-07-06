#
#               CrosSys
#
#               combine data, PCA, normalization etc. 
#



library(methylKit)
library(tidyverse)
library(ggplot2)
library(cowplot)


# Set working directory max's mac
setwd("/Users/maxullrich/OneDrive - UGent/CrosSys")


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

for (i in 1:length(samples)) {
        
        samples[[i]] <- samples[[i]] %>% 
                mutate(check = (numCs+numTs == coverage)) %>% 
                filter(check != FALSE) %>% 
                mutate(percent_meth = numCs/coverage) %>% 
                dplyr::select(chr, start, end, strand, coverage, numCs, numTs, percent_meth)
        print(i)
}


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










