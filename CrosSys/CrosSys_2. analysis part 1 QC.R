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

head(density_df)

density_df %>% 
        ggplot(aes(x = percent_meth, group = ID, color = Tissue))+
        geom_density()




##########################################################################################################
#########               normalize                                                     ####################
##########################################################################################################



# Combine into a methylRawList
# Convert data frames to methylRaw objects
methylRawList_obj <- lapply(seq_along(samples), function(i) {
        new("methylRaw",
                samples[[i]],
                sample.id = sample_names[i],
                assembly = "hg38",
                context = "CpG",
                resolution = "base"
        )
})


# Combine into a methylRawList
methylRawList_obj <- new("methylRawList", methylRawList_obj)

# Normalise coverage
norm.filt.dat <- normalizeCoverage(methylRawList_obj)

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
#########               PCA before normalization                                                     ####################
##########################################################################################################

# merge percent meth into one df

# filter probes for already IDd probes

united_ID <- rownames(percent_meth)

percent_meth_non.norm <- data.frame(ID = united_ID)
        
for (i in 1:length(samples)) {
        x <- samples[[i]] %>% 
                mutate(ID = paste0(chr,".",start,".",end)) %>% 
                dplyr::select(ID, percent_meth) %>% 
                filter(ID %in% united_ID)
        
        percent_meth_non.norm <- merge(percent_meth_non.norm, x, by = "ID")
        
        print(i)
        
}

# rename colnames of samples

colnames(percent_meth_non.norm) <- c('ID', names(samples))

non.norm_PCA <- prcomp(t(percent_meth_non.norm[,2:69]))

non.norm_PCA$x[,1:2] %>% 
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
        scale_fill_manual(values = c("#440154FF","#B396B9","#5DC863FF","#B3D7B1"))+
        theme_classic(base_size = 20)+
        labs(x = 'Principal Component 1',
             y = 'Principal Component 2')+
        theme(legend.position = 'top')+
        # geom_label_repel(aes(label = ID), max.overlaps = 68)+
        guides(fill = guide_legend(override.aes = list(size = 6)))

# PC1 perfectly separates QF and ASAT samples, as expected

percent_meth_norm <- data.frame(ID = united_ID)

for (i in 1:length(norm.filt.dat)) {
        
        x <- data.frame(ID = paste0(norm.filt.dat[[i]]$chr, ".",
                               norm.filt.dat[[i]]$start, ".",
                               norm.filt.dat[[i]]$end),
                   percent_meth = norm.filt.dat[[i]]$percent_meth) %>% 
                filter(ID %in% united_ID)
      
        
        percent_meth_norm <- merge(percent_meth_norm, x, by = "ID")
        
        print(i)
        
}

colnames(percent_meth_norm) <- c('ID', names(samples))

norm_PCA <- prcomp(t(percent_meth_norm[,2:69]))

norm_PCA$x[,1:2] %>% 
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
        scale_fill_manual(values = c("#440154FF","#B396B9","#5DC863FF","#B3D7B1"))+
        theme_classic(base_size = 20)+
        labs(x = 'Principal Component 1',
             y = 'Principal Component 2')+
        theme(legend.position = 'top')+
        # geom_label_repel(aes(label = ID), max.overlaps = 68)+
        guides(fill = guide_legend(override.aes = list(size = 6)))




##########################################################################################################
#########               separate samples into QF and ASAT and perform dmr analysis separate                                                    ####################
##########################################################################################################

# treatment groups
# 
# treatment <- c()
# 
# treatment[grep('PRE-ASAT', sample_names)] <- 0
# treatment[grep('POST-ASAT', sample_names)] <- 1
# treatment[grep('PRE-QF', sample_names)] <- 2
# treatment[grep('POST-QF', sample_names)] <- 3

# split normalized and filtered data by condition (pre-post | QF-ASAT)
split_list <- split(norm.filt.dat, norm.filt.dat@treatment)

# Access the subsets
ASAT_df <- unlist(split_list[1:2])  # Treatment group 0 and 1 (ASAT)
QF_df <- split_list[3:4] # Treatment group 2 and 3 (QF)

# Merge everything together based on common bases
# Only merge CpGs that are in at least 6 samples per group after discussion with Adam (11/03/24)
minPerGroup <- 15L
meth_ASAT <- methylKit::unite(ASAT_df, destrand = F, min.per.group = minPerGroup)

















# annotate probe info

# Annotate bases 
BiocManager::install("annotatr")
library(annotatr)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("Homo.sapiens")
BiocManager::install("ChIPpeakAnno", force = TRUE)


# Select annotations for intersection with regions
# Note inclusion of custom annotation, and use of shortcuts
annots = c('hg38_basicgenes')

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg38', annotations = annots)

loadings <- loadings %>% 
        as.data.frame() %>% 
        rowwise() %>% 
        mutate(chr = str_split(ID, pattern = "\\.")[[1]][1],
               start = str_split(ID, pattern = "\\.")[[1]][2],
               end = str_split(ID, pattern = "\\.")[[1]][3]) %>% 
        ungroup()
str(loadings)
loadings$chr <- as.numeric(loadings$chr)
loadings$start <- as.numeric(loadings$start)
loadings$end <- as.numeric(loadings$end)


############################################
#### create annotation set

BiocManager::install("rtracklayer", force = TRUE)
library(rtracklayer)


GRC38 <- readGFF("methylation_results/Homo_sapiens.GRCh38.112.chr.gtf")



##########################################################################################################
#########               Annotate with genomation                                                     ####################
##########################################################################################################

# # Install dependencies
# install.packages( c("data.table","plyr","reshape2","ggplot2","gridBase","devtools"))
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("GenomicRanges","rtracklayer","impute","Rsamtools"))
# BiocManager::install("impute", force = TRUE)
# 
# Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/rtools40/usr/bin", sep = ";"))
# Sys.setenv(BINPREF = "C:/rtools40/mingw_$(WIN)/bin/")
# 
# # install the packages
# library(devtools)
# devtools::install_github("BIMSBbioinfo/genomation",build_vignettes=FALSE)
build_annotations(GRC38)

library(genomation)

# Annotate with gene annotations
gene_anno <- annotateWithGeneParts(target = meth, feature = GRC38_granges)



# Extract and analyze annotation results
gene_anno_results <- gene_anno$annotation
summary(gene_anno_results)





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



##########################################################################################################
#########               save files                                              ####################
##########################################################################################################


save_RDS(percent_meth, "./methylation_results/percent_meth.RDATA")

saveRDS(norm.filt.dat, "./methylation_results/norm.filt.dat.rds")

saveRDS(samples, "./methylation_results/samples.RDATA")

norm.filt.dat <- readRDS(file =  "./methylation_results/norm.filt.dat.rds")

samples <- readRDS("./methylation_results/samples.RDATA")

percent_meth <- readRDS("./methylation_results/percent_meth.RDATA")
