#
#
#   Align and preprocess CrosSys RRBS data
#
#


# Script for analysis of RRBS data pre-processed using Bismark
# This script is specific for CrosSys, data originating from Finnish group

# Load packages
library(methylKit)
library(tidyverse)
library(ggplot2)

# Set working directory
setwd("/Users/maxul/Documents/Skole/Lab/CrosSys/mapping_results/")

# Locate all files
inFiles <- as.list(list.files(pattern = '*.cov.gz',
                              full.names = T))


# get names and make treatment vector

samples <- list.files(pattern = '*.cov.gz',
                      full.names = F)

treatments <- c()

treatments[grep('PRE-ASAT', samples)]  <- 0
treatments[grep('POST-ASAT', samples)] <- 1
treatments[grep('PRE-QF', samples)]    <- 2
treatments[grep('POST-QF', samples)]   <- 3




# Read in data using methRead
dat <- methRead(inFiles,
                sample.id = as.list(samples),
                assembly = 'hg38',
                pipeline = 'bismarkCoverage',
                header = F,
                treatment = treatments,
                mincov = 10)


# Filter low coverage (already filtered) and high coverage
filt.dat = filterByCoverage(
        dat,
        lo.count = 10,
        lo.perc = NULL,
        hi.count = NULL,
        hi.perc = 99.9
)

# Normalise coverage
norm.filt.dat <- normalizeCoverage(filt.dat)


# extract methylation data
filt.meth.data <- c()

for (i in 1:length(filt.dat)) {
        

        x <- data.frame(chr = norm.filt.dat[[i]]$chr,
                         start = norm.filt.dat[[i]]$start,
                         end = norm.filt.dat[[i]]$end, 
                         strand = norm.filt.dat[[i]]$strand,
                         coverage = norm.filt.dat[[i]]$coverage, 
                         numCs = norm.filt.dat[[i]]$numCs,
                         numTs = norm.filt.dat[[i]]$numTs) %>% 
                mutate(check = numTs+numCs==coverage) %>% 
                filter(check == TRUE) %>% 
                dplyr::select(1:7) %>% 
                mutate(percent_meth = numCs/coverage)
        
        filt.meth.data[[i]] <- x
        
        names(filt.meth.data)[i] <- sub("_R[12].*", "", samples[i])
        print(i)
        
}

#####################

### merge into one dataframe



meth_df <- data.frame("ID" = "NA")

for (i in 1:length(filt.meth.data)) {
        df <- data.frame(chr = filt.meth.data[[i]]$chr,
                         start = filt.meth.data[[i]]$start,
                         coverage = filt.meth.data[[i]]$coverage, 
                         numCs = filt.meth.data[[i]]$numCs,
                         numTs = filt.meth.data[[i]]$numTs) %>% 
                mutate(check = numTs+numCs==coverage) %>% 
                filter(check == TRUE) %>%
                mutate(percent_meth = numCs/coverage, 
                       ID = paste(chr, start, sep = "_")) %>% 
                dplyr::select(ID, percent_meth) %>% 
                distinct()
        
        names(df) <- c("ID", names(filt.meth.data[i]))
        
        meth_df <- merge(meth_df, df, by = "ID", all = TRUE)
        
        print(i)
        
}

meth_df %>% 
        pivot_longer(cols = 2:length(colnames(meth_df)), names_to = "Sample", values_to = "percent_meth") %>% 
        na.omit() %>% 
        group_by(Sample) %>% 
        summarize(number = length(percent_meth),
                  mean = mean(percent_meth),
                  min = min(percent_meth),
                  max = max(percent_meth),
                  sd = sd(percent_meth)) -> summary

summary


# plot density plot

meth_df %>% 
        pivot_longer(cols = 2:14, names_to = "Sample", values_to = "percent_meth") %>% 
        na.omit() %>% 
        ggplot(aes(x = percent_meth, group = Sample, color = Sample))+
        geom_density()+
        theme_classic(base_size = 20)+
        labs(title = "normalized percent meth distribution")


# count NA per row of meth_df

table(rowSums(is.na(meth_df)))

sum(dat[[15]]$coverage)



### save file to send

write.csv(df, "./first_sample_df.csv")






##########################################################################################################
#########           Jameses run on RMA samples             ###############################################
##########################################################################################################







# Locate all files
inFiles <- as.list(list.files(path = 'bismark_cov_files',
                              pattern = '*.cov.gz',
                              full.names = T))

# Create sample names based on filenames 
samples <-
        gsub(
                'bismark_cov_files/GC-AS-10308-',
                '',
                list.files(
                        path = 'bismark_cov_files',
                        pattern = '*.cov.gz',
                        full.names = T
                )
        )
samples <- gsub('_merged_R1_processed_val_1_bismark_bt2_pe_MOD.deduplicated.bismark.cov.gz', '', samples)

# Make a treatment vector
# 0 = control (baseline?) (B)
# 1 = atrophy (A1)
# 2 = recovery (R)
# 3 = repeated atrophy (A2)
# Base this on sample name
treatments <- c()
treatments[grep('B', samples)] <- 0
treatments[grep('A1', samples)] <- 1
treatments[grep('R', samples)] <- 2
treatments[grep('A2', samples)] <- 3

# Read in data using methRead
dat <- methRead(inFiles,
                sample.id = as.list(samples),
                assembly = 'hg38',
                pipeline = 'bismarkCoverage',
                header = F,
                treatment = treatments,
                mincov = 10)

# Set names of list elements
names(dat) <- samples

# Plot distribution of methylation per base
for (x in samples) {
        tiff(
                filename = paste0(x, '_hist.tiff'),
                height = 8,
                width = 8,
                units = 'in',
                res = 300
        )
        getMethylationStats(dat[[x]], plot = TRUE, both.strands = FALSE)
        dev.off()
}

# Filter low coverage (already filtered) and high coverage
filt.dat = filterByCoverage(
        dat,
        lo.count = 10,
        lo.perc = NULL,
        hi.count = NULL,
        hi.perc = 99.9
)

# Normalise coverage
norm.filt.dat <- normalizeCoverage(filt.dat)

# Merge everything together based on common bases
meth <- unite(norm.filt.dat, destrand = F)

# Cluster samples
dendrogram <- clusterSamples(meth,
                             dist = "correlation",
                             method = "average",
                             plot = F)
tiff(filename = 'dendrogram.tiff',
     height = 8,
     width = 8,
     units = 'in',
     res = 300)
plot(dendrogram)
dev.off()

# PCA samples
pca <- PCASamples(meth, obj.return = T)

# Extract first two PCs
pcaPlt <- as.data.frame(pca$x[, c(1,2)])

# Add meta
ID <- gsub('A1', '',
           gsub('A2', '',
                gsub('B', '', 
                     gsub('R', '', samples))))

treatmentVerb <- c()
treatmentVerb[treatments == 0] <- 'B'
treatmentVerb[treatments == 1] <- 'A1'
treatmentVerb[treatments == 2] <- 'R'
treatmentVerb[treatments == 3] <- 'A2'

pcaPlt$Treatment <- treatmentVerb
pcaPlt$ID <- ID

# Plot
ggplot(pcaPlt, aes(x = PC1, y = PC2, colour = ID, shape = Treatment)) +
        geom_point(size = 5)
ggsave(plot = last_plot(),
       filename = 'pca.tiff',
       height = 8,
       width = 8,
       units = 'in',
       dpi = 300)

# Subset data to perform pairwise comparisons:
# Make a list of subsets
comparisons <- list()

# a.  Atrophy compared to Baseline
comparisons[['A1-B']] <- reorganize(meth,
                                    sample.ids = meth@sample.ids[meth@treatment %in% c(1, 0)],
                                    treatment = meth@treatment[meth@treatment %in% c(1, 0)])
# b.  Recovery compared to Baseline
comparisons[['R-B']] <- reorganize(meth,
                                   sample.ids = meth@sample.ids[meth@treatment %in% c(2, 0)],
                                   treatment = meth@treatment[meth@treatment %in% c(2, 0)])
# c.  Repeated atrophy compared to Baseline
comparisons[['A2-B']] <- reorganize(meth,
                                    sample.ids = meth@sample.ids[meth@treatment %in% c(3, 0)],
                                    treatment = meth@treatment[meth@treatment %in% c(3, 0)])
# d.  Recovery compared to atrophy
comparisons[['R-A1']] <- reorganize(meth,
                                    sample.ids = meth@sample.ids[meth@treatment %in% c(2, 1)],
                                    treatment = meth@treatment[meth@treatment %in% c(2, 1)])
# e.  Repeated atrophy compared to Recovery
comparisons[['A2-R']] <- reorganize(meth,
                                    sample.ids = meth@sample.ids[meth@treatment %in% c(3, 2)],
                                    treatment = meth@treatment[meth@treatment %in% c(3, 2)])
# f.  Repeated atrophy compared to atrophy
comparisons[['A2-A1']] <- reorganize(meth,
                                     sample.ids = meth@sample.ids[meth@treatment %in% c(3, 1)],
                                     treatment = meth@treatment[meth@treatment %in% c(3, 1)])

# Calculate differential methylation for all comparisons
diffMethResults <- list()
for (comp in names(comparisons)) {
        diffMethResults[[comp]] <- calculateDiffMeth(comparisons[[comp]])
}

# Annotate bases 
library(annotatr)

# Select annotations for intersection with regions
# Note inclusion of custom annotation, and use of shortcuts
annots = c('hg38_basicgenes')

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg38', annotations = annots)

# Loop to go through all comparisons, annotate, export results
for (comp in names(diffMethResults)) {
        print(paste('Starting for compairson:', comp))
        # Intersect the regions we read in with the annotations
        test <- as(diffMethResults[[comp]], "GRanges")
        print('Annotating...')
        dm_annotated = annotate_regions(
                regions = test,
                annotations = annotations,
                ignore.strand = TRUE,
                quiet = FALSE
        )
        
        print('Converting to dataframe...')
        test_df = as(dm_annotated, "data.frame")
        
        print('Exporting at CSV...')
        write.csv(test_df, file = paste0(comp, '_DM_Result.csv'))
}





#############################################################################

### Max continuation

#############################################################################

library(methylKit)
library(tidyverse)



# create dataframe with one participant and annotate

# filter out probes where number of Ts and Cs does not equal coverage (eks. ~2000 in 12A1)

df <- data.frame(chr = norm.filt.dat[[1]]$chr,
                 start = norm.filt.dat[[1]]$start,
                 end = norm.filt.dat[[1]]$end, 
                 strand = norm.filt.dat[[1]]$strand,
                 coverage = norm.filt.dat[[1]]$coverage, 
                 numCs = norm.filt.dat[[1]]$numCs,
                 numTs = norm.filt.dat[[1]]$numTs) %>% 
        mutate(check = numTs+numCs==coverage) %>% 
        filter(check == TRUE) %>% 
        dplyr::select(1:7) %>% 
        mutate(percent_meth = numCs/coverage)

# try to annotate

BiocManager::install("annotatr")
library(annotatr)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

# Select annotations for intersection with regions
# Note inclusion of custom annotation, and use of shortcuts
annots = c('hg38_basicgenes')

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg38', annotations = annots)


# annotate df
dm_annotated = annotate_regions(
        regions = as(df, "GRanges"),
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = FALSE
)

# convert annotated data to data.frame
dm_annotated %>% as.data.frame() -> x



saveRDS(annotations, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/annotations.RDATA")

# create list of norm.filt.dat percent meth data

FP12 <- list()
FP15 <- list()
FP27 <- list()
FP3 <- list()
FP4 <- list()
FP5 <- list()
FP7 <- list()
FP8 <- list()
FP9 <- list()

for (i in 1:4) {
        
        print(paste("calculating % meth",names(norm.filt.dat[i])))
        
        df <- data.frame(chr = norm.filt.dat[[i]]$chr,
                         start = norm.filt.dat[[i]]$start,
                         end = norm.filt.dat[[i]]$end, 
                         strand = norm.filt.dat[[i]]$strand,
                         coverage = norm.filt.dat[[i]]$coverage, 
                         numCs = norm.filt.dat[[i]]$numCs,
                         numTs = norm.filt.dat[[i]]$numTs) %>% 
                mutate(check = numTs+numCs==coverage) %>% 
                filter(check == TRUE) %>% 
                dplyr::select(1:7) %>% 
                mutate(percent_meth = numCs/coverage)
        
        print(paste("annotating",names(norm.filt.dat[i])))
        
        # annotate df
        dm_annotated = annotate_regions(
                regions = as(df, "GRanges"),
                annotations = annotations,
                ignore.strand = TRUE,
                quiet = FALSE
        )
        
        # convert annotated data to data.frame
        dm_annotated %>% as.data.frame() -> x
        
        FP12[[names(norm.filt.dat[i])]] <- x
        
        print(paste("done with",i, "/ 36", "   --- sample:", names(norm.filt.dat[i]))) # tell me how far i have come
}

saveRDS(FP12, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP12_meth.RDATA")

FP15 <- list()

for (i in 5:8) {
        
        print(paste("calculating % meth",names(norm.filt.dat[i])))
        
        df <- data.frame(chr = norm.filt.dat[[i]]$chr,
                         start = norm.filt.dat[[i]]$start,
                         end = norm.filt.dat[[i]]$end, 
                         strand = norm.filt.dat[[i]]$strand,
                         coverage = norm.filt.dat[[i]]$coverage, 
                         numCs = norm.filt.dat[[i]]$numCs,
                         numTs = norm.filt.dat[[i]]$numTs) %>% 
                mutate(check = numTs+numCs==coverage) %>% 
                filter(check == TRUE) %>% 
                dplyr::select(1:7) %>% 
                mutate(percent_meth = numCs/coverage)
        
        print(paste("annotating",names(norm.filt.dat[i])))
        
        # annotate df
        dm_annotated = annotate_regions(
                regions = as(df, "GRanges"),
                annotations = annotations,
                ignore.strand = TRUE,
                quiet = FALSE
        )
        
        # convert annotated data to data.frame
        dm_annotated %>% as.data.frame() -> x
        
        FP15[[names(norm.filt.dat[i])]] <- x
        
        print(paste("done with",i, "/ 36", "   --- sample:", names(norm.filt.dat[i]))) # tell me how far i have come
}

saveRDS(FP15, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP15_meth.RDATA")

FP27 <- list()

for (i in 9:12) {
        
        print(paste("calculating % meth",names(norm.filt.dat[i])))
        
        df <- data.frame(chr = norm.filt.dat[[i]]$chr,
                         start = norm.filt.dat[[i]]$start,
                         end = norm.filt.dat[[i]]$end, 
                         strand = norm.filt.dat[[i]]$strand,
                         coverage = norm.filt.dat[[i]]$coverage, 
                         numCs = norm.filt.dat[[i]]$numCs,
                         numTs = norm.filt.dat[[i]]$numTs) %>% 
                mutate(check = numTs+numCs==coverage) %>% 
                filter(check == TRUE) %>% 
                dplyr::select(1:7) %>% 
                mutate(percent_meth = numCs/coverage)
        
        print(paste("annotating",names(norm.filt.dat[i])))
        
        # annotate df
        dm_annotated = annotate_regions(
                regions = as(df, "GRanges"),
                annotations = annotations,
                ignore.strand = TRUE,
                quiet = FALSE
        )
        
        # convert annotated data to data.frame
        dm_annotated %>% as.data.frame() -> x
        
        FP27[[names(norm.filt.dat[i])]] <- x
        
        print(paste("done with",i, "/ 36", "   --- sample:", names(norm.filt.dat[i]))) # tell me how far i have come
}


saveRDS(FP27, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP27_meth.RDATA")


FP3 <- list()

for (i in 13:16) {
        
        print(paste("calculating % meth",names(norm.filt.dat[i])))
        
        df <- data.frame(chr = norm.filt.dat[[i]]$chr,
                         start = norm.filt.dat[[i]]$start,
                         end = norm.filt.dat[[i]]$end, 
                         strand = norm.filt.dat[[i]]$strand,
                         coverage = norm.filt.dat[[i]]$coverage, 
                         numCs = norm.filt.dat[[i]]$numCs,
                         numTs = norm.filt.dat[[i]]$numTs) %>% 
                mutate(check = numTs+numCs==coverage) %>% 
                filter(check == TRUE) %>% 
                dplyr::select(1:7) %>% 
                mutate(percent_meth = numCs/coverage)
        
        print(paste("annotating",names(norm.filt.dat[i])))
        
        # annotate df
        dm_annotated = annotate_regions(
                regions = as(df, "GRanges"),
                annotations = annotations,
                ignore.strand = TRUE,
                quiet = FALSE
        )
        
        # convert annotated data to data.frame
        dm_annotated %>% as.data.frame() -> x
        
        FP3[[names(norm.filt.dat[i])]] <- x
        
        print(paste("done with",i, "/ 36", "   --- sample:", names(norm.filt.dat[i]))) # tell me how far i have come
}


saveRDS(FP3, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP3_meth.RDATA")


FP4 <- list()

for (i in 17:20) {
        
        print(paste("calculating % meth",names(norm.filt.dat[i])))
        
        df <- data.frame(chr = norm.filt.dat[[i]]$chr,
                         start = norm.filt.dat[[i]]$start,
                         end = norm.filt.dat[[i]]$end, 
                         strand = norm.filt.dat[[i]]$strand,
                         coverage = norm.filt.dat[[i]]$coverage, 
                         numCs = norm.filt.dat[[i]]$numCs,
                         numTs = norm.filt.dat[[i]]$numTs) %>% 
                mutate(check = numTs+numCs==coverage) %>% 
                filter(check == TRUE) %>% 
                dplyr::select(1:7) %>% 
                mutate(percent_meth = numCs/coverage)
        
        print(paste("annotating",names(norm.filt.dat[i])))
        
        # annotate df
        dm_annotated = annotate_regions(
                regions = as(df, "GRanges"),
                annotations = annotations,
                ignore.strand = TRUE,
                quiet = FALSE
        )
        
        # convert annotated data to data.frame
        dm_annotated %>% as.data.frame() -> x
        
        FP4[[names(norm.filt.dat[i])]] <- x
        
        print(paste("done with",i, "/ 36", "   --- sample:", names(norm.filt.dat[i]))) # tell me how far i have come
}
saveRDS(FP4, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP4_meth.RDATA")


FP5 <- list()

for (i in 21:24) {
        
        print(paste("calculating % meth",names(norm.filt.dat[i])))
        
        df <- data.frame(chr = norm.filt.dat[[i]]$chr,
                         start = norm.filt.dat[[i]]$start,
                         end = norm.filt.dat[[i]]$end, 
                         strand = norm.filt.dat[[i]]$strand,
                         coverage = norm.filt.dat[[i]]$coverage, 
                         numCs = norm.filt.dat[[i]]$numCs,
                         numTs = norm.filt.dat[[i]]$numTs) %>% 
                mutate(check = numTs+numCs==coverage) %>% 
                filter(check == TRUE) %>% 
                dplyr::select(1:7) %>% 
                mutate(percent_meth = numCs/coverage)
        
        print(paste("annotating",names(norm.filt.dat[i])))
        
        # annotate df
        dm_annotated = annotate_regions(
                regions = as(df, "GRanges"),
                annotations = annotations,
                ignore.strand = TRUE,
                quiet = FALSE
        )
        
        # convert annotated data to data.frame
        dm_annotated %>% as.data.frame() -> x
        
        FP5[[names(norm.filt.dat[i])]] <- x
        
        print(paste("done with",i, "/ 36", "   --- sample:", names(norm.filt.dat[i]))) # tell me how far i have come
}

saveRDS(FP5, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP5_meth.RDATA")


# FP7
names(norm.filt.dat)

for (i in 25:28) {
        
        print(paste("calculating % meth",names(norm.filt.dat[i])))
        
        df <- data.frame(chr = norm.filt.dat[[i]]$chr,
                         start = norm.filt.dat[[i]]$start,
                         end = norm.filt.dat[[i]]$end, 
                         strand = norm.filt.dat[[i]]$strand,
                         coverage = norm.filt.dat[[i]]$coverage, 
                         numCs = norm.filt.dat[[i]]$numCs,
                         numTs = norm.filt.dat[[i]]$numTs) %>% 
                mutate(check = numTs+numCs==coverage) %>% 
                filter(check == TRUE) %>% 
                dplyr::select(1:7) %>% 
                mutate(percent_meth = numCs/coverage)
        
        print(paste("annotating",names(norm.filt.dat[i])))
        
        # annotate df
        dm_annotated = annotate_regions(
                regions = as(df, "GRanges"),
                annotations = annotations,
                ignore.strand = TRUE,
                quiet = FALSE
        )
        
        # convert annotated data to data.frame
        dm_annotated %>% as.data.frame() -> x
        
        FP7[[names(norm.filt.dat[i])]] <- x
        
        print(paste("done with",i, "/ 36", "   --- sample:", names(norm.filt.dat[i]))) # tell me how far i have come
}

saveRDS(FP7, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP7_meth.RDATA")



# FP8

for (i in 29:32) {
        
        print(paste("calculating % meth",names(norm.filt.dat[i])))
        
        df <- data.frame(chr = norm.filt.dat[[i]]$chr,
                         start = norm.filt.dat[[i]]$start,
                         end = norm.filt.dat[[i]]$end, 
                         strand = norm.filt.dat[[i]]$strand,
                         coverage = norm.filt.dat[[i]]$coverage, 
                         numCs = norm.filt.dat[[i]]$numCs,
                         numTs = norm.filt.dat[[i]]$numTs) %>% 
                mutate(check = numTs+numCs==coverage) %>% 
                filter(check == TRUE) %>% 
                dplyr::select(1:7) %>% 
                mutate(percent_meth = numCs/coverage)
        
        print(paste("annotating",names(norm.filt.dat[i])))
        
        # annotate df
        dm_annotated = annotate_regions(
                regions = as(df, "GRanges"),
                annotations = annotations,
                ignore.strand = TRUE,
                quiet = FALSE
        )
        
        # convert annotated data to data.frame
        dm_annotated %>% as.data.frame() -> x
        
        FP8[[names(norm.filt.dat[i])]] <- x
        
        print(paste("done with",i, "/ 36", "   --- sample:", names(norm.filt.dat[i]))) # tell me how far i have come
}

saveRDS(FP8, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP8_meth.RDATA")


#FP9

for (i in 33:36) {
        
        print(paste("calculating % meth",names(norm.filt.dat[i])))
        
        df <- data.frame(chr = norm.filt.dat[[i]]$chr,
                         start = norm.filt.dat[[i]]$start,
                         end = norm.filt.dat[[i]]$end, 
                         strand = norm.filt.dat[[i]]$strand,
                         coverage = norm.filt.dat[[i]]$coverage, 
                         numCs = norm.filt.dat[[i]]$numCs,
                         numTs = norm.filt.dat[[i]]$numTs) %>% 
                mutate(check = numTs+numCs==coverage) %>% 
                filter(check == TRUE) %>% 
                dplyr::select(1:7) %>% 
                mutate(percent_meth = numCs/coverage)
        
        print(paste("annotating",names(norm.filt.dat[i])))
        
        # annotate df
        dm_annotated = annotate_regions(
                regions = as(df, "GRanges"),
                annotations = annotations,
                ignore.strand = TRUE,
                quiet = FALSE
        )
        
        # convert annotated data to data.frame
        dm_annotated %>% as.data.frame() -> x
        
        FP9[[names(norm.filt.dat[i])]] <- x
        
        print(paste("done with",i, "/ 36", "   --- sample:", names(norm.filt.dat[i]))) # tell me how far i have come
}

saveRDS(FP9, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP9_meth.RDATA")



########################################

### read annotated meth data

FP12 <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP12_meth.RDATA")
FP15 <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP15_meth.RDATA")
FP27 <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP27_meth.RDATA")
FP3 <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP3_meth.RDATA")
FP4 <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP4_meth.RDATA")
FP5 <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP5_meth.RDATA")
FP7 <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP7_meth.RDATA")
FP8 <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP8_meth.RDATA")
FP9 <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/FP9_meth.RDATA")





#######################################################

# check unique probes





FP8$`8A1` %>% 
        dplyr::select(1:2) %>% 
        mutate(ID = paste(seqnames, start, sep = "_")) %>% 
        distinct(ID)


#############################

# plot raw data

df_raw <- data.frame("ID" = "NA", 
                     "percent_meth" = "NA",
                     "Sample" = "NA")

for (i in 1:length(dat)) {
        df <- data.frame(chr = dat[[i]]$chr,
                         start = dat[[i]]$start,
                         coverage = dat[[i]]$coverage, 
                         numCs = dat[[i]]$numCs) %>% 
                #head(100) %>% 
                mutate(percent_meth = numCs/coverage, 
                       ID = paste(chr, start, sep = "_")) %>% 
                dplyr::select(ID, percent_meth) %>% 
                mutate(Sample = names(dat[i]))        
        
        df_raw <- rbind(df_raw, df)
        
        print(i)
        
}

df_raw <- df_raw[-1,]
df_raw$percent_meth <- as.numeric(df_raw$percent_meth)



ggplot(df_raw, aes(x = percent_meth, color = Sample)) +
        geom_density() +
        theme_minimal() +
        labs(title = "Raw Density Plot", x = "Percent_meth", y = "Density")




# plot filtered data

df_filt <- data.frame("ID" = "NA", 
                      "percent_meth" = "NA",
                      "Sample" = "NA")

for (i in 1:length(norm.filt.dat)) {
        df <- data.frame(chr = norm.filt.dat[[i]]$chr,
                         start = norm.filt.dat[[i]]$start,
                         coverage = norm.filt.dat[[i]]$coverage, 
                         numCs = norm.filt.dat[[i]]$numCs) %>% 
                #head(100) %>% 
                mutate(percent_meth = numCs/coverage, 
                       ID = paste(chr, start, sep = "_")) %>% 
                dplyr::select(ID, percent_meth) %>% 
                mutate(Sample = names(dat[i]))        
        
        df_filt <- rbind(df_filt, df)
        
        print(i)
        
}

df_filt <- df_filt[-1,]
df_filt$percent_meth <- as.numeric(df_filt$percent_meth)



ggplot(df_filt, aes(x = percent_meth, color = Sample)) +
        geom_density() +
        theme_minimal() +
        labs(title = "Filtered Density Plot", x = "Percent_meth", y = "Density")



#####################

### merge into one dataframe



meth_df <- data.frame("ID" = "NA")

for (i in 1:length(norm.filt.dat)) {
        df <- data.frame(chr = norm.filt.dat[[i]]$chr,
                         start = norm.filt.dat[[i]]$start,
                         coverage = norm.filt.dat[[i]]$coverage, 
                         numCs = norm.filt.dat[[i]]$numCs,
                         numTs = norm.filt.dat[[i]]$numTs) %>% 
                mutate(check = numTs+numCs==coverage) %>% 
                filter(check == TRUE) %>%
                mutate(percent_meth = numCs/coverage, 
                       ID = paste(chr, start, sep = "_")) %>% 
                dplyr::select(ID, percent_meth) %>% 
                distinct()
        
        names(df) <- c("ID", names(norm.filt.dat[i]))
        
        meth_df <- merge(meth_df, df, by = "ID", all = TRUE)
        
        print(i)
        
}



saveRDS(meth_df, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/math_df.RDATA")
meth_df <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/math_df.RDATA")

# filter out rows/probes missing more than 4 measurements in at least one timepoint






filtered_df <- meth_df %>% 
        dplyr::select("ID", "3B", "4B", "5B", "7B", "8B", "9B", "12B", "15B", "27B",
                      "3A1", "4A1", "5A1", "7A1", "8A1", "9A1", "12A1", "15A1", "27A1",
                      "3R", "4R", "5R", "7R", "8R", "9R", "12R", "15R", "27R",
                      "3A2", "4A2", "5A2", "7A2", "8A2", "9A2", "12A2", "15A2", "27A2") 


filtered_df$baseline <- rowSums(is.na(filtered_df[,2:10]))
filtered_df$atrophy1 <- rowSums(is.na(filtered_df[,11:19]))

# ran out of memory, remove preliminary rows, and continue

filtered_df <- filtered_df %>% 
        filter(baseline <= 4)

filtered_df <- filtered_df %>% 
        filter(atrophy1 <= 4)


filtered_df$recovery <- rowSums(is.na(filtered_df[,20:28]))
filtered_df$atrophy2 <- rowSums(is.na(filtered_df[,29:37]))

filtered_df <- filtered_df %>% 
        filter(recovery <= 4)

filtered_df <- filtered_df %>% 
        filter(atrophy2 <= 4)




##########################################

### impute remaining NAs





# split data into timepoints, impute mean, and recombine data

baseline <- filtered_df[,2:10]
atrophy1 <- filtered_df[,11:19]
recovery <- filtered_df[,20:28]
atrophy2 <- filtered_df[,29:37]



baseline_imputeNA <- apply(baseline, 1, function(row) {
        if (any(is.na(row))) {
                non_na_values <- na.omit(as.numeric(row))
                mean_val <- sum(non_na_values) / length(non_na_values)
                row[is.na(row)] <- mean_val
        }
        return(row)
})

baseline_imputeNA <- as.data.frame(t(baseline_imputeNA))


atrophy1_imputeNA <- apply(atrophy1, 1, function(row) {
        if (any(is.na(row))) {
                non_na_values <- na.omit(as.numeric(row))
                mean_val <- sum(non_na_values) / length(non_na_values)
                row[is.na(row)] <- mean_val
        }
        return(row)
})

atrophy1_imputeNA <- as.data.frame(t(atrophy1_imputeNA))


recovery_imputeNA <- apply(recovery, 1, function(row) {
        if (any(is.na(row))) {
                non_na_values <- na.omit(as.numeric(row))
                mean_val <- sum(non_na_values) / length(non_na_values)
                row[is.na(row)] <- mean_val
        }
        return(row)
})

recovery_imputeNA <- as.data.frame(t(recovery_imputeNA))




atrophy2_imputeNA <- apply(atrophy2, 1, function(row) {
        if (any(is.na(row))) {
                non_na_values <- na.omit(as.numeric(row))
                mean_val <- sum(non_na_values) / length(non_na_values)
                row[is.na(row)] <- mean_val
        }
        return(row)
})

atrophy2_imputeNA <- as.data.frame(t(atrophy2_imputeNA))


### 

### merge dataframes

meth_imputeNA = data.frame(ID = filtered_df$ID)

meth_imputeNA <- cbind(meth_imputeNA, baseline_imputeNA)

meth_imputeNA <- cbind(meth_imputeNA, atrophy1_imputeNA)

meth_imputeNA <- cbind(meth_imputeNA, recovery_imputeNA)

meth_imputeNA <- cbind(meth_imputeNA, atrophy2_imputeNA)

saveRDS(meth_imputeNA, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/meth_imputeNA_mean.RDATA")

meth_imputeNA <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/meth_imputeNA_mean.RDATA")

#########################################

### plot PCA


pca.out <- prcomp(t(meth_imputeNA[2:37]), scale. = FALSE)

saveRDS(pca.out, "/Users/maxul/Documents/Skole/Lab/RMA_RRBS/pca.out.RDATA")

pca.out <- readRDS("/Users/maxul/Documents/Skole/Lab/RMA_RRBS/pca.out.RDATA")

# Extract the proportion of variance explained by each principal component
var_explained <- pca.out$sdev^2 / sum(pca.out$sdev^2)

# Create a data frame for plotting
df <- data.frame(Component = 1:length(var_explained), Variance = var_explained)

# Keep only the first 50 principal components
df <- df[1:50,]

# Plot

ggplot(df, aes(x = Component, y = Variance)) +
        geom_bar(stat = "identity") +
        scale_x_continuous(breaks = 1:50) +
        labs(x = "Principal Component", y = "Proportion of Variance Explained",
             title = "Variance Explained by Principal Components")+
        theme_bw()

library(ggrepel)

ggplot(data.frame(pca.out$x), aes(x = PC1, y = PC2)) +
        geom_point() +
        geom_text_repel(data = data.frame(pca.out$x), aes(label = rownames(data.frame(pca.out$x)))) +
        labs(x = "PC1", y = "PC2", title = "First Two Principal Components")+
        theme_bw()


# check for outliers

# Compute Mahalanobis distances for the first two principal components
distances <- mahalanobis(pca.out$x[,1:2], colMeans(pca.out$x[,1:2]), cov(pca.out$x[,1:2]))

# Identify outliers as samples with a Mahalanobis distance greater than a certain threshold
# Here, I'm using the 97.5 percentile of the Chi-square distribution with 3 degrees of freedom as the threshold
outliers <- which(distances > qchisq(0.975, df = 3))

# Print the row names of the outliers
print(colnames(filtered_df[2:37])[outliers])

# Plot the first two principal components, highlighting the outliers

ggplot(data.frame(pca.out$x), aes(x = PC1, y = PC2)) +
        geom_point() +
        geom_text_repel(data = data.frame(pca.out$x)[outliers,], aes(label = rownames(data.frame(pca.out$x))[outliers])) +
        labs(x = "PC1", y = "PC2", title = "First Two Principal Components with Outliers Highlighted")






# density plot

library(reshape2)

melt(meth_imputeNA[2:37],variable.name = "Sample", value.name = "Value") %>% 
        
        ggplot(aes(x = Value, color = Sample)) +
        geom_density() +
        theme_minimal() +
        labs(title = "ImputedNA Density Plot", x = "Percent_meth", y = "Density")



###########################################

### remove non-variable CpGs to greatly reduce total number of diff-methylation tests

# example code that needs to be adapted from https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationSeq/Seq_Tutorial.html

# get percent methylation matrix
pm=percMethylation(meth)

# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)

# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
hist(sds, breaks = 100)

# keep only CpG with standard deviations larger than 2%
meth <- meth[sds > 2]

# This leaves us with this number of CpG sites
nrow(meth)



#########################################

### remove known SNPs










##########################################

### run linear mixed effects model on dataframe with missing values

library(lme4)
library(emmeans)


colnames(filtered_df)



# add pheno data

FP = c(3,4,5,7,8,9,12,15,27)
sex = c("female","male","male","male","male","female","male","male","male")
age = c(21,21,30,26,33,27,25,31,28)

lmm_data <- data.frame("FP" = as.factor(FP),
                       "sex" = as.factor(sex),
                       "age" = age)

# transform dataframe to long format


lmm_data$baseline <- as.numeric(filtered_df[1,2:10])

lmm_data$atrphy1 <- as.numeric(filtered_df[1,11:19])

lmm_data$recovery <- as.numeric(filtered_df[1,20:28])

lmm_data$atrophy2 <- as.numeric(filtered_df[1,29:37])


# change to long format

lmm_data %>% 
        pivot_longer(cols = 4:7, names_to = "timepoint", values_to = "percent_meth") -> x

#x$timepoint = factor(x$timepoint, levels = c("baseline", "atrophy1", "recovery","atrophy2"))


# run mixed model

model <- lmer(percent_meth ~ timepoint + age + sex + (1|FP), data = x)

summary(model)


# post-hoc analysis for timepoint

post_hoc_results <- emmeans(model, pairwise ~ timepoint)

print(post_hoc_results)



plot(resid(model))
qqnorm(resid(model))
qqline(resid(model))



# run on all probes

lmm_models <- list()


for (i in 1:100) {
        
        
        lmm_data$baseline <- as.numeric(filtered_df[i,2:10])
        
        lmm_data$atrphy1 <- as.numeric(filtered_df[i,11:19])
        
        lmm_data$recovery <- as.numeric(filtered_df[i,20:28])
        
        lmm_data$atrophy2 <- as.numeric(filtered_df[i,29:37])
        
        # change to long format
        
        lmm_data %>% 
                pivot_longer(cols = 4:7, names_to = "timepoint", values_to = "percent_meth") -> x
        
        # run mixed model
        
        model <- lmer(percent_meth ~ timepoint + age + sex + (1|FP), data = x)
        
        
        lmm_models[filtered_df$ID[i]] <- model
        
        print(i)
        
        
        
}


emmeans_models <- list()

for (i in 1:length(lmm_models)) {
        
        try({
                
                post_hoc_results <- emmeans(lmm_models[[i]], pairwise ~ timepoint)
                
        }, silent = TRUE)
        
        if (any(as.data.frame(post_hoc_results[[2]])$p.value < 0.05)) {
                
                emmeans_models[[names(lmm_models[i])]] <- post_hoc_results
        }
        
        print(i)
        
        
}


summary(emmeans_models[[1]])

summary(lmm_models[["chr1_1000135"]])

plot(resid(lmm_models[["chr1_1000135"]]))
qqnorm(resid(lmm_models[["chr1_1000135"]]))
qqline(resid(lmm_models[["chr1_1000135"]]))

plot(emmeans_models[[1]])

































