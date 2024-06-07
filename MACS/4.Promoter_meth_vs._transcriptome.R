###
###
###     Baseline gene expression for MACS comparison
###
###






BiocManager::install("biomaRt")
library(biomaRt)

library(DESeq2)
library(tidyverse)
library(EpiSCORE)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(cowplot)



# Select the dataset corresponding to your organism, e.g., Homo sapiens
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene information including gene length
genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position"), mart = ensembl)

# Calculate gene lengths
genes$length <- genes$end_position - genes$start_position + 1

genes


norm_counts <- counts(dds, normalize = TRUE)


metadata <- colData(dds)

###
#   get only the baseline samples
###

pre_smaples <- metadata %>%
    as.data.frame() %>%
    filter(Time == "Pre") %>%
    pull(RNA_number)

norm_baseline_counts <- norm_counts[, colnames(norm_counts) %in% pre_smaples]

###
#   filter data so that i have gene length for all genes
###



baseline_genes <- rownames(norm_baseline_counts) %>% length()

genes <- genes %>%
    as.data.frame() %>%
    filter(external_gene_name %in% baseline_genes)%>%
    distinct(external_gene_name, .keep_all = TRUE)

filt_norm_baseline_counts <- norm_baseline_counts[rownames(norm_baseline_counts)%in% genes$external_gene_name,]

df <- filt_norm_baseline_counts %>%
    as.data.frame() %>%
    rownames_to_column(var = "external_gene_name") %>%
    merge(., as.data.frame(genes), by = "external_gene_name")



# Calculate Counts Per Kilobase (CPK)
countsPerKb <- sweep(df[,2:12], 1, df[,18] / 1000, FUN="/")

# Sum CPKs for each sample to get total counts per sample
totalCountsPerSample <- colSums(countsPerKb)

# Calculate TPM
tpm <- sweep(countsPerKb, 2, totalCountsPerSample, FUN="/") * 1e6

# add genes as rownames

rownames(tpm) <- rownames(filt_norm_baseline_counts)

tpm$rowMean <- rowSums(tpm)/ncol(tpm)


# find the most expressed genes

tpm %>%
    arrange(-rowMean)








######################################################################################################
###         Get the average promoter methylation from MACS samples          ##########################
######################################################################################################


# DNA methylation profile
beta <- readRDS("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/Skole/M.Sc/Master 21-22/Master/DATA/Epigenetics/beta.RDATA")



# average the promoter methylation in my samples
average_promoter.m <- constAvBetaTSS(beta.m = beta, type = "850k")


average_promoter.m <- average_promoter.m %>%
    as.data.frame() %>%
    dplyr::select("1BH","2BH","4BH","5BH","6BH","7BH","8BH","12BH",
                  "1BM","2BM","4BM","5BM","6BM","7BM","8BM","12BM",
                  "1PH","2PH","4PH","5PH","6PH","7PH","8PH","12PH",
                  "1PM","2PM","4PM","5PM","6PM","7PM","8PM","12PM") %>%
    mutate(BH = (`1BH`+`2BH`+`4BH`+`5BH`+`6BH`+`7BH`+`8BH`+`12BH`)/8,
           BM = (`1BM`+`2BM`+`4BM`+`5BM`+`6BM`+`7BM`+`8BM`+`12BM`)/8,
           PH = (`1PH`+`2PH`+`4PH`+`5PH`+`6PH`+`7PH`+`8PH`+`12PH`)/8,
           PM = (`1PM`+`2PM`+`4PM`+`5PM`+`6PM`+`7PM`+`8PM`+`12PM`)/8)


# annotate with SYMBOL

gene_symbols <- select(org.Hs.eg.db, keys = rownames(average_promoter.m), columns = "SYMBOL", keytype = "ENTREZID")


# add gene names and save file

BH_average_promoter.m <- average_promoter.m %>%
    rownames_to_column(var = "ENTREZID") %>%
    merge(., gene_symbols, by = "ENTREZID") %>%
    dplyr::select("external_gene_name" = SYMBOL, BH)


BM_average_promoter.m <- average_promoter.m %>%
    rownames_to_column(var = "ENTREZID") %>%
    merge(., gene_symbols, by = "ENTREZID") %>%
    dplyr::select("external_gene_name" = SYMBOL, BM)





######################################################################################################
###         Merge and plot bulk transcriptome and methylome.                ##########################
######################################################################################################

# BH vs. bulk transcriptome

p1 <- tpm %>%
    rownames_to_column(var = "external_gene_name") %>%
    dplyr::select(external_gene_name, "TPM" = rowMean) %>%
    merge(., BH_average_promoter.m, by = "external_gene_name") %>%
    mutate(external_gene_name = ifelse(TPM < 1000, "", external_gene_name)) %>%
    arrange(-TPM) %>%
    head(1000) %>%
    filter(external_gene_name != "ACTA1") %>%
    ggplot(aes(x = BH, y = TPM, color = TPM))+
    geom_point(size = 2) +
    scale_y_sqrt(n.breaks = 6)+
    scale_color_viridis_b()+
    theme_classic(base_size = 14)+
    geom_label_repel(aes(label = external_gene_name), box.padding = 0.4, max.overlaps = 10, force = 3, point.padding = 1, color = "black")+
    labs(y = "TPM (Baseline Bulk Transcriptome)",
         x = "Baseline MYO+INT")+
    theme(legend.position = "none")


# BM vs. bulk transcripome

p2 <- tpm %>%
    rownames_to_column(var = "external_gene_name") %>%
    dplyr::select(external_gene_name, "TPM" = rowMean) %>%
    merge(., BM_average_promoter.m, by = "external_gene_name") %>%
    mutate(external_gene_name = ifelse(TPM < 1000, "", external_gene_name)) %>%
    arrange(-TPM) %>%
    head(1000) %>%
    filter(external_gene_name != "ACTA1") %>%
    ggplot(aes(x = BM, y = TPM, color = TPM))+
    geom_point(size = 2) +
    scale_y_sqrt(n.breaks = 6)+
    scale_color_viridis_b()+
    theme_classic(base_size = 14)+
    geom_label_repel(aes(label = external_gene_name), box.padding = 0.4, max.overlaps = 10, force = 3, point.padding = 1, color = "black")+
    labs(y = "TPM (Baseline Bulk Transcriptome)",
         x = "Baseline MYO")+
    theme(legend.position = "none")


######################################################################################################
###         get single fiber TPM.                                           ##########################
######################################################################################################

sf_counts <- readRDS("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/Ph.D/deconvolution/data/singlefiber_counts.RDATA")

sf_metadata <- readRDS("/Users/maxullrich/Library/CloudStorage/OneDrive-UGent/Ph.D/deconvolution/data/singlefiber_metadata.RDATA")


# get only baseline single fibers

pre_samples <- sf_metadata %>%
    filter(time == "Pre") %>%
    rownames_to_column(var = "sample") %>%
    pull(sample)


sf_counts_baseline <- sf_counts[,colnames(sf_counts) %in% pre_samples]

# get gene length for all genes in sf data

# Select the dataset corresponding to your organism, e.g., Homo sapiens
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene information including gene length
genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position"), mart = ensembl)

# Calculate gene lengths
genes$length <- genes$end_position - genes$start_position + 1



sf_baseline_genes <- rownames(sf_counts_baseline)

genes <- genes %>%
    as.data.frame() %>%
    filter(external_gene_name %in% sf_baseline_genes)%>%
    distinct(external_gene_name, .keep_all = TRUE)

sf_counts_baseline_flt <- sf_counts_baseline[rownames(sf_counts_baseline)%in% genes$external_gene_name,]

sf_df <- sf_counts_baseline_flt %>%
    as.data.frame() %>%
    rownames_to_column(var = "external_gene_name") %>%
    merge(., as.data.frame(genes), by = "external_gene_name")



# Calculate Counts Per Kilobase (CPK)
countsPerKb <- sweep(sf_df[,2:304], 1, sf_df[,309] / 1000, FUN="/")

# Sum CPKs for each sample to get total counts per sample
totalCountsPerSample <- colSums(countsPerKb)

# Calculate TPM
sf_tpm <- sweep(countsPerKb, 2, totalCountsPerSample, FUN="/") * 1e6

# add genes as rownames

rownames(sf_tpm) <- rownames(sf_counts_baseline_flt)

sf_tpm$rowMean <- rowSums(sf_tpm)/ncol(sf_tpm)


# find the most expressed genes

sf_tpm %>%
    arrange(-rowMean)


######################################################################################################
###         Merge and plot sf transcriptome and methylome.                  ##########################
######################################################################################################


# BH vs. sf transcriptome

p3 <- sf_tpm %>%
    rownames_to_column(var = "external_gene_name") %>%
    dplyr::select(external_gene_name, "TPM" = rowMean) %>%
    merge(., BH_average_promoter.m, by = "external_gene_name") %>%
    mutate(external_gene_name = ifelse(TPM < 1000, "", external_gene_name)) %>%
    arrange(-TPM) %>%
    head(1000) %>%
    filter(external_gene_name != "HP1BP3") %>%
    ggplot(aes(x = BH, y = TPM, color = TPM))+
    geom_point(size = 2) +
    scale_y_sqrt(n.breaks = 6)+
    scale_color_viridis_b()+
    theme_classic(base_size = 14)+
    geom_label_repel(aes(label = external_gene_name), box.padding = 0.4, max.overlaps = 10, force = 3, point.padding = 1, color = "black")+
    labs(y = "TPM (Baseline Single Fiber Transcriptome)",
         x = "Baseline MYO+INT")+
    theme(legend.position = "none")

# BM vs. sf transcripome

p4 <- sf_tpm %>%
    rownames_to_column(var = "external_gene_name") %>%
    dplyr::select(external_gene_name, "TPM" = rowMean) %>%
    merge(., BM_average_promoter.m, by = "external_gene_name") %>%
    mutate(external_gene_name = ifelse(TPM < 1000, "", external_gene_name)) %>%
    arrange(-TPM) %>%
    head(1000) %>%
    filter(external_gene_name != "HP1BP3") %>%
    ggplot(aes(x = BM, y = TPM, color = TPM))+
    geom_point(size = 2) +
    scale_y_sqrt(n.breaks = 6)+
    theme_classic(base_size = 14)+
    scale_color_viridis_b()+
    geom_label_repel(aes(label = external_gene_name), box.padding = 0.4, max.overlaps = 10, force = 3, point.padding = 1, colour = "black")+
    labs(y = "TPM (Baseline Single Fiber Transcriptome)",
         x = "Baseline MYO")+
    theme(legend.position = "none")




tpm["HP1BP3",]

sf_tpm["HP1BP3",]

sf_counts["MIR29A",]
sf_counts["MIR29B1",]

norm_counts["MIR29A",]
norm_counts["MIR29B1",]


# plot together

plot_grid(p1, p2, p3, p4, ncol = 2, labels = "AUTO", label_size = 20)


