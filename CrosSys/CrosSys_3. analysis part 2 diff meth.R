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
                filter(ID %in% united_ID_qf) 
        
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





############################################################################################
############################################################################################
##########              run differential methylation analysis        #######################
############################################################################################
############################################################################################

# Load necessary libraries
library(limma)
library(ENmix)

# Extract relevant information
extract_info <- function(name) {
        parts <- strsplit(name, "-")[[1]]
        twin_id <- substr(parts[2], 1, 2)
        twin_status <- substr(parts[2], 3, 3)
        time_point <- tolower(parts[3])
        return(c(twin_id, twin_status, time_point))
}

###################################
###             ASAT            ###
###################################

# Example column names
column_names <- colnames(percent_meth_ASAT[,2:36])


info_matrix <- t(sapply(column_names, extract_info))
info_df <- as.data.frame(info_matrix, stringsAsFactors = FALSE)
colnames(info_df) <- c("twin_id", "twin_status", "time_point")

# Convert to appropriate types
info_df$twin_id <- as.factor(info_df$twin_id)
info_df$twin_status <- as.factor(info_df$twin_status)
info_df$time_point <- as.factor(info_df$time_point)

# Add a sample name column
info_df$sample_name <- column_names

# Ensure percent_meth_ASAT is a matrix of methylation values
percent_meth_matrix <- as.matrix(percent_meth_ASAT[, 2:36])

# Create the design matrix for the model
design <- model.matrix(~ twin_id + twin_status * time_point, data = info_df)

# View the design matrix
print(design)


# Adjust the design matrix with valid column names
colnames(design) <- make.names(colnames(design))

# View the design matrix
print(design)

# Fit the linear model
fit <- lmFit(percent_meth_matrix, design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Adjust the design matrix with valid column names
colnames(design) <- make.names(colnames(design))

contrast_matrix <- makeContrasts(
        # Pre-post difference within twin_status 2
        pre_post_diff_within_twin_status2 = (X.Intercept. + twin_status2) - (X.Intercept. + twin_status2 + time_pointpre + twin_status2.time_pointpre),
        
        # Pre-post difference within twin_status 1 (reference level)
        pre_post_diff_within_twin_status1 = X.Intercept. - (X.Intercept. + time_pointpre),
        
        # Difference in pre-post changes between twin statuses
        pre_post_diff_between_twin_statuses = ((X.Intercept. + twin_status2) - (X.Intercept. + twin_status2 + time_pointpre + twin_status2.time_pointpre)) - 
                (X.Intercept. - (X.Intercept. + time_pointpre)),
        levels = design
)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get the results for the contrasts
results_ASAT_pre_post_diff_within_twin_status2 <- topTable(fit2, adjust.method = "fdr", number = Inf,coef = 1, sort.by = "none")
results_ASAT_pre_post_diff_within_twin_status1 <- topTable(fit2, adjust.method = "fdr", number = Inf,coef = 2, sort.by = "none")
results_ASAT_pre_post_diff_between_twin_statuses <- topTable(fit2, adjust.method = "fdr", number = Inf,coef = 3, sort.by = "none")


#######################################################################################
##########              Check assumptions for ASAT                #######################
#######################################################################################

# Obtain residuals from the fit object
residuals <- residuals(fit, percent_meth_matrix)

# Convert residuals to a data frame for plotting
residuals_df <- as.data.frame(residuals)
colnames(residuals_df) <- colnames(percent_meth_matrix)

# Convert percent_meth_matrix and residuals_df to data frames and remove NAs
percent_meth_df <- as.data.frame(percent_meth_matrix) %>% na.omit()
residuals_df <- residuals_df %>% na.omit()

# Check the dimensions of the data frames
dim(percent_meth_df)
dim(residuals_df)

# Pivot the data frames longer for plotting
percent_meth_long <- pivot_longer(percent_meth_df , cols = everything(), names_to = "Sample", values_to = "percent_meth")
residuals_long <- pivot_longer(residuals_df , cols = everything(), names_to = "Sample", values_to = "residuals")

# Ensure both data frames are of the same length
if (nrow(percent_meth_long) == nrow(residuals_long)) {
        combined_df <- cbind(percent_meth_long, residuals = residuals_long$residuals)
        
        # Plot residuals vs percent methylation
        ggplot(combined_df, aes(x = percent_meth, y = residuals)) +
                geom_point() +
                geom_hline(yintercept = 0, col = "red") +
                ggtitle("Residuals vs Percent Methylation") +
                xlab("Percent Methylation") +
                ylab("Residuals") +
                theme_minimal()
} else {
        stop("The lengths of the percent methylation data and residuals do not match.")
}

###################################
###             QF              ###
###################################

# Example column names
column_names <- colnames(percent_meth_QF[,2:34])

info_matrix <- t(sapply(column_names, extract_info))
info_df <- as.data.frame(info_matrix, stringsAsFactors = FALSE)
colnames(info_df) <- c("twin_id", "twin_status", "time_point")

# Convert to appropriate types
info_df$twin_id <- as.factor(info_df$twin_id)
info_df$twin_status <- as.factor(info_df$twin_status)
info_df$time_point <- as.factor(info_df$time_point)

# Add a sample name column
info_df$sample_name <- column_names

# Ensure percent_meth_ASAT is a matrix of methylation values
percent_meth_matrix <- as.matrix(percent_meth_QF[, 2:34])

# Create the design matrix for the model
design <- model.matrix(~ twin_id + twin_status * time_point, data = info_df)

# View the design matrix
print(design)


# Adjust the design matrix with valid column names
colnames(design) <- make.names(colnames(design))

# View the design matrix
print(design)

# Fit the linear model
fit <- lmFit(percent_meth_matrix, design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Adjust the design matrix with valid column names
colnames(design) <- make.names(colnames(design))


contrast_matrix <- makeContrasts(
        # Pre-post difference within twin_status 2
        pre_post_diff_within_twin_status2 = (X.Intercept. + twin_status2) - (X.Intercept. + twin_status2 + time_pointpre + twin_status2.time_pointpre),
        
        # Pre-post difference within twin_status 1 (reference level)
        pre_post_diff_within_twin_status1 = X.Intercept. - (X.Intercept. + time_pointpre),
        
        # Difference in pre-post changes between twin statuses
        pre_post_diff_between_twin_statuses = ((X.Intercept. + twin_status2) - (X.Intercept. + twin_status2 + time_pointpre + twin_status2.time_pointpre)) - 
                (X.Intercept. - (X.Intercept. + time_pointpre)),
        levels = design
)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

colnames(coef(fit2))

# Get the results for the contrasts
results_QF_pre_post_diff_within_twin_status2 <- topTable(fit2, adjust.method = "fdr", number = Inf,coef = "pre_post_diff_within_twin_status2", sort.by = "none") 
results_QF_pre_post_diff_within_twin_status1 <- topTable(fit2, adjust.method = "fdr", number = Inf,coef = "pre_post_diff_within_twin_status1", sort.by = "none")
results_QF_pre_post_diff_between_twin_statuses <- topTable(fit2, adjust.method = "fdr", number = Inf,coef = "pre_post_diff_between_twin_statuses", sort.by = "none") 


#######################################################################################
##########              Check assumptions for QF                #######################
#######################################################################################

# Obtain residuals from the fit object
residuals <- residuals(fit, percent_meth_matrix)

# Convert residuals to a data frame for plotting
residuals_df <- as.data.frame(residuals)
colnames(residuals_df) <- colnames(percent_meth_matrix)

# Convert percent_meth_matrix and residuals_df to data frames and remove NAs
percent_meth_df <- as.data.frame(percent_meth_matrix) %>% na.omit()
residuals_df <- residuals_df %>% na.omit()

# Check the dimensions of the data frames
dim(percent_meth_df)
dim(residuals_df)

# Pivot the data frames longer for plotting
percent_meth_long <- pivot_longer(percent_meth_df , cols = everything(), names_to = "Sample", values_to = "percent_meth")
residuals_long <- pivot_longer(residuals_df , cols = everything(), names_to = "Sample", values_to = "residuals")

# Ensure both data frames are of the same length
if (nrow(percent_meth_long) == nrow(residuals_long)) {
        combined_df <- cbind(percent_meth_long, residuals = residuals_long$residuals)
        
        # Plot residuals vs percent methylation
        ggplot(combined_df, aes(x = percent_meth, y = residuals)) +
                geom_point() +
                geom_hline(yintercept = 0, col = "red") +
                ggtitle("Residuals vs Percent Methylation") +
                xlab("Percent Methylation") +
                ylab("Residuals") +
                theme_minimal()
} else {
        stop("The lengths of the percent methylation data and residuals do not match.")
}


#######################################################################################
##########              save diff meth before annotation        #######################
#######################################################################################

write.csv(results_ASAT_pre_post_diff_between_twin_statuses, "./methylation_results/results_ASAT_pre_post_diff_between_twin_statuses.csv")
write.csv(results_ASAT_pre_post_diff_within_twin_status1, "./methylation_results/results_ASAT_pre_post_diff_within_twin_status1.csv")
write.csv(results_ASAT_pre_post_diff_within_twin_status2, "./methylation_results/results_ASAT_pre_post_diff_within_twin_status2.csv")

write.csv(results_QF_pre_post_diff_between_twin_statuses, "./methylation_results/results_QF_pre_post_diff_between_twin_statuses.csv")
write.csv(results_QF_pre_post_diff_within_twin_status1, "./methylation_results/results_QF_pre_post_diff_within_twin_status1.csv")
write.csv(results_QF_pre_post_diff_within_twin_status2, "./methylation_results/results_QF_pre_post_diff_within_twin_status2.csv")


#######################################################################################
##########              Annotation                              #######################
#######################################################################################

library(annotatr)

annots = c('hg38_basicgenes')

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg38', annotations = annots)

annotation_df <- as.data.frame(annotations) 

i =1 

for (i in 1:length(unique(annotation_df$symbol))) {
        
        annotation_df %>% 
                filter(symbol == unique(annotation_df$symbol)[i])
        
        
}


percent_meth_ASAT_ID <- percent_meth_ASAT %>% 
        as.data.frame() %>% 
        dplyr::select(ID) %>% 
        rowwise() %>% 
        mutate(chr = str_split(ID, pattern = "\\.")[[1]][1],
               start = str_split(ID, pattern = "\\.")[[1]][2],
               end = (str_split(ID, pattern = "\\.")[[1]][3])) %>% 
        ungroup() %>% 
        dplyr::select(chr, start, end)

percent_meth_ASAT_ID$start <- as.numeric(percent_meth_ASAT_ID$start)
percent_meth_ASAT_ID$end <- as.numeric(percent_meth_ASAT_ID$end)

# Convert methylation data frame to GRanges object
meth_ASAT_gr <- makeGRangesFromDataFrame(percent_meth_ASAT_ID, 
                                    keep.extra.columns = TRUE, 
                                    seqnames.field = "chr", 
                                    start.field = "start", 
                                    end.field = "end")

# Check sequence levels in methylation data
seqlevels(meth_ASAT_gr)

# Check sequence levels in annotation data
seqlevels(annotations)

# Remove "chr" prefix if needed
seqlevelsStyle(meth_ASAT_gr) <- "NCBI"  # This will set the style to NCBI (e.g., "1")
seqlevelsStyle(annotations) <- "NCBI"   # Ensure the annotation data also uses UCSC style


ASAT_annotations <- GenomicRanges::findOverlaps(meth_ASAT_gr, annotations)

# Extract the annotation data for the overlapping ranges
annotation_hits <- annotations[subjectHits(ASAT_annotations)]
methylation_hits <- meth_ASAT_gr[queryHits(ASAT_annotations)]

# Combine methylation data with annotation
annotated_meth_df <- as.data.frame(mcols(methylation_hits))
annotation_data <- as.data.frame(mcols(annotation_hits))

# Combine methylation and annotation data
annotated_meth_df <- cbind(as.data.frame(methylation_hits), annotation_data)

annotated_ASAT_df <- annotated_meth_df %>% 
        mutate(ID = paste(seqnames, start, end, sep = ".")) %>% 
        distinct(ID, .keep_all = TRUE) 


###########
#       annotate QF
####


percent_meth_QF_ID <- percent_meth_QF %>% 
        as.data.frame() %>% 
        dplyr::select(ID) %>% 
        rowwise() %>% 
        mutate(chr = str_split(ID, pattern = "\\.")[[1]][1],
               start = str_split(ID, pattern = "\\.")[[1]][2],
               end = (str_split(ID, pattern = "\\.")[[1]][3])) %>% 
        ungroup() %>% 
        dplyr::select(chr, start, end)

percent_meth_QF_ID$start <- as.numeric(percent_meth_QF_ID$start)
percent_meth_QF_ID$end <- as.numeric(percent_meth_QF_ID$end)

# Convert methylation data frame to GRanges object
meth_QF_gr <- makeGRangesFromDataFrame(percent_meth_QF_ID, 
                                         keep.extra.columns = TRUE, 
                                         seqnames.field = "chr", 
                                         start.field = "start", 
                                         end.field = "end")

# Check sequence levels in methylation data
seqlevels(meth_QF_gr)

# Check sequence levels in annotation data
seqlevels(annotations)

# Remove "chr" prefix if needed
seqlevelsStyle(meth_QF_gr) <- "NCBI"  # This will set the style to NCBI (e.g., "1")
seqlevelsStyle(annotations) <- "NCBI"   # Ensure the annotation data also uses UCSC style


QF_annotations <- GenomicRanges::findOverlaps(meth_QF_gr, annotations)

# Extract the annotation data for the overlapping ranges
annotation_hits <- annotations[subjectHits(QF_annotations)]
methylation_hits <- meth_QF_gr[queryHits(QF_annotations)]

# Combine methylation data with annotation
annotated_meth_df <- as.data.frame(mcols(methylation_hits))
annotation_data <- as.data.frame(mcols(annotation_hits))

# Combine methylation and annotation data
annotated_meth_df <- cbind(as.data.frame(methylation_hits), annotation_data)

annotated_QF_df <- annotated_meth_df %>% 
        mutate(ID = paste(seqnames, start, end, sep = ".")) %>% 
        distinct(ID, .keep_all = TRUE) 


######
#       merge annotations together
#####

annotated_df <- rbind(annotated_ASAT_df, annotated_QF_df) %>% 
        distinct(ID, .keep_all = TRUE)

# add ID to meth analysis and merge with annotation

results_ASAT_pre_post_diff_between_twin_statuses$ID <- percent_meth_ASAT$ID
results_ASAT_pre_post_diff_within_twin_status1$ID <- percent_meth_ASAT$ID
results_ASAT_pre_post_diff_within_twin_status2$ID <- percent_meth_ASAT$ID

results_QF_pre_post_diff_between_twin_statuses$ID <- percent_meth_QF$ID
results_QF_pre_post_diff_within_twin_status1$ID <- percent_meth_QF$ID
results_QF_pre_post_diff_within_twin_status2$ID <- percent_meth_QF$ID

results_ASAT_pre_post_diff_between_twin_statuses <- left_join(results_ASAT_pre_post_diff_between_twin_statuses, annotated_df, by = "ID")
results_ASAT_pre_post_diff_within_twin_status1 <- left_join(results_ASAT_pre_post_diff_within_twin_status1, annotated_df, by = "ID")
results_ASAT_pre_post_diff_within_twin_status2 <- left_join(results_ASAT_pre_post_diff_within_twin_status2, annotated_df, by = "ID")

results_QF_pre_post_diff_between_twin_statuses <- left_join(results_QF_pre_post_diff_between_twin_statuses, annotated_df, by = "ID")
results_QF_pre_post_diff_within_twin_status1 <- left_join(results_QF_pre_post_diff_within_twin_status1, annotated_df, by = "ID")
results_QF_pre_post_diff_within_twin_status2 <- left_join(results_QF_pre_post_diff_within_twin_status2, annotated_df, by = "ID")


#######################################################################################
##########              filter diff meth                        #######################
#######################################################################################


results_ASAT_pre_post_diff_between_twin_statuses_flt <- results_ASAT_pre_post_diff_between_twin_statuses %>% 
        filter(P.Value < 0.05) 

results_ASAT_pre_post_diff_within_twin_status1_flt <- results_ASAT_pre_post_diff_within_twin_status1 %>% 
        filter(P.Value < 0.05) 

results_ASAT_pre_post_diff_within_twin_status2_flt <- results_ASAT_pre_post_diff_within_twin_status2 %>% 
        filter(P.Value < 0.05) 


results_QF_pre_post_diff_between_twin_statuses_flt <- results_QF_pre_post_diff_between_twin_statuses %>% 
        filter(P.Value < 0.05) 

results_QF_pre_post_diff_within_twin_status1_flt <- results_QF_pre_post_diff_within_twin_status1 %>% 
        filter(P.Value < 0.05) 

results_QF_pre_post_diff_within_twin_status2_flt <- results_QF_pre_post_diff_within_twin_status2 %>% 
        filter(P.Value < 0.05) 



#######################################################################################
##########              save diff meth                          #######################
#######################################################################################

write.csv(results_ASAT_pre_post_diff_between_twin_statuses_flt, "./methylation_results/results_ASAT_pre_post_diff_between_twin_statuses_flt.csv")
write.csv(results_ASAT_pre_post_diff_within_twin_status1_flt, "./methylation_results/results_ASAT_pre_post_diff_within_twin_status1_flt.csv")
write.csv(results_ASAT_pre_post_diff_within_twin_status2_flt, "./methylation_results/results_ASAT_pre_post_diff_within_twin_status2_flt.csv")

write.csv(results_QF_pre_post_diff_between_twin_statuses_flt, "./methylation_results/results_QF_pre_post_diff_between_twin_statuses_flt.csv")
write.csv(results_QF_pre_post_diff_within_twin_status1_flt, "./methylation_results/results_QF_pre_post_diff_within_twin_status1_flt.csv")
write.csv(results_QF_pre_post_diff_within_twin_status2_flt, "./methylation_results/results_QF_pre_post_diff_within_twin_status2_flt.csv")


