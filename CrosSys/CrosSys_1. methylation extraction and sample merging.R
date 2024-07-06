#
#               CrosSys
#
#               Load and merge data
#

BiocManager::install('methylKit')

# Load packages
library(methylKit)
library(tidyverse)
library(ggplot2)


# Set working directory max's mac
setwd("/Users/maxullrich/OneDrive - UGent/CrosSys/")

# Set working directory NIH pc
setwd("D:/OneDrive - UGent/CrosSys/")

# load and combine all three files, filter for mincov = 10, and then save on harddisk. Remove file from R, and continue to next file pair. 



##########################################################################################################
#########               get all files                          ###########################################
##########################################################################################################


# Locate all files
meth_files_1 <- as.list(list.files(pattern = '*.cov.gz', 
                                      path = "./methylation_results/methylation_results_1st_seq", full.names = TRUE))

# remove nedativa and positive control
meth_files_1 <- meth_files_1[1:135]

meth_files_2 <- as.list(list.files(pattern = '*.cov.gz', 
                                   path = "./methylation_results/methylation_results_2nd_seq", full.names = TRUE))

meth_files <- c(meth_files_1, meth_files_2)


##########################################################################################################
#########               identify sample name                   ###########################################
##########################################################################################################

# get the list of unique file names

sample_names <- sub(".*220048-II-(.*?)_S.*", "\\1", unlist(meth_files)) %>% unique()



##########################################################################################################
#########               loop that:                                                      ##################
#########               1. identifies matching samples                                  ##################
#########               2. gets the data                                                ##################
#########               3. merges and filters for minimum coverage of 10                ##################
#########               4. saves the merged filtered file in folder "./merged_meth_files/"   ##################
##########################################################################################################

treatment <- c()

treatment[grep('PRE-ASAT', meth_files)]  <- 0
treatment[grep('POST-ASAT', meth_files)] <- 1
treatment[grep('PRE-QF', meth_files)]    <- 2
treatment[grep('POST-QF', meth_files)]   <- 3


for (i in 1:length(sample_names)) {
        
        # Check if the file already exists
        if (file.exists(paste0("./methylation_results/merged_meth_files/", sample_names[i], ".RDATA"))) {
                print(paste("File already exists for sample:", sample_names[i], "- skipping."))
                next
        }
        
        
        
        index_1 <- grep(sample_names[i], meth_files)
        
        # get the file paths
        
        sample_files <- c(meth_files[index_1])
        

        # Read in data using methRead
        sample_data <- methRead(sample_files,
                                sample.id = as.list(sample_files),
                                assembly = 'hg38',
                                pipeline = 'bismarkCoverage',
                                header = F,
                                treatment = treatment[index_1],
                                mincov = 1)
        
        print(paste("Read data for sample:", sample_names[i]))
        
        # merge the files together
        
        # Initialize merged_data with the first data set
        merged_data <- getData(sample_data[[1]])
        
        # Loop through the rest of the sample_data starting from the second file
        for (n in 2:length(sample_data)) {
                merged_data <- merge(merged_data,
                                     getData(sample_data[[n]]),
                                     by = c("chr", "start", "end", "strand"),
                                     all = TRUE)
                
        }
        
        print(paste("Merged data for sample:", sample_names[i]))
        
        
        # Replace NA with 0 in coverage, numCs, and numTs columns
        merged_data <- merged_data %>%
                mutate(across(colnames(merged_data)[5:length(colnames(merged_data))], ~if_else(is.na(.), 0, .))) %>%
                mutate(coverage = rowSums(select(., starts_with("coverage")), na.rm = TRUE),
                       numCs = rowSums(select(., starts_with("numCs")), na.rm = TRUE),
                       numTs = rowSums(select(., starts_with("numTs")), na.rm = TRUE)) %>%
                select(chr, start, end, strand, coverage, numCs, numTs) %>%
                filter(coverage > 9) %>% 
                as.data.frame()
        
        
        saveRDS(merged_data, paste0("./methylation_results/merged_meth_files/", sample_names[i], ".RDATA"))
        
        print(paste("Finished with file", i, "/", length(sample_names), " filename:", sample_names[i]))
        
}


