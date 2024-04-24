#
#
#    CrosSys RRBS data
#
#    quality control script


# This script is specific for CrosSys, data originating from Finnish group



library(tidyverse)


# set working directory to the mapping results folder

setwd("/Users/maxul/Documents/Skole/Lab/CrosSys/mapping_results2/")



# get list of reposrts in the folder

reports <- list.files(pattern = "PE_report.txt$")


# extract all the data

reports_list <- list()

for (i in 1:length(reports)) {
        
        x <- read.delim(reports[i], sep = ":", header = FALSE)
        
        colnames(x) <- c("output", sub("_R[12].*", "", reports[i]))
        reports_list[[i]] <- x
        
        names(reports_list)[i] <- sub("_R[12].*", "", reports[i])
        
        print(i)
        
}


# extract mapping efficiency and 


reports_df <- reports_list[[1]][6:8,1:2] 

for (i in 2:length(reports_list)) {
        
        
        reports_df <- left_join(reports_df, reports_list[[i]][6:8,1:2], by = "output")
        
        print(i)
}


# plot mapping efficiency and sequencing depth

# plot ASAT


mapping_efficiency <-  reports_df %>%
        pivot_longer(names_to = "Sample", values_to = "val", cols = 2:(length(colnames(reports_df)))) %>%
        mutate(
                val = gsub("\t", "", val),
                val = as.numeric(gsub("%", "", val)),
                output = factor(output, levels = c("Sequence pairs analysed in total",
                                                   "Number of paired-end alignments with a unique best hit",
                                                   "Mapping efficiency"))) %>%
        filter(output == "Mapping efficiency") %>% 
        mutate(Tissue = ifelse(grepl("ASAT", Sample), "ASAT", "QF")) %>% 
        filter(Tissue == "ASAT")
        


reports_df %>%
        pivot_longer(names_to = "Sample", values_to = "val", cols = 2:(length(colnames(reports_df)))) %>%
        mutate(
                val = gsub("\t", "", val),
                val = as.numeric(gsub("%", "", val)),
                output = factor(output, levels = c("Sequence pairs analysed in total",
                                                   "Number of paired-end alignments with a unique best hit",
                                                   "Mapping efficiency"))) %>%
        filter(output != "Mapping efficiency") %>% 
        mutate(Tissue = ifelse(grepl("ASAT", Sample), "ASAT", "QF")) %>% 
        filter(Tissue == "ASAT") %>% 
        ggplot(aes(x = Sample, y = val, fill = output)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.6) +
        scale_y_continuous(n.breaks = 10) +
        scale_fill_manual(values = c("red","#440154FF", "#5DC863FF")) +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              axis.title.x = element_blank()) +
        geom_text(data = mapping_efficiency, aes(label = sprintf("%.1f%%", val), y = 12000000),
                  position = position_dodge(width = 0.6),
                  size = 5, color = "red",
                  angle = 90)+
        theme(legend.position = "top",
              legend.title = element_blank())+
        labs(y = "Sequences")


# plot GF


mapping_efficiency <-  reports_df %>%
        pivot_longer(names_to = "Sample", values_to = "val", cols = 2:(length(colnames(reports_df)))) %>%
        mutate(
                val = gsub("\t", "", val),
                val = as.numeric(gsub("%", "", val)),
                output = factor(output, levels = c("Sequence pairs analysed in total",
                                                   "Number of paired-end alignments with a unique best hit",
                                                   "Mapping efficiency"))) %>%
        filter(output == "Mapping efficiency") %>% 
        mutate(Tissue = ifelse(grepl("ASAT", Sample), "ASAT", "QF")) %>% 
        filter(Tissue == "QF")



qf_plot <- reports_df %>%
        pivot_longer(names_to = "Sample", values_to = "val", cols = 2:(length(colnames(reports_df)))) %>%
        mutate(
                val = gsub("\t", "", val),
                val = as.numeric(gsub("%", "", val)),
                output = factor(output, levels = c("Sequence pairs analysed in total",
                                                   "Number of paired-end alignments with a unique best hit",
                                                   "Mapping efficiency"))) %>%
        filter(output != "Mapping efficiency") %>% 
        mutate(Tissue = ifelse(grepl("ASAT", Sample), "ASAT", "QF")) %>% 
        filter(Tissue == "QF") %>% 
        ggplot(aes(x = Sample, y = val, fill = output)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.6) +
        scale_y_continuous(n.breaks = 10) +
        scale_fill_manual(values = c("red","#440154FF", "#5DC863FF")) +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              axis.title.x = element_blank()) +
        geom_text(data = mapping_efficiency, aes(label = sprintf("%.1f%%", val), y = 12000000),
                  position = position_dodge(width = 0.6),
                  size = 5, color = "red",
                  angle = 90)+
        theme(legend.position = "top",
              legend.title = element_blank())+
        labs(y = "Sequences")



# save
ggsave("/Users/maxul/Documents/Skole/Lab/CrosSys/QF_QC.emf", plot = qf_plot, width = 12, height = 8)

##############################################################################3
####### time of allignment     ###############################################
################################################################################3

# Regular expression to extract hours, minutes, and seconds
pattern <- "(\\d+)h (\\d+)m (\\d+)s"

# Use str_extract to pull out the time components
library(stringr)

time <- data.frame("Time" = NA, "Hours"= NA, "Minutes"= NA, "Seconds"= NA)

for (i in 1:length(reports_list)) {
        x <- as.data.frame(str_match(reports_list[[i]][33,1], pattern)) 
        colnames(x) <- c("Time", "Hours", "Minutes", "Seconds")
        time[i,] <- x
        
}

# total time in hours spent on bismark allign (not including trimming or methylation extraction)
(sum(as.numeric(time$Seconds))/60 + sum(as.numeric(time$Minutes)))/60



