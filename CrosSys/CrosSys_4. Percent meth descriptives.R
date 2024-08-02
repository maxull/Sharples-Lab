#
#               CrosSys
#
#               descriptive methylation stats
#



library(methylKit)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggvenn)

# Set working directory max's mac
setwd("/Users/maxullrich/OneDrive - UGent/CrosSys")

# Set working directory NIH pc
setwd("D:/OneDrive - UGent/CrosSys/")

##########################################################################################################
#########               get all files                          ###########################################
##########################################################################################################

percent_meth <- readRDS("methylation_results/percent_meth.RDATA")

results_ASAT_pre_post_diff_between_twin_statuses_flt <- read.csv("./methylation_results/results_ASAT_pre_post_diff_between_twin_statuses_flt.csv")
results_ASAT_pre_post_diff_within_twin_status1_flt   <- read.csv("./methylation_results/results_ASAT_pre_post_diff_within_twin_status1_flt.csv")
results_ASAT_pre_post_diff_within_twin_status2_flt   <- read.csv("./methylation_results/results_ASAT_pre_post_diff_within_twin_status2_flt.csv")

results_QF_pre_post_diff_between_twin_statuses_flt   <- read.csv("./methylation_results/results_QF_pre_post_diff_between_twin_statuses_flt.csv")
results_QF_pre_post_diff_within_twin_status1_flt     <- read.csv("./methylation_results/results_QF_pre_post_diff_within_twin_status1_flt.csv")
results_QF_pre_post_diff_within_twin_status2_flt     <- read.csv("./methylation_results/results_QF_pre_post_diff_within_twin_status2_flt.csv")



##########################################################################################################
#########               hypo vs. hyper methyaltion                          ##############################
##########################################################################################################

ASAT_1 <- results_ASAT_pre_post_diff_within_twin_status1_flt %>% 
        mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>% 
        count(meth_direction) %>%
        mutate(percentage = n / sum(n) * 100) %>%
        ggplot(aes(x = meth_direction, fill = meth_direction, y = n)) +
        geom_bar(stat = "identity") +
        geom_text(aes(y = n / 2, label = paste0(n, "\n", round(percentage, 1), "%")), 
                  color = "white", 
                  fontface = "bold", 
                  size = 15) +
        theme_classic(base_size = 20) + 
        scale_fill_manual(values = c("#F5C242", "#353F4F")) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(axis.title.x = element_blank(),
              legend.position = "none") +
        labs(y = "Number of DMPs", 
             title = "ASAT exercise response", 
             subtitle = "Diff methylation heavier twin")

ASAT_2 <- results_ASAT_pre_post_diff_within_twin_status2_flt %>% 
        mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>% 
        count(meth_direction) %>%
        mutate(percentage = n / sum(n) * 100) %>%
        ggplot(aes(x = meth_direction, fill = meth_direction, y = n)) +
        geom_bar(stat = "identity") +
        geom_text(aes(y = n / 2, label = paste0(n, "\n", round(percentage, 1), "%")), 
                  color = "white", 
                  fontface = "bold", 
                  size = 15) +
        theme_classic(base_size = 20) + 
        scale_fill_manual(values = c("#F5C242", "#353F4F")) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(axis.title.x = element_blank(),
              legend.position = "none") +
        labs(y = "Number of DMPs", 
             title = "ASAT exercise response", 
             subtitle = "Diff methylation lean twin")


ASAT_3 <- results_ASAT_pre_post_diff_between_twin_statuses_flt %>% 
        mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>% 
        count(meth_direction) %>%
        mutate(percentage = n / sum(n) * 100) %>%
        ggplot(aes(x = meth_direction, fill = meth_direction, y = n)) +
        geom_bar(stat = "identity") +
        geom_text(aes(y = n / 2, label = paste0(n, "\n", round(percentage, 1), "%")), 
                  color = "white", 
                  fontface = "bold", 
                  size = 15) +
        theme_classic(base_size = 20) + 
        scale_fill_manual(values = c("#F5C242", "#353F4F")) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(axis.title.x = element_blank(),
              legend.position = "none") +
        labs(y = "Number of DMPs", 
             title = "ASAT exercise response", 
             subtitle = "Diff methylation lean vs. heavy twin")

plot_grid(ASAT_1, ASAT_2, ASAT_3, nrow = 1)


QF_1 <- results_QF_pre_post_diff_within_twin_status1_flt %>% 
        mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>% 
        count(meth_direction) %>%
        mutate(percentage = n / sum(n) * 100) %>%
        ggplot(aes(x = meth_direction, fill = meth_direction, y = n)) +
        geom_bar(stat = "identity") +
        geom_text(aes(y = n / 2, label = paste0(n, "\n", round(percentage, 1), "%")), 
                  color = "white", 
                  fontface = "bold", 
                  size = 15) +
        theme_classic(base_size = 20) + 
        scale_fill_manual(values = c("#F5C242", "#353F4F")) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(axis.title.x = element_blank(),
              legend.position = "none") +
        labs(y = "Number of DMPs", 
             title = "QF exercise response", 
             subtitle = "Diff methylation heavier twin")

QF_2 <- results_QF_pre_post_diff_within_twin_status2_flt %>% 
        mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>% 
        count(meth_direction) %>%
        mutate(percentage = n / sum(n) * 100) %>%
        ggplot(aes(x = meth_direction, fill = meth_direction, y = n)) +
        geom_bar(stat = "identity") +
        geom_text(aes(y = n / 2, label = paste0(n, "\n", round(percentage, 1), "%")), 
                  color = "white", 
                  fontface = "bold", 
                  size = 15) +
        theme_classic(base_size = 20) + 
        scale_fill_manual(values = c("#F5C242", "#353F4F")) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(axis.title.x = element_blank(),
              legend.position = "none") +
        labs(y = "Number of DMPs", 
             title = "QF exercise response", 
             subtitle = "Diff methylation lean twin")


QF_3 <- results_QF_pre_post_diff_between_twin_statuses_flt %>% 
        mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>% 
        count(meth_direction) %>%
        mutate(percentage = n / sum(n) * 100) %>%
        ggplot(aes(x = meth_direction, fill = meth_direction, y = n)) +
        geom_bar(stat = "identity") +
        geom_text(aes(y = n / 2, label = paste0(n, "\n", round(percentage, 1), "%")), 
                  color = "white", 
                  fontface = "bold", 
                  size = 15) +
        theme_classic(base_size = 20) + 
        scale_fill_manual(values = c("#F5C242", "#353F4F")) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(axis.title.x = element_blank(),
              legend.position = "none") +
        labs(y = "Number of DMPs", 
             title = "QF exercise response", 
             subtitle = "Diff methylation lean vs. heavy twin")

plot_grid(QF_1, QF_2, QF_3, nrow = 1)


##########################################################################################################
#########               venn diagrams of hypo and hyper meth probes                        ###############
##########################################################################################################

####
#       HYPO
####

ASAT_HYPO_DMPs <- list(Lean_HYPO = results_ASAT_pre_post_diff_within_twin_status2_flt %>%
             mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>%
             filter(meth_direction == "HYPO") %>% 
             pull(ID),
     Heavy_HYPO = results_ASAT_pre_post_diff_within_twin_status1_flt %>%
             mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>%
             filter(meth_direction == "HYPO") %>% 
             pull(ID),
     DIFF_HYPO = results_ASAT_pre_post_diff_between_twin_statuses_flt %>%
             mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>%
             filter(meth_direction == "HYPO") %>% 
             pull(ID))


ggvenn(ASAT_HYPO_DMPs, set_name_size = 10, stroke_size = 1, 
               fill_color = c("#F5C242", "#353F4F", "Red"),
               text_size = 10,
               stroke_alpha = 0.8,
               set_name_color = c("Black", "Black", "Black"))+
        labs(title = "ASAT exercise response")



QF_HYPO_DMPs <- list(Lean_HYPO = results_QF_pre_post_diff_within_twin_status2_flt %>%
                               mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>%
                               filter(meth_direction == "HYPO") %>% 
                               pull(ID),
                       Heavy_HYPO = results_QF_pre_post_diff_within_twin_status1_flt %>%
                               mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>%
                               filter(meth_direction == "HYPO") %>% 
                               pull(ID),
                       DIFF_HYPO = results_QF_pre_post_diff_between_twin_statuses_flt %>%
                               mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>%
                               filter(meth_direction == "HYPO") %>% 
                               pull(ID))


ggvenn(QF_HYPO_DMPs, set_name_size = 10, stroke_size = 1, 
       fill_color = c("#F5C242", "#353F4F", "Red"),
       text_size = 10,
       stroke_alpha = 0.8,
       set_name_color = c("Black", "Black", "Black"))+
        labs(title = "QF exercise response")


####
#       HYPER
####


ASAT_HYPER_DMPs <- list(Lean_HYPER = results_ASAT_pre_post_diff_within_twin_status2_flt %>%
                               mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>%
                               filter(meth_direction == "HYPER") %>% 
                               pull(ID),
                       Heavy_HYPER = results_ASAT_pre_post_diff_within_twin_status1_flt %>%
                               mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>%
                               filter(meth_direction == "HYPER") %>% 
                               pull(ID),
                       DIFF_HYPER = results_ASAT_pre_post_diff_between_twin_statuses_flt %>%
                               mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>%
                               filter(meth_direction == "HYPER") %>% 
                               pull(ID))


ggvenn(ASAT_HYPER_DMPs, set_name_size = 10, stroke_size = 1, 
       fill_color = c("#F5C242", "#353F4F", "Red"),
       text_size = 10,
       stroke_alpha = 0.8,
       set_name_color = c("Black", "Black", "Black"))+
        labs(title = "ASAT exercise response")



QF_HYPER_DMPs <- list(Lean_HYPER = results_QF_pre_post_diff_within_twin_status2_flt %>%
                             mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>%
                             filter(meth_direction == "HYPER") %>% 
                             pull(ID),
                     Heavy_HYPER = results_QF_pre_post_diff_within_twin_status1_flt %>%
                             mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>%
                             filter(meth_direction == "HYPER") %>% 
                             pull(ID),
                     DIFF_HYPER = results_QF_pre_post_diff_between_twin_statuses_flt %>%
                             mutate(meth_direction = factor(ifelse(logFC < 0, "HYPO", "HYPER"), levels = c("HYPO", "HYPER"))) %>%
                             filter(meth_direction == "HYPER") %>% 
                             pull(ID))


ggvenn(QF_HYPER_DMPs, set_name_size = 10, stroke_size = 1, 
       fill_color = c("#F5C242", "#353F4F", "Red"),
       text_size = 10,
       stroke_alpha = 0.8,
       set_name_color = c("Black", "Black", "Black"))+
        labs(title = "QF exercise response")


