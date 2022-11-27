# hjelpe ba studentene

library(ggplot2); library(tidyverse); library(readxl);library(cowplot); library(doBy); library(ggsignif); library(scales)



lab_df <- read_excel("C:/Users/maxul/Documents/lab.xlsx")




###




lab_df %>% 
        mutate(Timepoint = factor(Timepoint, levels = c("pre", "post")),
               Condition = factor(Condition)) %>% 
        group_by(Gene, Timepoint,Condition) %>% 
        summarise(mean = mean(foldchange),
                  sd = sd(foldchange)) %>% 
        ggplot(aes(x = Timepoint, y = mean, group = Condition, color = Condition))+
        geom_point(position = position_dodge(0.4), size = 3)+
        geom_errorbar(aes(ymin = (mean - sd),
                          ymax = (mean + sd)), position = position_dodge(0.4), width = 0.2, size = 1)+
        facet_grid(~Gene)+
        theme_bw()+
        theme(axis.title.x = element_blank())+
        labs(y = "foldchange")+
        geom_signif(comparisons = list(c("pre", "post")),
                    annotations = "****",
                    color = "black",
                    tip_length = 0.01)



ar <- lab_df %>% 
        mutate(Timepoint = factor(Timepoint, levels = c("pre", "post")),
               Condition = factor(Condition)) %>% 
        group_by(Gene, Timepoint,Condition) %>% 
        summarise(mean = mean(foldchange),
                  sd = sd(foldchange)) %>%
        filter(Gene == "AR",
               Timepoint == "post") %>% 
        mutate(mean = mean - 1) %>% 
        ggplot(aes(x = Condition, y = mean))+
        geom_point(size = 2)+
        geom_text(aes(label = c("+121%", "+117%")), hjust = -0.5)+
        geom_errorbar(aes(ymin = (mean - sd),
                          ymax = (mean + sd)), width = 0.2)+
        scale_y_continuous(labels = percent,
                           limits = c(0,1.5))+
        geom_hline(yintercept = 0, linetype = 2, size = 1.5)+
        theme_classic(base_size = 15)+
        theme(axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.text.x = element_text(size = 15),
              plot.title = element_text(size = 15))+
        labs(y = "AR % change",
             title = "AR gene expression after acute resistance exercise")+
        geom_signif(comparisons = list(c("40% - 20 Rep", "80% - 10 Rep")),
                    annotations = "p = 0.95",
                    tip_length = 0.01)
        


igf <- 
lab_df %>% 
        mutate(Timepoint = factor(Timepoint, levels = c("pre", "post")),
               Condition = factor(Condition)) %>% 
        group_by(Gene, Timepoint,Condition) %>% 
        summarise(mean = mean(foldchange),
                  sd = sd(foldchange)) %>%
        filter(Gene == "IGF-1",
               Timepoint == "post") %>% 
        mutate(mean = mean - 1) %>% 
        ggplot(aes(x = Condition, y = mean))+
        geom_point(size = 2)+
        geom_text(aes(label = c("+154%", "+226%")), hjust = -0.5)+
        geom_errorbar(aes(ymin = (mean - sd),
                          ymax = (mean + sd)), width = 0.2)+
        scale_y_continuous(labels = percent,
                           limits = c(0,2.5))+
        geom_hline(yintercept = 0, linetype = 2, size = 1.5)+
        theme_classic(base_size = 15)+
        theme(axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.text.x = element_text(size = 15),
              plot.title = element_text(size = 15))+
        labs(y = "IGF-1 % change",
             title = "IGF-1 gene expression after acute resistance exercise")+
        geom_signif(comparisons = list(c("40% - 20 Rep", "80% - 10 Rep")),
                    annotations = "p = 0.0006",
                    tip_length = 0.01)



plot_grid(igf, ar, ncol = 2)







