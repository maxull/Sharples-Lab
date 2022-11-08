#----------------------------------------------------------------
#
#   Analysis and plotting CSA data from MACS project
#
#
#----------------------------------------------------------------


library(ggplot2); library(tidyverse); library(readxl);library(cowplot)

csa_data <- read_excel("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/CSA/MACS_CSA.xlsx")


### plot RF

csa_data %>%
        filter(muscle == "RF") %>% 
        ggplot(aes(x = timepoint, y = cm2))+ 
        geom_boxplot()+
        expand_limits(y = 0)

csa_data %>%
        filter(muscle == "VL") %>% 
        ggplot(aes(x = timepoint, y = cm2))+ 
        geom_boxplot()+
        geom_point()+
        expand_limits(y = 0)

#### create mean and sd line and point chart for participant MACS_001
        
p1 <- csa_data %>%
        filter(FP == "MACS_001") %>% 
        group_by(timepoint, muscle) %>% 
        summarise(mean = mean(cm2),
                  sd = sd(cm2)) %>% 
        filter(muscle == "RF") %>% 
        ggplot(aes(x = timepoint, y = mean, group = muscle))+
        geom_errorbar(aes(ymin = (mean-sd),
                          ymax = (mean+sd)),
                      width = 0.2)+
        geom_point()+
        geom_line()+
        theme_classic()+
        labs(y = "cm^2 (mean ± sd)",
             x = "",
             title = "MACS_001 RF")+
        theme(axis.title.y = element_text(face = "bold"))
        
p2 <- csa_data %>%
        filter(FP == "MACS_001") %>% 
        group_by(timepoint, muscle) %>% 
        summarise(mean = mean(cm2),
                  sd = sd(cm2)) %>% 
        filter(muscle == "VL") %>% 
        ggplot(aes(x = timepoint, y = mean, group = muscle))+
        geom_errorbar(aes(ymin = (mean-sd),
                          ymax = (mean+sd)),
                      width = 0.2)+
        geom_point()+
        geom_line()+
        theme_classic()+
        labs(y = "",
             x = "",
             title = "MACS_001 VL")

plot_grid(p1, p2, ncol = 2)





