#----------------------------------------------------------------
#
#   Analysis and plotting Dexa data from MACS project
#
#
#----------------------------------------------------------------


library(ggplot2); library(tidyverse); library(readxl);library(cowplot)

dexa_data <- read_excel("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/DEXA_data.xlsx")


dexa_data %>% 
        na.omit() %>% 
        ggplot(aes(Height, Weight, color = FP))+
        geom_point(size = 2)+
        theme_bw()

dexa_data %>% 
        na.omit() %>% 
        ggplot(aes(x = Timepoint, y = `Mager(g)_Total`, color = FP))+
        geom_point()+
        theme_bw()


dexa_data %>% 
        na.omit() %>% 
        ggplot(aes(x = Timepoint, y = `Mager(g)_Legs`, color = FP))+
        geom_point()+
        theme_bw()


###### plot left and right leg on same plot

dexa_data %>% 
        na.omit() %>% 
        pivot_longer(cols = 7:16, names_to = "Measurement")-> t_dexa


t_dexa %>% 
        ggplot(aes(Measurement, value, group = FP, color = FP))+
        geom_point(position = position_dodge(width = 0.2))+
        geom_line(position = position_dodge(width = 0.2))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

### filter just legs, lean mass

t_dexa %>% 
        filter(Measurement %in% c("Mager(g)_Left_Leg",
                                "Mager(g)_Right_Leg")) %>% 
        ggplot(aes(Measurement, value, group = FP, color = FP))+
        geom_point(position = position_dodge(width = 0.2))+
        geom_line(position = position_dodge(width = 0.2))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
        scale_x_discrete(labels = c("Left Leg", "Right Leg"))+
        ylab("Lean mass (g)")

