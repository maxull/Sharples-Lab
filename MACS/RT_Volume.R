#----------------------------------------------------------------
#
#   Analysis and plotting RT volume data from MACS project
#
#
#----------------------------------------------------------------


library(ggplot2); library(tidyverse); library(readxl);library(cowplot)

RT_data <- read_excel("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/RT/MACS_001.xlsx")

RT_data %>% 
        ggplot(aes(x = Week, y = Knebøy))+
        geom_point(size=2)+
        geom_line()+
        geom_smooth(method = "loess")

### without familiarization week

p1 <- RT_data %>% 
        tail(14) %>% 
        ggplot(aes(x = Week, y = Knebøy))+
        geom_point(size=1.5)+
        geom_line()+
        geom_smooth(method = "loess")+
        labs(y = "KG (weight*reps*sets)",
             title = "Squat volume")+
        theme_classic()+
        theme(axis.title.x = element_blank())
        

p2 <- RT_data %>% 
        tail(14) %>% 
        ggplot(aes(x = Week, y = Beinpress))+
        geom_point(size=1.5)+
        geom_line()+
        geom_smooth(method = "loess")+
        labs(title = "Legpress volume")+
        theme_classic()+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())

p3 <- RT_data %>% 
        tail(14) %>% 
        ggplot(aes(x = Week, y = `Leg ekstensjon`))+
        geom_point(size=1.5)+
        geom_line()+
        geom_smooth(method = "loess")+
        labs(y = "KG (weight*reps*sets)",
             title = "Leg extension volume",
             x = "Session")+
        theme_classic()

p4 <- RT_data %>% 
        tail(14) %>% 
        ggplot(aes(x = Week, y = `Leg curl`))+
        geom_point(size=1.5)+
        geom_line()+
        geom_smooth(method = "loess")+
        labs(title = "Leg curl volume",
             x = "Session")+
        theme(axis.title.y = element_blank())+
        theme_classic()


p5 <- RT_data %>% 
        tail(14) %>% 
        ggplot(aes(x = Week, y = v_total))+
        geom_point(size=2)+
        geom_line()+
        geom_smooth(method = "loess")+
        labs(y = "KG (weight*reps*sets)",
             title = "Total volume legs",
             x = "Session")+
        theme_classic()

pgrid <- plot_grid(p1,p2,p3,p4, nrow = 2, ncol = 2)

figure1 <- plot_grid(pgrid, p5, ncol = 1)

ggsave("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/RT/MACS_001(2).png", plot = figure1)
