#----------------------------------------------------------------
#
#   Analysis and plotting RT volume data from MACS project
#
#
#----------------------------------------------------------------


library(ggplot2); library(tidyverse); library(readxl);library(cowplot); library(doBy)

RT_data <- read_excel("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/RT/MACS_001.xlsx")

RT_data %>% 
        tail(15) %>% 
        ggplot(aes(x = Week, y = `Leg ekstensjon`))+
        geom_point(size=2)+
        geom_line()+
        geom_smooth(method = "lm")

### without familiarization week

p1 <- 
RT_data %>% 
        tail(14) %>% 
        ggplot(aes(x = Week, y = Knebøy))+
        geom_point(size=1.5)+
        geom_line()+
        geom_smooth(method = "lm", se = FALSE)+
        labs(y = "KG (weight*reps*sets)",
             title = "Squat volume")+
        theme_classic()+
        theme(axis.title.x = element_blank())+
        theme(axis.title = element_text(size = rel(0.8)),
              title = element_text(size = rel(0.7)))

p2 <- 
        RT_data %>% 
        tail(15) %>% 
        ggplot(aes(x = Week, y = Knebøy/Average_RPE))+
        geom_point(size=1.5)+
        geom_line()+
        geom_smooth(method = "lm", se = FALSE)+
        labs(y = "KG (weight*reps*sets)",
             title = "Squat volume/RPE")+
        theme_classic()+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        theme(axis.title = element_text(size = rel(0.8)),
              title = element_text(size = rel(0.7)))  

p3 <- RT_data %>% 
        tail(15) %>% 
        ggplot(aes(x = Week, y = Beinpress))+
        geom_point(size=1.5)+
        geom_line()+
        geom_smooth(method = "lm", se = FALSE)+
        labs(title = "Legpress volume",
             y = "KG (weight*reps*sets)")+
        theme_classic()+
        theme(axis.title.x = element_blank())+
        theme(axis.title = element_text(size = rel(0.8)),
              title = element_text(size = rel(0.7)))

p4 <- RT_data %>% 
        tail(15) %>% 
        ggplot(aes(x = Week, y = Beinpress/Average_RPE))+
        geom_point(size=1.5)+
        geom_line()+
        geom_smooth(method = "lm", se = FALSE)+
        labs(title = "Legpress volume/RPE")+
        theme_classic()+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        theme(axis.title = element_text(size = rel(0.8)),
              title = element_text(size = rel(0.7)))

p5 <- RT_data %>% 
        tail(15) %>% 
        ggplot(aes(x = Week, y = `Leg ekstensjon`))+
        geom_point(size=1.5)+
        geom_line()+
        geom_smooth(method = "lm", se = FALSE)+
        theme_classic()+
        labs(y = "KG (weight*reps*sets)",
             title = "Leg extension volume")+
        theme(axis.title.x = element_blank())+
        theme(axis.title = element_text(size = rel(0.8)),
              title = element_text(size = rel(0.7)))
        

p6 <- RT_data %>% 
        tail(15) %>% 
        ggplot(aes(x = Week, y = `Leg ekstensjon`/Average_RPE))+
        geom_point(size=1.5)+
        geom_line()+
        geom_smooth(method = "lm", se = FALSE)+
        theme_classic()+
        labs(title = "Leg extension volume/RPE")+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        theme(axis.title = element_text(size = rel(0.8)),
              title = element_text(size = rel(0.7)))
        



p7 <- RT_data %>% 
        tail(15) %>% 
        ggplot(aes(x = Week, y = v_total))+
        geom_point(size=2)+
        geom_line()+
        geom_smooth(method = "lm", se = FALSE)+
        labs(y = "KG (weight*reps*sets)",
             title = "Total volume legs",
             x = "Session")+
        theme_classic()+
        theme(axis.title = element_text(size = rel(0.8)),
              title = element_text(size = rel(0.7)))


p8 <- RT_data %>% 
        tail(15) %>% 
        ggplot(aes(x = Week, y = v_total/Average_RPE))+
        geom_point(size=2)+
        geom_line()+
        geom_smooth(method = "lm", se = FALSE)+
        theme_classic()+
        labs(title = "Total volume legs/RPE",
             x = "Session")+
        theme(axis.title.y = element_blank())+
        theme(axis.title = element_text(size = rel(0.8)),
              title = element_text(size = rel(0.7)))
       

#pgrid <- 
figure1 <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 2)


ggsave("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/RT/MACS_001(3).png", 
       plot = figure1,
       units = "px",
       height = 2000,
       width = 1800)

ggsave("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/isom_plot.png", 
       plot = isom_plot,
       units = "px",
       dpi = 600,
       bg = "white")
###############################################################3

### plot total volume for all participants with mean and sd

### get data from RT sheets and remove familiarization timepoints

# plot only leg volume

# mean and sd is mean across all timepoints... fix this

mean_RT <- summaryBy(v_total ~Week, data = RT_data, FUN = c(mean, sd)) %>% 
        filter(Week >= 3)

mean_RT <- RT_data %>% 
        group_by(Week) %>% 
        summarize(m = mean(v_total))

R_data <- RT_data %>% 
        filter(Week >= 3)
        
ggplot()+
        geom_point(aes(color = FP), alpha = 0.2)+
        geom_line(aes(color = FP), alpha = 0.2)+
        geom_point(mean_RT, aes( y = m), size = 2)
        geom_line(aes(x = Week, y = mean(v_total)))+
        geom_errorbar(aes(x = Week, ymin = mean(v_total) - sd(v_total), ymax = mean(v_total) + sd(v_total)), width = .2)+
        labs(y = "KG (weight*reps*sets)",
             title = "Total volume legs",
             x = "Session")+
        theme_classic()


RT_data %>% 
                tail(14) %>% 
                ggplot(aes(x = Week, y = v_total/Average_RPE))+
                geom_point(size=2)+
                geom_line()+
                geom_smooth(method = "lm")+
                labs(y = "KG (weight*reps*sets)",
                     title = "Total volume legs",
                     x = "Session")+
                theme_classic()



RT_data2 <- read_excel("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/RT/MACS_001.xlsx", sheet = "Sheet2")

RT_data2 %>% 
        ggplot(aes(y = Weight, x = Session))+
        geom_point()+
        geom_line()+
        geom_smooth(method = "lm", se = FALSE)

m = lm(RT_data2$Weight ~RT_data2$Session)        
m


t_RT <- RT_data %>% 
        select(!2)

RT_df <- t(t_RT) %>% 
        as.data.frame() %>% 
        select(!1) %>% 
        mutate(v_change = V16-V2,
               "%_change" = ((v_change/V2)*100))






















