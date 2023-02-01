##############################################################################################
###                                                                                        ###
###     Data manipulation -> plots                                                         ###
###                                                                                        ###
##############################################################################################

# https://dhammarstrom.github.io/IDR4000/lesson_6_make_summaries.html

# topic: Data manipulation to gain insight and make plots


library(tidyverse); library(ggplot2); library(readxl)


#  load data

rt_data <- read_excel("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/RT/MACS_001.xlsx")

# open data understand the variables

view(rt_data)

summary(rt_data)

rt_data


# lets gain som insight



rt_data %>% 
        na.omit() %>% 
        group_by(Session) %>% 
        summarise(mean. = mean(v_total),
                  sd = sd(v_total)) %>% 
        ggplot(aes(x = Session, y = mean.))+
        geom_point()+
        scale_y_continuous(limits = c(0,18000), n.breaks = 10)+
        scale_x_continuous(n.breaks = 16)+
        geom_smooth(method = "loess")




# calculate change and % change from first real session. aka. session 3, session 1 and 2 were familiarization sessions

rt_data %>% 
        na.omit() %>% 
        filter(Session >= 3) %>% 
        mutate(Session = as.factor(Session)) %>% 
        select(FP, Session, v_total) %>% 
        pivot_wider(names_from = Session, values_from = v_total) %>% 
        mutate(change = '15'-'3')


rt_data %>% 
        na.omit() %>% 
        filter(Session >= 3) %>% 
        mutate(Session = paste("Session", Session, sep = "_")) %>% 
        select(FP, Session, v_total) %>% 
        pivot_wider(names_from = Session, values_from = v_total) %>% 
        mutate(Change = Session_15 - Session_3,
               percent_change = (Change/Session_3)*100) %>% 
        select(FP, Session_15, percent_change) %>% 
        mutate(Session = x) -> rt_data2

x = 15


rt_data %>% 
        na.omit() %>% 
        group_by(Session) %>% 
        summarise(mean. = mean(v_total),
                  sd = sd(v_total)) -> rt_mean

rt_data %>% 
        filter(v_total > 0) %>% 
        ggplot(aes(x= Session, y = v_total))+
        geom_point(aes(color = FP), size = 2, alpha = 0.3)+
        geom_line(aes(color = FP), size = 2, alpha = 0.3)+
        geom_point(data=rt_mean, aes(y=mean.), size = 3)+
        geom_errorbar(data = rt_mean, aes(ymin=mean.-sd, ymax=mean.+sd, x= Session), inherit.aes = FALSE, width =0.2, size = 1.2)+
        geom_line(data = rt_mean, aes(y=mean.), size = 2)+
        scale_y_continuous(limits = c(0,28000), n.breaks = 14)
