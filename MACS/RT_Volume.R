#----------------------------------------------------------------
#
#   Analysis and plotting RT volume data from MACS project
#
#
#----------------------------------------------------------------


library(ggplot2); library(tidyverse); library(readxl);library(cowplot); library(doBy)

RT_data <- read_excel("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/RT/RT_data.xlsx", na = "NA")

RT_data %>% 
        ggplot(aes(x = Session, y = Leg_ekstensjon))+
        geom_point(size=2, aes(color = FP))+
        geom_smooth(method = "lm")


RT_mean <- RT_data %>% 
        na.omit() %>% 
        filter(v_total != 0) %>% 
        group_by(Session) %>% 
        summarise(mean = mean(v_total),
                  sd = sd(v_total))

RT_data %>% 
        filter(v_total != 0) %>% 
        na.omit() -> RT_df
        


ggplot(data = RT_df, aes(x = Session-2, y  = v_total))+
        geom_point(alpha = 0.2, size = 2)+
        geom_line(data = RT_df, aes(y = v_total, group = FP, color = FP), alpha = 0.2, size = 1.2)+
        geom_point(data = RT_mean, aes(x = Session-2, y = mean), size = 3)+
        geom_errorbar(data = RT_mean, aes(x = Session-2, ymin = (mean-sd), ymax = (mean+sd)), inherit.aes = FALSE, width = 0.2, size = 1.1)+
        geom_line(data = RT_mean, aes(y = mean), size = 1.2)+
        scale_y_continuous(expand = c(0,0),
                           limits = c(0,(max(RT_df$v_total))),
                           n.breaks = 10)+
        scale_x_continuous(n.breaks = 16)+
        theme_classic()+
        labs(y = "Volume load (KG*REPS*SETS)",
             x = "Training session")+
        geom_text(data = RT_mean, aes(label = as.integer(mean), x = Session-2, y = mean-sd), inherit.aes = FALSE, vjust = 3)



###############################################################################################

### individual plots

for (i in (unique(RT_data$FP))) {


RT_data %>% 
        filter(v_total != 0) %>% 
        na.omit() %>% 
        filter(FP == (i)) -> x

        p1 <- ggplot(data = x, aes(x = Session, y = v_total))+
        geom_point(size = 2)+
        geom_line(size = 2)+
        scale_y_continuous(expand = c(0,0),
                           limits = c(0,(max(RT_df$v_total)+100)),
                           n.breaks = 10)+
        scale_x_continuous(n.breaks = 16)+
        theme_classic()+
        labs(y = "Volume load (KG*REPS*SETS)",
             x = "Training session",
             title = "Totalt treningsvolum")
        
        p2 <- ggplot(data = x, aes(x = Session, y = Knebøy))+
                geom_point(size = 2)+
                geom_line(size = 2)+
                scale_y_continuous(expand = c(0,0),
                                   limits = c(0,(max(RT_df$Knebøy)+100)),
                                   n.breaks = 10)+
                scale_x_continuous(n.breaks = 16)+
                theme_classic()+
                labs(y = "",
                     x = "Training session",
                     title = "Knebøy volum")+
                geom_text(data = x, aes(label = as.integer(Knebøy/40)), vjust = 3)
        
        p3 <- ggplot(data = x, aes(x = Session, y = Beinpress))+
                geom_point(size = 2)+
                geom_line(size = 2)+
                scale_y_continuous(expand = c(0,0),
                                   limits = c(0,(max(RT_df$Beinpress)+100)),
                                   n.breaks = 10)+
                scale_x_continuous(n.breaks = 16)+
                theme_classic()+
                labs(y = "",
                     x = "Training session",
                     title = "Beinpress volum")+
                geom_text(data = x, aes(label = as.integer(Beinpress/40)), vjust = 3)
        
        p4 <- ggplot(data = x, aes(x = Session, y = Leg_ekstensjon))+
                geom_point(size = 2)+
                geom_line(size = 2)+
                scale_y_continuous(expand = c(0,0),
                                   limits = c(0,(max(RT_df$Leg_ekstensjon)+100)),
                                   n.breaks = 10)+
                scale_x_continuous(n.breaks = 16)+
                theme_classic()+
                labs(y = "",
                     x = "Training session",
                     title = "Leg ekstensjon volum")+
                geom_text(data = x, aes(label = as.integer(Leg_ekstensjon/40)), vjust = 3)
        
        
        p5 <- ggplot(data = x, aes(x = Session, y = Leg_curl))+
                geom_point(size = 2)+
                geom_line(size = 2)+
                scale_y_continuous(expand = c(0,0),
                                   limits = c(0,(max(RT_df$Leg_curl)+100)),
                                   n.breaks = 10)+
                scale_x_continuous(n.breaks = 16)+
                theme_classic()+
                labs(y = "",
                     x = "Training session",
                     title = "Leg curl volum")+
                geom_text(data = x, aes(label = as.integer(Leg_curl/40)), vjust = 3)
        
        p6 <- ggplot(data = x, aes(x = Session, y = Utfall))+
                geom_point(size = 2)+
                geom_line(size = 2)+
                scale_y_continuous(expand = c(0,0),
                                   limits = c(0,(max(RT_df$Utfall)+100)),
                                   n.breaks = 10)+
                scale_x_continuous(n.breaks = 16)+
                theme_classic()+
                labs(y = "",
                     x = "Training session",
                     title = "Utfall volum")+
                geom_text(data = x, aes(label = as.integer(Utfall/40)), vjust = 3)
        
        p7 <- ggplot(data = x, aes(x = Session, y = Tåhev))+
                geom_point(size = 2)+
                geom_line(size = 2)+
                scale_y_continuous(expand = c(0,100),
                                   limits = c(0,(max(RT_df$Tåhev)+100)),
                                   n.breaks = 10)+
                scale_x_continuous(n.breaks = 16)+
                theme_classic()+
                labs(y = "",
                     x = "Training session",
                     title = "Tåhev volum")+
                geom_text(data = x, aes(label = as.integer(Tåhev/40)), vjust = 3)
        
        
        p8 <- plot_grid(p2,p3,p4,p5,p6,p7, ncol = 2)
        
        p9 <- plot_grid(p1,p8, ncol = 1)
        
        ggsave(filename  = paste("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/RT/",(i),".png", sep = ""), 
               plot = p9,
               units = "cm",
               height = 27,
               width = 25)
        print(i)
}


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

# mean and sd is mean across all timepoints



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






















