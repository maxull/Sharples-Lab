#----------------------------------------------------------------
#
#   Analysis and plotting Dexa data from MACS project
#
#
#----------------------------------------------------------------


library(ggplot2); library(tidyverse); library(readxl);library(cowplot); library(doBy); library(ggsignif); library(scales); library(gtsummary)

dexa_data <- read_excel("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/DEXA/DEXA_data.xlsx")

###### plot left and right leg on same plot

dexa_data %>% 
        na.omit() %>% 
        pivot_longer(cols = 7:16, names_to = "Measurement")-> t_dexa



df <- dexa_data %>% 
        #filter(FP == "MACS_001") %>% 
        dplyr::select(1,2,6:15) %>% 
        pivot_wider(names_from = Timepoint, values_from = 3:12) %>% 
        mutate(lean_total_change = (.[[3]]-.[[2]]),
               lean_arms_change = (.[[5]]-.[[4]]),
               lean_legs_change = (.[[7]]-.[[6]]),
               lean_left_leg_change = (.[[9]]-.[[8]]),             #this variable is change
               lean_right_leg_change = (.[[11]]-.[[10]]),
               "lean_right_leg_change_%" =  (((.[[11]]-.[[10]])/.[[10]])*100),
               "lean_total_change_%" =  (((.[[3]]-.[[2]])/.[[2]])*100),
               "lean_arms_change_%" =  (((.[[5]]-.[[4]])/.[[4]])*100),
               "lean_legs_change_%" =  (((.[[7]]-.[[6]])/.[[6]])*100),
               "lean_left_leg_change_%" =  (((.[[9]]-.[[8]])/.[[8]])*100)) %>% 
        na.omit()

# t test for change from baseline

t.test(df$`Mager(g)_Total_Baseline`,df$`Mager(g)_Total_Post`, paired = TRUE)

t.test(df$`Mager(g)_Left_Leg_Post`, df$`Mager(g)_Left_Leg_Baseline`, paired = TRUE)

t.test(df$`Mager(g)_Legs_Baseline`, df$`Mager(g)_Legs_Post`, paired = TRUE)

df %>% 
        summarise(m = mean(df$'lean_legs_change_%'),
                  s = sd(df$'lean_legs_change_%'))
              
df2 <- df %>% 
        select(1, 27:31) %>% 
        pivot_longer(names_to = "percent", cols = 2:6) %>% 
        mutate("%" = gsub('.{9}$', "", percent)) %>% 
        select(1, measure = 4, percent_change = 3)

df3 <- df %>% 
        select(1, 22:26) %>% 
        pivot_longer(names_to = "percent", cols = 2:6) %>% 
        mutate("%" = gsub('.{7}$', "", percent)) %>% 
        select(1, measure = 4, change_gram = 3)

lean_dexa_change <- merge(df3, df2, by = c("FP", "measure"))




#########################################################################3
#option to use tbl_summary to create table


lean_dexa_change %>%  
        filter(FP == "MACS_001") %>% 
        pivot_longer(names_to = "change", values_to = "value", cols = 3:4) %>% 
        tbl_summary(by = change,
                    statistic = value~"{mean}")




tbl_1 <- data1 %>% 
        tbl_summary(by = "Gender",
                    statistic = c(Stature, Body_Mass, Est_VO2_max) ~"{mean} ({sd})", 
                    digits = ~2) %>% 
        add_p(vars = TRUE) %>% 
        add_ci()



#################################################################################

# plot dexa data change %
df_change <- lean_dexa_change %>% 
        na.omit()


#absolute_change <- 
        lean_dexa_change %>%
        na.omit() %>% 
        group_by(measure) %>% 
        summarise(mean_change_gram = mean(change_gram),
                  sd_change_gram = sd(change_gram),
                  mean_percent = mean(percent_change),
                  sd_percent = sd(percent_change)) %>% 
        ggplot(aes(x = measure, y = mean_change_gram))+
        geom_point(size = 3)+
        geom_errorbar(aes(ymin = (mean_change_gram - sd_change_gram),
                          ymax = (mean_change_gram + sd_change_gram)), width = 0.2)+
        geom_point(data = lean_dexa_change, aes(x = measure, y = change_gram, color = FP), size = 3, alpha = 0.5)+
        theme_classic(base_size = 20)


#percent_change <- 
df_percent <- lean_dexa_change %>%
        na.omit() %>% 
        group_by(measure) %>% 
        summarise(mean_change_gram = mean(change_gram)/100,
                  sd_change_gram = sd(change_gram)/100,
                  mean_percent = mean(percent_change)/100,
                  sd_percent = sd(percent_change)/100) %>% 
        mutate(measure = factor(measure, levels = c("lean_arms", "lean_legs", "lean_left_leg", "lean_right_leg", "lean_total")))

        
dexa_plot <- 
ggplot(df_percent, aes(x = measure, y = mean_percent))+
        geom_point(size = 4)+
        geom_errorbar(data = df_percent, aes(ymin = (mean_percent - sd_percent),
                          ymax = (mean_percent + sd_percent)), width = 0.2, size = 1.2)+
        geom_text(data = df_percent, aes(label = paste0(format(round(mean_percent*100, 2), nsmall = 2),"%"), fontface = "bold", hjust = -0.3), size = 6)+
        geom_point(data = df_change, aes(x = measure, y = percent_change/100, color = FP), size = 3)+
        theme_classic(base_size = 20)+
        scale_y_continuous(labels = percent,
                           limits = c(-0.005,0.13),
                           expand = c(0,0))+
        scale_x_discrete(labels = c("Arms","Legs","Left Leg","Right Leg", "Whole Body"))+
        theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 0.9),
              axis.title.x = element_blank(),
              plot.title = element_text(size = 15),
              legend.position = "none")+
        labs(y = "% Change")


plot_grid(isom_plot, dexa_plot, ncol = 1)
