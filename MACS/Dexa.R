#----------------------------------------------------------------
#
#   Analysis and plotting Dexa data from MACS project
#
#
#----------------------------------------------------------------


library(ggplot2); library(tidyverse); library(readxl);library(cowplot); library(gtsummary)

dexa_data <- read_excel("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/DEXA/DEXA_data.xlsx")


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

dexa_data %>% 
        pivot_wider(names_from = Timepoint, values_from = 7:16)
    


df <- dexa_data %>% 
        #filter(FP == "MACS_001") %>% 
        select(1,3,7:16) %>% 
        pivot_wider(names_from = Timepoint, values_from = 3:12) %>% 
        mutate(lean_total_change = (.[[3]]-.[[2]]),
               lean_arms_change = (.[[5]]-.[[4]]),
               lean_legs_change = (.[[7]]-.[[6]]),
               lean_left_leg_change = (.[[9]]-.[[8]]),
               lean_right_leg_change = (.[[11]]-.[[10]]),
               "lean_right_leg_change_%" =  (((.[[11]]-.[[10]])/.[[10]])*100),
               "lean_total_change_%" =  (((.[[3]]-.[[2]])/.[[2]])*100),
               "lean_arms_change_%" =  (((.[[5]]-.[[4]])/.[[4]])*100),
               "lean_legs_change_%" =  (((.[[7]]-.[[6]])/.[[6]])*100),
               "lean_left_leg_change_%" =  (((.[[9]]-.[[8]])/.[[8]])*100))
              
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


df4 <- t_dexa %>% 
        select(1,3,7,8) %>% 
        pivot_wider(names_from = Timepoint, values_from = value) %>% 
        filter(FP == "MACS_001") %>% 
        head(5) %>% 
        slice(2:5,1) %>% 
        select(1,2,3,4)

df5 <- lean_dexa_change %>% 
        head(5)


#########################################################################3
#option to use tbl_summary to create table


t_dexa %>% 
        select(1,3,7,gram = 8) %>% 
        filter(FP == "MACS_001",
               Measurement == m_df[1:5]) %>% 
        tbl_summary(by = "Timepoint")

dexa_data %>% 
        select(1,3,7:11) %>% 
        filter(FP == "MACS_001") %>% 
        tbl_summary(by = "Timepoint")




tbl_1 <- data1 %>% 
        tbl_summary(by = "Gender",
                    statistic = c(Stature, Body_Mass, Est_VO2_max) ~"{mean} ({sd})", 
                    digits = ~2) %>% 
        add_p(vars = TRUE) %>% 
        add_ci()



#################################################################################




t_dexa %>% 
        ggplot(aes(Measurement, value, group = FP, color = FP))+
        geom_point(position = position_dodge(width = 0.2))+
        geom_line(position = position_dodge(width = 0.2))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

### filter just legs, lean mass

t_dexa %>% 
        filter(Measurement %in% c("Mager(g)_Left_Leg",
                                "Mager(g)_Right_Leg")) %>% 
        ggplot(aes(Measurement, value, group = Timepoint, color = FP))+
        geom_point(position = position_dodge(width = 0.2))+
        geom_line(position = position_dodge(width = 0.2))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
        scale_x_discrete(labels = c("Left Leg", "Right Leg"))+
        ylab("Lean mass (g)")

