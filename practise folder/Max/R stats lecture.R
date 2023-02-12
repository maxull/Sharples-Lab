##############################################################################################
###                                                                                        ###
###     Statistics                                                                         ###
###                                                                                        ###
##############################################################################################



# topic: R stats


library(tidyverse); library(ggplot2); library(readxl)




# https://crumplab.com/rstatsmethods/articles/Stats1/Lab10_ttest.html

# t.test

?t.test

# lets extract the baseline mean and post mean from MACS RT data and treat them as a one sample t.test

rt_data <- read_excel("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/RT/RT_data.xlsx", na = "NA")


# filter values and timepoints

rt_data %>% 
        na.omit() %>% 
        filter(Session %in% c(1,16)) %>% 
        group_by(Session) %>% 
        summarise(m = mean(v_total)) -> Rt_means


t.test(Rt_means$m)


# paired

rt_data %>% 
        filter(Session == 1) -> rt1

rt_data %>% 
        filter(Session == 16)-> rt2


t.test(rt1$v_total, rt2$v_total, paired = TRUE)

# two sample

t.test(rt1$v_total, rt2$v_total, var.equal = TRUE)


# write as formula


rt_data %>% 
        na.omit() %>% 
        filter(Session %in% c(1,16)) -> rt_df

t.test(rt_df$v_total~rt_df$Session, paired = TRUE)





# https://crumplab.com/rstatsmethods/articles/Stats1/Lab11_Correlation.html

# correlation
rt_df %>% 
        select(FP, Session, v_total) %>% 
        mutate(Session = paste("Session", Session, sep = "_")) %>% 
        pivot_wider(names_from = Session, values_from = v_total) %>% 
        na.omit() -> rt_df2
ggplot(data = rt_df2, aes(x = Session_1, Session_16))+
        geom_point()+
        geom_smooth(method = "lm", se = FALSE)


cor(rt_df2$Session_1, rt_df2$Session_16)
