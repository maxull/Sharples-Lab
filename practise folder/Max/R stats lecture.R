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



# https://dhammarstrom.github.io/IDR4000/lesson_12_regressionModels.html

# simple regression

# lets look at increase in total volume across sessions

# filter out 0 values and first two sessions

rt_data <- rt_data %>% 
        filter(v_total != 0, 
               Session > 2)

m1 <- lm(v_total~Session, data = rt_data)

summary(m1)
coef(m1)
confint(m1)

ggplot(data = rt_data, aes(x= Session, y = v_total))+
        geom_smooth(method = "lm")+
        geom_point()

# assumption checks by plotting the residuals agains the fitted values

rt_data$resid <- residuals(m1)
rt_data$fitted <- fitted(m1)

rt_data %>%
        ggplot(aes(fitted, resid)) + 
        geom_hline(yintercept = 0) +
        geom_point(size = 3, fill = "lightblue", shape = 21) +
        theme_minimal()


# ideal normal distribution

set.seed(1)
ggplot(data.frame(y = rnorm(100, 0, 1)), aes(sample = y)) + stat_qq(size = 3, fill = "lightblue", shape = 21) + stat_qq_line() + theme_minimal()

rt_data %>%
        mutate(st.resid = resid/sd(resid)) %>%
        ggplot(aes(sample = st.resid)) +
        stat_qq(size = 3, fill = "lightblue", shape = 21) + 
        stat_qq_line() +
        theme_minimal()

library(knitr); library(broom)

tidy(m1) %>%
        kable(col.names = c("", "Estimate", "SE", "t-statistic", "p-value"), 
              digits = c(NA, 1, 1, 2, 5))



# Multiple regression


# manipulate dataset wider

rt_change <- rt_data %>% 
        select(FP, Session, v_total) %>% 
        filter(Session %in% c(3,16)) %>% 
        mutate(Session = paste("Session", Session, sep = "_")) %>% 
        pivot_wider(names_from = Session, values_from = v_total) %>% 
        na.omit() %>% 
        mutate(change = Session_16-Session_3,
               percent_change = change/Session_3*100)

# load dexa data

dexa_data <- read_excel("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/DEXA/DEXA_data.xlsx", na = "NA")

dexa <- dexa_data %>% 
        filter(Timepoint == "Baseline") %>% 
        na.omit()

merged_df <- merge(rt_change, dexa, by = "FP")

df <- merge(rt_data, dexa, by = "FP")

m2 <- lm(v_total ~Session + Weight*Height, data = df)

summary(m2)

resid(m2)
fitted(m2)

plot(fitted(m2),resid(m2))



library(lme4)

lmer(v_total ~Session + Weight + 1|FP, data = df)
