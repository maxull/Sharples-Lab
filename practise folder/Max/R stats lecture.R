##############################################################################################
###                                                                                        ###
###     Statistics                                                                         ###
###                                                                                        ###
##############################################################################################



# topic: R stats


library(tidyverse); library(ggplot2); library(readxl)


# t.test
################################################################################## 
# https://crumplab.com/rstatsmethods/articles/Stats1/Lab10_ttest.html



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

rt_data %>% 
        na.omit() %>% 
        filter(Session %in% c(1,16)) %>% 
        t.test(v_total~Session, paired = TRUE, data = .)

rt_data %>% 
        na.omit() %>% 
        filter(Session %in% c(1,16)) %>%
        ggplot(aes(x = as.factor(Session), y = v_total, group = as.factor(Session)))+
        geom_boxplot()

#####
# correlation
###################################################################################
# https://crumplab.com/rstatsmethods/articles/Stats1/Lab11_Correlation.html


rt_df %>% 
        select(FP, Session, v_total) %>% 
        mutate(Session = paste("Session", Session, sep = "_")) %>% 
        pivot_wider(names_from = Session, values_from = v_total) %>% 
        na.omit() -> rt_df2
ggplot(data = rt_df2, aes(x = Session_1, Session_16))+
        geom_point()+
        geom_smooth(method = "lm", se = FALSE)


r2 <- (cor(rt_df2$Session_1, rt_df2$Session_16))^2


#####
# simple regression
###################################################################################
# https://dhammarstrom.github.io/IDR4000/lesson_12_regressionModels.html



# lets look at increase in total volume across sessions

# filter out 0 values and first two sessions

rt_data <- rt_data %>% 
        filter(v_total != 0, 
               Session > 2)

?lm
m1 <- lm(v_total~Session, data = rt_data)

attributes(m1)


summary(m1)
coef(m1)
confint(m1)

ggplot(data = rt_data, aes(x= Session, y = v_total))+
        geom_smooth(method = "lm", se = FALSE)+
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


#####
# Multiple regression
###################################################################################


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





#####
# anova
###################################################################################

# load data with multiple groups and timepoints

contratrain <- read_excel("/Users/maxul/Documents/Skole/Master 21-22/Master/Contratrain/tr014_ultrasound.xlsx", na = "NA")

# extract only set and change from T1 to T4

contratrain <- contratrain %>% 
        mutate(change = T4-T1)

# anova change per excercise sets

mod <- aov(change~as.factor(set), data = contratrain)

summary(mod)

plot(fitted(mod), resid(mod))

hist(resid(mod))

attributes(mod)

# which group comparisons are significant

TukeyHSD(mod)


#####
# ancova
###################################################################################

# https://www.datanovia.com/en/lessons/ancova-in-r/#:~:text=The%20Analysis%20of%20Covariance%20(ANCOVA,two%20or%20more%20independent%20groups.

# An ancova
#m2 <- lm(change ~ pre + group, data = change.data)

contratrain <- contratrain %>% 
        mutate(set = as.factor(set))

mod2 <- lm(change ~T1 + set, data = contratrain)

summary(mod2)

hist(resid(mod2))

#####
# repeated measures anova
###################################################################################

# continue with dataset from anova

# repeated measures

# visualize first

contratrain %>% 
        select(1:9) %>% 
        pivot_longer(names_to = "timepoint", values_to = "value", cols = 6:9) %>% 
        na.omit() %>% 
        ggplot(aes(x = timepoint, y = value, color = set))+
        geom_point()+
        geom_smooth(method = lm,aes(group = set), se = FALSE)+
        facet_wrap(~set)

contratrain <- contratrain %>% 
        mutate(set = as.factor(set))

contratrain %>% 
        select(1:9) %>% 
        pivot_longer(names_to = "timepoint", values_to = "value", cols = 6:9) %>% 
        na.omit() %>%
        lm(value~timepoint + set, data = .) -> mod3 

summary(mod3)



#####
# linear mixed models
###################################################################################
library(lme4); library(nlme); library(emmeans)

lmer(v_total ~Session + Weight + 1|FP, data = df)
