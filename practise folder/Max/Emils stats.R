# emils stats

library(readxl); library(ggplot2); library(multcomp); library(tidyverse)


kroppsvekt <- read_excel("/Users/maxul/Downloads/Pre_kroppsvekt.xlsx")

abs_df <- read_excel("/Users/maxul/Downloads/Pre_absolutt.xlsx")

plot_supine <- abs_df %>% 
        ggplot(aes(x = Age, y = `Total supine Force (N)`))+
        geom_point()
plot_supine

model <- aov(`Total supine Force (N)`~Age, data = abs_df)

summary(model)

# Tukey HSD test:
post_test <- glht(model,
                  linfct = mcp(as.factor(Age) = "Tukey"))

summary(post_test)


mod <- lm(`Total supine Force (N)`~Age + Vekt, data = abs_df)

summary(mod)


mod2 <- lm(abs_df$`Total 90deg Force (N)` ~ Age + Vekt, data = abs_df)

summary(mod2)


mod3 <- lm(abs_df$`Total Force CF (N)` ~ Age + Vekt, data = abs_df)

summary(mod3)

mod3 <- lm(abs_df$`Total Force CF (N)` ~ Position, data = abs_df)


abs_df %>% 
        ggplot(aes(x = Position, y = `Total Force CF (N)`))+
        geom_point()

r2 <- (cor(x = abs_df$Vekt, abs_df$`Total Force CF (N)`)^2)

