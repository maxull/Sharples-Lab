# emils stats

library(readxl); library(ggplot2); library(multcomp)


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


abs_df %>% 
        ggplot(aes(x = Vekt, y = `Total supine Force (N)`))+
        geom_point()
