library(tidyverse)
library(multcomp)
library(readxl)
library(ggplot2)
library(tidyr)
library(psych)
library(dplyr)
library(dunn.test)
library(rlang)


dataset <- read_excel("practise folder/Chris/Data/thesisdatasum.xlsx", na ="NA")

#------------------------------------------------


# SUMMARY FOR DIAMETER + HISTOGRAM + STATISTICAL TEST

diameter_data <- subset(dataset, morph_type == "diameter")

# Calculate summary statistics
summary_table <- diameter_data %>%
        group_by(condition, morph_type) %>%
        summarize(
                Frequency = n(),
                Mean = mean(value),
                SD = sd(value),
                Min = min(value),
                Max = max(value)
        )

# Print the summary table
print(summary_table)

normality_test <- shapiro.test(diameter_data$value)
cat("Normality Test for Diameter:\n")
print(normality_test)


ggplot(diameter_data, aes(x = value)) +
        geom_histogram(fill = "steelblue", color = "white", bins = 20) +
        labs(title = "Histogram of Diameter") +
        xlab("Diameter") +
        ylab("Frequency")


#KRUSKAL-Wallis test
kruskal_test <- kruskal.test(value ~ condition, data = diameter_data)
cat("Kruskal-Wallis Test for Diameter:\n")
print(kruskal_test)


# POST HOC PAIRWIZE TESTS:

        #DUNN TEST
        pairwise_dunn_test <- dunn.test(diameter_data$value, diameter_data$condition, method = "bonferroni")
        print(pairwise_dunn_test)



        #Mann-Whitney U tests
        pairwise_mannwhitney_test <- pairwise.wilcox.test(diameter_data$value, diameter_data$condition,
                                      p.adjust.method = "bonferroni")
        print(pairwise_mannwhitney_test)
        

#Pooling MA + ME and LT + LTA to one condition and running the same analysis:
        
        diameter_pooled <- diameter_data
        diameter_pooled$condition[diameter_pooled$condition %in% c("LT", "LTA")] <- "LATE"
        diameter_pooled$condition[diameter_pooled$condition %in% c("MA", "ME")] <- "MEM"
        
        diameter_pooled$condition <- factor(diameter_pooled$condition,
                                        levels = c("ET", "LATE", "MEM", "CT"),
                                        labels = c("ET", "LATE", "MEM", "CT"))
        
        # Below is another inferior way to pool the data:
        # diameter_data$group <- ifelse(diameter_data$condition == "ET", "ET",
          #                            ifelse(diameter_data$condition %in% c("LTA", "LT"), "LATE",
           #                                  ifelse(diameter_data$condition %in% c("MA", "ME"), "MEM", "CT")))
        
        #unique(diameter_data$group) #checking that the previous grouping worked
        
        #KRUSKAL-WALLIS
        kruskal_test_diameter_pooled <- kruskal.test(value ~ condition, data = diameter_pooled)
        print(kruskal_test_diameter_pooled)
        
        #DUNN TEST
        pairwise_dunn_test_dp <- dunn.test(diameter_pooled$value, diameter_pooled$condition, method = "bonferroni")
        print(pairwise_dunn_test_dp)
        
        #MANN-WHITNEY
        pairwise_mannwhitney_test_dp <- pairwise.wilcox.test(diameter_pooled$value, diameter_pooled$condition, p.adjust.method = "bonferroni")
        print(pairwise_mannwhitney_test_dp)
        
        #BOXPLOT
        ggplot(diameter_pooled, aes(x = condition, y = value, fill = condition)) +
                geom_boxplot() +
                stat_summary(fun.data = mean_sdl, geom = "errorbar", width = 0.2) +
                labs(x = "Pooled Conditions", y = "microns", title = "Boxplot of Diameter by Pooled Condition") +
                scale_fill_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        
        #BARPLOT
        ggplot(diameter_pooled, aes(x = condition, y = value, fill = condition)) +
                geom_bar(stat = "identity", color = "black") +
                labs(x = "Pooled Conditions", y = "microns", title = "Histogram of Diameter by Pooled Condition") +
                scale_fill_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        
        #DENSITYPLOT
        ggplot(area_pooled, aes(x = value, fill = condition)) +
                geom_density(alpha = 0.5) +
                labs(x = "microns", y = "Density", title = "Density Plot of Diameter by Pooled Condition") +
                scale_fill_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        
        #SWARM PLOT
        ggplot(diameter_pooled, aes(x = condition, y = value, color = condition)) +
                geom_point(position = position_jitterdodge()) +
                labs(x = "Pooled Conditions", y = "microns", title = "Swarm Plot of Diameter by Pooled Condition") +
                scale_color_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        
      
        
        
#SUMMARY FOR AREA + HISTOGRAM + STATISTICAL TEST
area_data <- subset(dataset, morph_type == "area")

area_summary <- summary(area_data$value)
area_sd <- sd(area_data$value)

# Calculate summary statistics
area_summary_table <- area_data %>%
        group_by(condition, morph_type) %>%
        summarize(
                Frequency = n(),
                Mean = mean(value),
                SD = sd(value),
                Min = min(value),
                Max = max(value)
        )

cat("Summary Statistics for Area:\n")
print(area_summary_table)

normality_test <- shapiro.test(area_data$value)
cat("Normality Test for Area:\n")
print(normality_test)


ggplot(area_data, aes(x = value)) +
        geom_histogram(fill = "steelblue", color = "white", bins = 20) +
        labs(title = "Histogram of Area") +
        xlab("Area") +
        ylab("Frequency")



#KRUSKAL-Wallis test
kruskal_test_area <- kruskal.test(value ~ condition, data = area_data)
cat("Kruskal-Wallis Test for Area:\n")
print(kruskal_test)


        # POST HOC PAIRWIZE TESTS:

        #DUNN TEST
        pairwise_dunn_test_area <- dunn.test(area_data$value, area_data$condition, method = "bonferroni")
        print(pairwise_dunn_test_area)



        #Mann-Whitney U tests
        pairwise_mannwhitney_test_area <- pairwise.wilcox.test(area_data$value, area_data$condition,
                                                  p.adjust.method = "bonferroni")
        print(pairwise_mannwhitney_test_area)

        
        
#Pooling MA + ME and LT + LTA to one condition and running the same analysis:
        
        area_pooled <- area_data
        area_pooled$condition[area_pooled$condition %in% c("LT", "LTA")] <- "LATE"
        area_pooled$condition[area_pooled$condition %in% c("MA", "ME")] <- "MEM"
        
        area_pooled$condition <- factor(area_pooled$condition,
                                        levels = c("ET", "LATE", "MEM", "CT"),
                                        labels = c("ET", "LATE", "MEM", "CT"))
        
        #The following is another way to group but the above is simpler and more tidy.
        # area_data$group <- ifelse(area_data$condition == "ET", "ET",
                                     # ifelse(area_data$condition %in% c("LTA", "LT"), "LATE",
                                          #   ifelse(area_data$condition %in% c("MA", "ME"), "MEM", "CT")))
        
        #unique(area_data$group) #checking that the previous grouping worked
        
        #BARPLOT
        ggplot(area_pooled, aes(x = condition, y = value/1000000, fill = condition)) +
                geom_bar(stat = "identity", color = "black") +
                labs(x = "Pooled Conditions", y = "mm^2", title = "Histogram of Values by Pooled Condition") +
                scale_fill_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        
        #BOXPLOT
        ggplot(area_pooled, aes(x = condition, y = value/1000000, fill = condition)) +
                geom_boxplot() +
               # stat_summary(fun.data = mean_sdl, geom = "errorbar", width = 0.2) +
                labs(x = "Pooled Conditions", y = "mm^2", title = "Boxplot of Values by Pooled Condition") +
                scale_fill_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        
        #DENSITYPLOT
        ggplot(area_pooled, aes(x = value/1000000, fill = condition)) +
                geom_density(alpha = 0.5) +
                labs(x = "mm^2", y = "Density", title = "Density Plot of Values by Pooled Condition") +
                scale_fill_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        
        #SWARM POLOT
        ggplot(area_pooled, aes(x = condition, y = value/1000000, color = condition)) +
                geom_point(position = position_jitterdodge()) +
                labs(x = "Pooled Conditions", y = "mm^2", title = "Swarm Plot of Area by Pooled Condition") +
                scale_color_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        
        # QQ plot
        qqnorm(residuals(lm(value ~ 1, data = area_pooled)))
        qqline(residuals(lm(value ~ 1, data = area_pooled)))
        
        # Residual plot
        plot(fitted(lm(value ~ 1, data = area_pooled)), residuals(lm(value ~ 1, data = area_pooled)),
             xlab = "Fitted values", ylab = "Residuals", main = "Residual Plot")
        
        
        #KRUSKAL-WALLIS
        kruskal_test_area_pooled <- kruskal.test(value ~ group, data = area_pooled)
        print(kruskal_test_area_pooled)
        
        #DUNN TEST
        pairwise_dunn_test_area_pooled <- dunn.test(area_pooled$value, area_pooled$group, method = "bonferroni")
        print(pairwise_dunn_test_area_pooled)
        
        #MANN-WHITNEY
        pairwise_mannwhitney_test_area_pooled <- pairwise.wilcox.test(area_data$value, area_data$group, p.adjust.method = "bonferroni")
        print(pairwise_mannwhitney_test_area_pooled)
        
     

#Pooling MYOBRANCHES: MA + ME and LT + LTA to one condition and running the same analysis:
        
        ## NORMALITY TEST:¨
        
        # Calculate summary statistics
        myobranch_summary_table <- area_pooled %>%
                group_by(condition) %>%
                summarize(
                        Frequency = sum(!is.na(mb_area)),
                        Mean = mean(mb_area, na.rm = TRUE),
                        SD = sd(mb_area, na.rm = TRUE),
                        Min = min(mb_area, na.rm = TRUE),
                        Max = max(mb_area, na.rm = TRUE)
                )
        
        print(myobranch_summary_table)
        
        cat("Summary Statistics for Myobranch Area:\n")
        print(myobranch_summary_table)
        
        # Below is another way to do the summary, but the above is better.
        summary(area_pooled$mb_area)
        sd(area_pooled$mb_area, na.rm = TRUE)
        
        
        #NORMALITY TEST
        normality_test_MB <- shapiro.test(area_data$mb_area)
        
        cat("Normality Test for MyobranchArea:\n")
        print(normality_test_MB)
        
        ### NB: NON NORMAL DISTRIBUTION
        
        # POOLING AND STATISTICAL ANALYSIS + BOXPLOT:
      
        
        #KRUSKAL-WALLIS
        kruskal_test_area_pooled_myobranch <- kruskal.test(mb_area ~ group, data = area_pooled)
        print(kruskal_test_area_pooled_myobranch)
        
        
        #DUNN TEST
        pairwise_dunn_test_area_pooled_myobranch <- dunn.test(area_pooled$mb_area, area_pooled$group, method = "bonferroni")
        print(pairwise_dunn_test_area_pooled_myobranch)
        
        #MANN-WHITNEY
        pairwise_mannwhitney_test_area_pooled_myobranch <- pairwise.wilcox.test(area_pooled$mb_area, area_pooled$group, p.adjust.method = "bonferroni")
        print(pairwise_mannwhitney_test_area_pooled_myobranch)
        
        
        #BOXPLOT
        ggplot(area_pooled, aes(x = condition, y = mb_area/1000000, fill = condition)) +
                geom_boxplot() +
                stat_summary(fun.data = mean_sdl, geom = "errorbar", width = 0.2) +
                labs(x = "Pooled Conditions", y = "mm^2", title = "Boxplot of Myobranch Area by Pooled Condition") +
                scale_fill_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        
        #BARPLOT
        ggplot(area_pooled, aes(x = condition, y = mb_area/1000000, fill = condition)) +
                geom_bar(stat = "identity", color = "black") +
                labs(x = "Pooled Conditions", y = "mm^2", title = "Histogram of Myobranch Area by Pooled Condition") +
                scale_fill_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        
        #DENSITYPLOT
        ggplot(area_pooled, aes(x = mb_area/1000000, fill = condition)) +
                geom_density(alpha = 0.5) +
                labs(x = "mm^2", y = "Density", title = "Density Plot of Myobranch Area by Pooled Condition") +
                scale_fill_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        
        
        #SWARM PLOT
        ggplot(area_pooled, aes(x = condition, y = mb_area/1000000, color = condition)) +
                geom_point(position = position_jitterdodge()) +
                labs(x = "Pooled Conditions", y = "mm^2", title = "Swarm Plot of Myobranch Area by Pooled Condition") +
                scale_color_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        

# ANALYSIS WITH MYOBRANCHES EXCLUDED
        
        area_mb_excl <- area_data[area_data$myo_branch != "yes", ]  # Excluding values with corresponding "yes" in myo_branch
        
        area_mb_excl$group <- ifelse(area_mb_excl$condition == "ET", "ET",
                            ifelse(area_mb_excl$condition %in% c("LTA", "LT"), "LATE",
                                   ifelse(area_mb_excl$condition %in% c("MA", "ME"), "MEM", "CT")))

        normality_test <- shapiro.test(area_mb_excl$value)
        print(normality_test)
        
        means <- area_mb_excl %>%
                group_by(condition) %>%
                summarise(mean_value = mean(value))
        
           
        
        area_mb_excl_pooled <- area_mb_excl
        area_mb_excl_pooled$condition[area_mb_excl_pooled$condition %in% c("LT", "LTA")] <- "LATE"
        area_mb_excl_pooled$condition[area_mb_excl_pooled$condition %in% c("MA", "ME")] <- "MEM"
        
        area_mb_excl_pooled$condition <- factor(area_mb_excl_pooled$condition,
                                                levels = c("ET", "LATE", "MEM", "CT"),
                                                labels = c("ET", "LATE", "MEM", "CT"))
        
        ggplot(area_mb_excl_pooled, aes(x = condition, y = value/1000000, fill = condition)) +
                geom_boxplot() +
                labs(x = "Pooled Conditions", y = "mm^2", title = "Boxplot of Values by Pooled Condition") +
                scale_fill_manual(values = c("ET" = "red", "LATE" = "green", "MEM" = "blue", "CT" = "orange"))
        
        
        # QQ plot
        qqnorm(residuals(lm(value ~ 1, data = area_mb_excl)))
        qqline(residuals(lm(value ~ 1, data = area_mb_excl)))
        
        # Residual plot
        plot(fitted(lm(value ~ 1, data = area_mb_excl)), residuals(lm(value ~ 1, data = area_mb_excl)),
             xlab = "Fitted values", ylab = "Residuals", main = "Residual Plot")
        
        # Kruskal-Wallis Test
        kruskal_test_area_pooled_myobranch <- kruskal.test(value ~ group, data = area_mb_excl)
        print(kruskal_test_area_pooled_myobranch)
        
        # Dunn Test
        pairwise_dunn_test_area_pooled_myobranch <- dunn.test(area_mb_excl$value, area_mb_excl$group, method = "bonferroni")
        print(pairwise_dunn_test_area_pooled_myobranch)
        
        # Mann-Whitney Test
        pairwise_mannwhitney_test_area_pooled_myobranch <- pairwise.wilcox.test(area_mb_excl$value, area_mb_excl$group, p.adjust.method = "bonferroni")
        print(pairwise_mannwhitney_test_area_pooled_myobranch)
        
        
        

#SUMMARY FOR NUMBER + HISTOGRAM
cat("\nSummary Statistics for Number:\n")
print(number_summary)

number_data <- subset(dataset, morph_type == "number")
ggplot(number_data, aes(x = value)) +
        geom_histogram(fill = "steelblue", color = "white", bins = 20) +
        labs(title = "Histogram of Number") +
        xlab("Number") +
        ylab("Frequency")
