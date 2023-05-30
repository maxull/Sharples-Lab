##################################################################################
###
###  Histological analysis, MyoVision data
###
################################################################################

library(R.matlab)

MACS_002_2.3 <- readMat("C:/Users/maxul/Downloads/MACS_002_2.3_dys.mv.prj.mat")
BiocManager::install("rhdf5")


library(rhdf5)

MACS_002_2.3 <-h5read("C:/Users/maxul/Downloads/MACS_002_2.3_dys.mv.prj.mat")

MACS_002_2.3 <- h5dump("C:/Users/maxul/Downloads/MACS_002_2.3_dys.mv.prj.mat")



# load excel spread sheets

library(readxl)
library(tidyverse)
library(ggpubr)

directory <- "/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/IHC/Myovision_Mark Viggars/"

files <- list.files(directory)


ihc_data <- lapply(files, function(x) read_excel(paste0(directory, x), skip = 1))

names(ihc_data) = files

ihc_df = data.frame()

for (i in 1:length(ihc_data)) {
        
        x = str_split(strin = files[i], pattern = "_",simplify = TRUE)
        
        ihc_data[[i]] %>% 
                dplyr::select(1:5) %>% 
                mutate(FP = x[2],
                       Timepoint = ifelse(x[3] == "2.1", "Baseline", "Post")) %>% 
                rbind(ihc_df, .) -> ihc_df
        
        
        
        
}

ihc_df %>% 
        group_by(FP, Timepoint) %>% 
        summarise(mean_fCSA = mean(Fiber_CSA),
                  fiber_number = sum(Fiber_CSA > 0),
                  MHC1 = sum(FT__FT1 == "TRUE"),
                  myo_per_fiber = mean(Myonuclei)) %>% 
        mutate(MHC2 = fiber_number-MHC1) -> ihc

ihc %>% group_by(Timepoint) %>% summarise(m = mean(mean_fCSA),
                  s = sd(mean_fCSA)) -> df

p1 <- ihc %>% 
        ggplot(aes(mean_fCSA/1000, x = as.factor(Timepoint), color = FP))+
        geom_point(size = 2)+
        geom_line(aes(group = as.factor(FP)))+
        geom_point(data = df,aes(y = m/1000), color = "black", size = 3)+
        geom_errorbar(data = df, aes(ymin = (m-s)/1000,
                                     ymax = (m+s)/1000, x = Timepoint),inherit.aes = FALSE, width = 0.2)+
        theme_classic(base_size = 20)+
        theme(axis.title.x = element_blank())+
        labs(y = "mean fCSA (mm^2")


ihc %>% group_by(Timepoint) %>% summarise(m = mean(myo_per_fiber),
                                          s = sd(myo_per_fiber)) -> df2

p2 <- ihc %>% 
        ggplot(aes(myo_per_fiber, x = as.factor(Timepoint), color = FP))+
        geom_point(size = 2)+
        geom_line(aes(group = as.factor(FP)))+
        geom_point(data = df2,aes(y = m), color = "black", size = 3)+
        geom_errorbar(data = df2, aes(ymin = m-s,
                                     ymax = m+s, x = Timepoint),inherit.aes = FALSE, width = 0.2)+
        theme_classic(base_size = 20)+
        theme(axis.title.x = element_blank())+
        labs(y = "myonuclei/myofiber")



ggarrange(p1,p2, nrow = 1, labels = c("A","B"), common.legend = TRUE, legend = "right")





# replot with fibertype specific data

ihc_df %>% 
        group_by(FP, Timepoint, FT__FT1) %>% 
        summarise(mean_fCSA = mean(Fiber_CSA),
                  fiber_number = sum(Fiber_CSA > 0),
                  MHC1 = sum(FT__FT1 == "TRUE"),
                  myo_per_fiber = mean(Myonuclei)) %>% 
        mutate(MHC2 = fiber_number-MHC1,
               fiber_type = as.factor(ifelse(FT__FT1 == "TRUE", "MHC1", "MHC2"))) -> ihc2


ihc2 %>% 
        ggplot(aes(mean_fCSA, x = as.factor(Timepoint), color = FP, group = fiber_type))+
        geom_point(size = 2)+
        geom_line(aes(group = as.factor(FP)))+
        geom_point(data = df,aes(y = m), color = "black", size = 3)+
        geom_errorbar(data = df, aes(ymin = m-s,
                                     ymax = m+s, x = Timepoint),inherit.aes = FALSE, width = 0.2)+
        theme_classic()+
        theme(axis.title.x = element_blank())


# t test on fCSA

t.test(mean_fCSA~Timepoint, paired = TRUE, data = ihc) 

# Load necessary libraries
library(knitr)
library(kableExtra)

# Convert Timepoint to factor
ihc$Timepoint <- as.factor(ihc$Timepoint)

# Conduct t-test
result <- t.test(mean_fCSA~Timepoint, paired = TRUE, data = ihc)

# Calculate mean, standard deviation, and percent change for each time point
baseline <- subset(ihc, Timepoint == "Baseline")
post <- subset(ihc, Timepoint == "Post")

mean_baseline <- round(mean(baseline$mean_fCSA), 2)
sd_baseline <- round(sd(baseline$mean_fCSA), 2)

mean_post <- round(mean(post$mean_fCSA), 2)
sd_post <- round(sd(post$mean_fCSA), 2)

percent_change <- round((mean_post - mean_baseline) / mean_baseline * 100, 2)

# Create a data frame for the APA style table
apa_table <- data.frame(
        Timepoint = c("Baseline", "Post", "T-test Results"),
        Mean = c(mean_baseline, mean_post, paste("t(", result$parameter, ") =", round(result$statistic, 2))),
        SD = c(sd_baseline, sd_post, paste("p =", round(result$p.value, 4))),
        PercentChange = c("", percent_change, ""),
        stringsAsFactors = FALSE
)

# Create APA style table using kable
kable(apa_table, format = "html", col.names = c("Timepoint", "Mean", "Standard Deviation", "Percent Change"),
      align = c("l", "r", "r", "r")) %>%
        kable_styling("striped", full_width = F)

# Load necessary libraries
library(knitr)
library(kableExtra)

# Convert Timepoint to factor
ihc$Timepoint <- as.factor(ihc$Timepoint)



     
# Conduct t-test
result <- t.test(myo_per_fiber~Timepoint, paired = TRUE, data = ihc)

# Calculate mean, standard deviation, and percent change for each time point
baseline <- subset(ihc, Timepoint == "Baseline")
post <- subset(ihc, Timepoint == "Post")

mean_baseline <- round(mean(baseline$myo_per_fiber), 2)
sd_baseline <- round(sd(baseline$myo_per_fiber), 2)

mean_post <- round(mean(post$myo_per_fiber), 2)
sd_post <- round(sd(post$myo_per_fiber), 2)

percent_change <- round((mean_post - mean_baseline) / mean_baseline * 100, 2)

# Create a data frame for the APA style table
apa_table <- data.frame(
        Timepoint = c("Baseline", "Post", "T-test Results"),
        Mean = c(mean_baseline, mean_post, paste("t(", result$parameter, ") =", round(result$statistic, 2))),
        SD = c(sd_baseline, sd_post, paste("p =", round(result$p.value, 4))),
        PercentChange = c("", percent_change, ""),
        stringsAsFactors = FALSE
)

# Create APA style table using kable
kable(apa_table, format = "html", col.names = c("Timepoint", "Mean", "Standard Deviation", "Percent Change"),
      align = c("l", "r", "r", "r")) %>%
        kable_styling("striped", full_width = F)

