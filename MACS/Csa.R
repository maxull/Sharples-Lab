#----------------------------------------------------------------
#
#   Analysis and plotting CSA data from MACS project
#
#
#----------------------------------------------------------------


library(ggplot2); library(tidyverse); library(readxl);library(cowplot); library(lubridate)

getwd()

# load all sheet names

sheets <- excel_sheets("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/CSA/MACS_CSA.xlsx")

# create dataframe of all csa data

csa_df <- list()

for (i in 1:length(sheets)) {
        x <- read_excel("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/CSA/MACS_CSA.xlsx", sheet = i)
        csa_df[[i]] <- as.data.frame(x)
        print(i)
}


df <- bind_rows(csa_df) %>% 
        mutate(acquisition_date = dmy(acquisition_date))


# plot absolute change

df %>% 
        group_by(FP, timepoint, muscle) %>% 
        summarise(mean = mean(cm2),
                  sd = sd(cm2)) %>% 
        ggplot(aes(x = timepoint, y = mean, group = FP, color = FP))+
        geom_point()+
        geom_line(aes(group = FP))+
        geom_errorbar(aes(ymin = mean-sd,
                          ymax = mean+sd), width = 0.2)+
        facet_grid(~muscle)
                

# plot % change

df2 <- df %>% 
        group_by(FP, timepoint, muscle) %>% 
        summarise(mean = mean(cm2)) %>% 
        pivot_wider(names_from = timepoint, values_from = mean) %>% 
        mutate(baseline = (T1+T2)/2,
               change = T3-baseline,
               percent_change = (change/baseline)*100) %>% 
        pivot_longer(cols = 6:8, names_to = "comparison", values_to = "change") %>% 
        filter(comparison =="percent_change")

mean_df2 <- df2 %>% 
        group_by(muscle) %>% 
        summarise(mean = mean(change),
                  sd = sd(change))

ggplot(data = mean_df2, aes(x = muscle, y = mean/100))+
        geom_point(size = 2)+
        geom_errorbar(data = mean_df2, aes(ymin = ((mean/100)-(sd/100)),
                                           ymax = ((mean/100)+(sd/100))), width = 0.2)+
        geom_text(data = mean_df2, aes(label = paste0(format(round(mean, 2), nsmall = 2),"%")), fontface = "bold", hjust = -0.3)+
        geom_point(data = df2, aes(y = change/100, color = FP), size = 2)+
        theme_classic(base_size = 15)+
        scale_y_continuous(,
                           limits = c(0,0.3),
                           expand = c(0,0))+
        geom_hline(yintercept = 0, alpha = 0.5, size = 1)+
        theme(axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.text.x = element_text(size = 15, vjust = 8),
              plot.title = element_text(size = 15),
              legend.position = "none")+
        labs(y = "% change",
             title = "CSA")
        












################################################################################

### plotted MACS_001 data

csa_data <- read_excel("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/CSA/MACS_CSA.xlsx")


### plot RF

csa_data %>%
        filter(muscle == "RF") %>% 
        ggplot(aes(x = timepoint, y = cm2))+ 
        geom_boxplot()+
        expand_limits(y = 0)

csa_data %>%
        filter(muscle == "VL") %>% 
        ggplot(aes(x = timepoint, y = cm2))+ 
        geom_boxplot()+
        geom_point()+
        expand_limits(y = 0)

#### create mean and sd line and point chart for participant MACS_001
        
p1 <- csa_data %>%
        filter(FP == "MACS_001") %>% 
        group_by(timepoint, muscle) %>% 
        summarise(mean = mean(cm2),
                  sd = sd(cm2)) %>% 
        filter(muscle == "RF") %>% 
        ggplot(aes(x = timepoint, y = mean, group = muscle))+
        geom_errorbar(aes(ymin = (mean-sd),
                          ymax = (mean+sd)),
                      width = 0.2)+
        geom_point()+
        geom_line()+
        theme_classic()+
        labs(y = "cm^2 (mean ± sd)",
             x = "",
             title = "MACS_001 RF")+
        theme(axis.title.y = element_text(face = "bold"))
        
p2 <- csa_data %>%
        filter(FP == "MACS_001") %>% 
        group_by(timepoint, muscle) %>% 
        summarise(mean = mean(cm2),
                  sd = sd(cm2)) %>% 
        filter(muscle == "VL") %>% 
        ggplot(aes(x = timepoint, y = mean, group = muscle))+
        geom_errorbar(aes(ymin = (mean-sd),
                          ymax = (mean+sd)),
                      width = 0.2)+
        geom_point()+
        geom_line()+
        theme_classic()+
        labs(y = "",
             x = "",
             title = "MACS_001 VL")

plot_grid(p1, p2, ncol = 2)





