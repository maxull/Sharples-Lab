#----------------------------------------------------------------
#
#   Analysis and plotting ISOM data from MACS project
#
#
#----------------------------------------------------------------


library(ggplot2); library(tidyverse); library(readxl);library(cowplot); library(tidyr); library(lubridate); library(ggpubr)

# get a list of files in the end folder


isom_list <- list.files("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/ISOM/")

#get list of csv files in woring directory, so you have to add filepath or set working directory to the folder of interest
#temp = list.files(pattern="*.csv")

wd <- getwd()
setwd("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/ISOM/")
getwd()
setwd(wd)


myfiles = lapply(isom_list, read.csv)


#create dataframe with file names which contain test info

col_names <- c("FP", "date", "time", "angle", "reps", "leg", "set")

isom_df <- as.data.frame(isom_list) %>% 
        separate(isom_list, sep = " ", into = col_names) %>% 
        mutate(FP = substr(FP, 1, 8),
               set = substr(set, 1,1)) 


#extract peak torque from myfiles
t_max <- list()

for (i in 1: length(myfiles)){
        t_max[[i]] = max(myfiles[[i]]$Torque.Newton.Meters.)
}

# now add the string of torque max data to a column in isom_df

isom_df <- isom_df %>% 
        mutate(t_max = as.numeric(t_max))


#change date to date format

isom_df <- isom_df %>% 
        mutate(date = dmy(date))


# recode each participants unique date to timepoint


        
for (i in unique(isom_df$FP)) {

a <-unique(isom_df$date[isom_df$FP==i])
        
for (l in 1:length(a)) {
        logical <- isom_df$FP == i & isom_df$date == a[l]
        isom_df$timepoint[logical] = l
        
}


}  


################################################################################3

### add a list of your timepoints

timepoints = c("Familiarization","Baseline","Post")



# plot torque max grouped by participant, date, angle for extension


#p1 <- 
isom_df %>% 
        filter(reps == "3" & angle == "30") %>% 
        ggplot(aes(x = as.factor(timepoint), y = t_max))+
        geom_boxplot()+
        geom_point(aes(group = FP, color = FP),show.legend = FALSE)+
        geom_line(aes(group = FP, color = FP),show.legend = FALSE)+
        scale_x_discrete(labels=timepoints)+
        theme_bw()+
        theme(axis.title.x = element_blank())+
        labs(title = "30° peak torque (ext)",
             y = "Torque (Nm)")
        

#p2 <- 
isom_df %>% 
        filter(reps == "3" & angle == "60") %>% 
        ggplot(aes(x = as.factor(timepoint), y = t_max))+
        geom_boxplot()+
        geom_point(aes(group = FP, color = FP),show.legend = FALSE)+
        geom_line(aes(group = FP, color = FP),show.legend = FALSE)+
                scale_x_discrete(labels=timepoints)+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        labs(title = "60° peak torque (ext)",
             y = "Torque (Nm)")
        

#p3 <- 
isom_df %>% 
        filter(reps == "3" & angle == "90") %>% 
        ggplot(aes(x = as.factor(timepoint), y = t_max, group = FP, color = FP))+
        geom_point()+
        geom_line()+
        scale_x_discrete(labels=timepoints)+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        labs(title = "90° peak torque (ext)",
             y = "Torque (Nm)")

#p4 <- 
isom_df %>% 
        filter(reps == "1" & angle == "30") %>% 
        ggplot(aes(x = as.factor(timepoint), y = t_max, group = FP, color = FP))+
        geom_point(show.legend = FALSE)+
        geom_line(show.legend = FALSE)+
        scale_x_discrete(labels=timepoints)+
        theme(axis.title.x = element_blank())+
        labs(title = "30° peak torque (flex)",
             y = "Torque (Nm)")

#p5 <- 
isom_df %>% 
        filter(reps == "1" & angle == "60") %>% 
        ggplot(aes(x = as.factor(timepoint), y = t_max, group = FP, color = FP))+
        geom_point(show.legend = FALSE)+
        geom_line(show.legend = FALSE)+
        scale_x_discrete(labels=timepoints)+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        labs(title = "60° peak torque (flex)",
             y = "Torque (Nm)")

#p6 <- 
isom_df %>% 
        filter(reps == "1" & angle == "90") %>% 
        ggplot(aes(x = as.factor(timepoint), y = t_max, group = FP, color = FP))+
        geom_point()+
        geom_line()+
        scale_x_discrete(labels=timepoints)+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        labs(title = "90° peak torque (flex)",
             y = "Torque (Nm)")

isom_plot <- ggarrange(p1,p2,p3,p4,p5,p6, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")

plot_grid(p1,p2,p3,p4,p5,p6, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")

### save the plot

ggsave("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/isom_plot.png", 
       plot = isom_plot,
       units = "px",
       dpi = 600,
       bg = "white")

############################################################################################33

### plot with facet grid instead



isom_df %>% 
        ggplot(aes(x = as.factor(timepoint), y = t_max, group = FP, color = FP))+
        geom_point()+
        geom_line()+
        scale_x_discrete(labels=timepoints)+
        labs(title = "90° peak torque (flex)",
             y = "Torque (Nm)")+
        facet_grid(reps~angle)



# plot 1 with mean and error bars


# i like this plot stile, but error in mean calculation

isom_df %>% 
        filter(reps == "3" & angle == "30") %>% 
        ggplot(aes(x = as.factor(timepoint), y = t_max))+
        geom_point(position = position_dodge(width = 0.1),
                    aes(color = FP))+
        geom_line(aes(group = FP, color = FP),
                  alpha = 0.5, position = position_dodge(width = 0.1))+
        geom_errorbar(aes(ymin = mean(t_max) - sd(t_max), ymax = mean(t_max) + sd(t_max)), width = .2)+
        geom_point(aes(y = mean(t_max), group = timepoint),
                   size = 3)+
        geom_line(aes(y = mean(t_max), group = mean(t_max)))+
        scale_x_discrete(labels = timepoints)+
        #theme_classic()+
        theme(axis.title.x = element_blank())+
        labs(title = "30° peak torque (ext)",
             y = "Torque (Nm)")



################################################################################################33
################################################################################################3
##############################################################################################33

# create dataset with filtered extension and angle

ext_30 <- isom_df %>% 
        filter(reps == "3" & angle == "30")

ext_30_m <- isom_df %>% 
        filter(reps == "3" & angle == "30") %>% 
        group_by(timepoint) %>% 
        summarize(m = mean(t_max),
                  sd = sd(t_max))

ext_60 <- isom_df %>% 
        filter(reps == "3" & angle == "60")

ext_60_m <- isom_df %>% 
        filter(reps == "3" & angle == "60") %>% 
        group_by(timepoint) %>% 
        summarize(m = mean(t_max),
                  sd = sd(t_max))

ext_90 <- isom_df %>% 
        filter(reps == "3" & angle == "90")

ext_90_m <- isom_df %>% 
        filter(reps == "3" & angle == "90") %>% 
        group_by(timepoint) %>% 
        summarize(m = mean(t_max),
                  sd = sd(t_max))

flex_30 <- isom_df %>% 
        filter(reps == "1" & angle == "30")

flex_30_m <- isom_df %>% 
        filter(reps == "1" & angle == "30") %>% 
        group_by(timepoint) %>% 
        summarize(m = mean(t_max),
                  sd = sd(t_max))

flex_60 <- isom_df %>% 
        filter(reps == "1" & angle == "60")

flex_60_m <- isom_df %>% 
        filter(reps == "1" & angle == "60") %>% 
        group_by(timepoint) %>% 
        summarize(m = mean(t_max),
                  sd = sd(t_max))

flex_90 <- isom_df %>% 
        filter(reps == "1" & angle == "90")

flex_90_m <- isom_df %>% 
        filter(reps == "1" & angle == "90") %>% 
        group_by(timepoint) %>% 
        summarize(m = mean(t_max),
                  sd = sd(t_max))



p30 <- ggplot()+
        geom_point(data = ext_30, aes(x = as.factor(timepoint), y = t_max, color = FP))+
        geom_line(data = ext_30, aes(x = timepoint, y = t_max, color = FP))+
        geom_point(data = ext_30_m, aes(x = as.factor(timepoint), y = m))+
        geom_line(data = ext_30_m, aes(x = as.factor(timepoint), y = m, group = 1))+
        geom_errorbar(data = ext_30_m, aes(x = as.factor(timepoint), ymin = m-sd, ymax = m+sd), width = 0.2)+
        scale_x_discrete(labels = timepoints)+
        theme(axis.title.x = element_blank())+
        labs(title = "30° peak torque (ext)",
             y = "Torque (Nm)")

p60 <- ggplot()+
        geom_point(data = ext_60, aes(x = as.factor(timepoint), y = t_max, color = FP))+
        geom_line(data = ext_60, aes(x = timepoint, y = t_max, color = FP))+
        geom_point(data = ext_60_m, aes(x = as.factor(timepoint), y = m))+
        geom_line(data = ext_60_m, aes(x = as.factor(timepoint), y = m, group = 1))+
        geom_errorbar(data = ext_60_m, aes(x = as.factor(timepoint), ymin = m-sd, ymax = m+sd), width = 0.2)+
        scale_x_discrete(labels = timepoints)+
        theme(axis.title.x = element_blank())+
        labs(title = "60° peak torque (ext)",
             y = "Torque (Nm)")

p90 <- ggplot()+
        geom_point(data = ext_90, aes(x = as.factor(timepoint), y = t_max, color = FP))+
        geom_line(data = ext_90, aes(x = timepoint, y = t_max, color = FP))+
        geom_point(data = ext_90_m, aes(x = as.factor(timepoint), y = m))+
        geom_line(data = ext_90_m, aes(x = as.factor(timepoint), y = m, group = 1))+
        geom_errorbar(data = ext_90_m, aes(x = as.factor(timepoint), ymin = m-sd, ymax = m+sd), width = 0.2)+
        scale_x_discrete(labels = timepoints)+
        theme(axis.title.x = element_blank())+
        labs(title = "90° peak torque (ext)",
             y = "Torque (Nm)")

p30e <- ggplot()+
        geom_point(data = flex_30, aes(x = as.factor(timepoint), y = t_max, color = FP))+
        geom_line(data = flex_30, aes(x = timepoint, y = t_max, color = FP))+
        geom_point(data = flex_30_m, aes(x = as.factor(timepoint), y = m))+
        geom_line(data = flex_30_m, aes(x = as.factor(timepoint), y = m, group = 1))+
        geom_errorbar(data = flex_30_m, aes(x = as.factor(timepoint), ymin = m-sd, ymax = m+sd), width = 0.2)+
        scale_x_discrete(labels = timepoints)+
        theme(axis.title.x = element_blank())+
        labs(title = "30° peak torque (flex)",
             y = "Torque (Nm)")

p60e <- ggplot()+
        geom_point(data = flex_60, aes(x = as.factor(timepoint), y = t_max, color = FP))+
        geom_line(data = flex_60, aes(x = timepoint, y = t_max, color = FP))+
        geom_point(data = flex_60_m, aes(x = as.factor(timepoint), y = m))+
        geom_line(data = flex_60_m, aes(x = as.factor(timepoint), y = m, group = 1))+
        geom_errorbar(data = flex_60_m, aes(x = as.factor(timepoint), ymin = m-sd, ymax = m+sd), width = 0.2)+
        scale_x_discrete(labels = timepoints)+
        theme(axis.title.x = element_blank())+
        labs(title = "60° peak torque (flex)",
             y = "Torque (Nm)")

p90e <- ggplot()+
        geom_point(data = flex_90, aes(x = as.factor(timepoint), y = t_max, color = FP))+
        geom_line(data = flex_90, aes(x = timepoint, y = t_max, color = FP))+
        geom_point(data = flex_90_m, aes(x = as.factor(timepoint), y = m))+
        geom_line(data = flex_90_m, aes(x = as.factor(timepoint), y = m, group = 1))+
        geom_errorbar(data = flex_90_m, aes(x = as.factor(timepoint), ymin = m-sd, ymax = m+sd), width = 0.2)+
        scale_x_discrete(labels = timepoints)+
        theme(axis.title.x = element_blank())+
        labs(title = "90° peak torque (flex)",
             y = "Torque (Nm)")

plot_grid(p30,p60,p90,p30e,p60e,p90e, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")

# Wrong label e = flex
ggarrange(p30,p60,p90,p30e,p60e,p90e, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")




# 3 participants had setup errors in the dynamometer, therefore have to use T1 as baseline!!! 




################################################################################################33
################################################################################################3
##############################################################################################33

# % change from baseline

df_ext <- isom_df %>%
        select(1,4,5,6,8,9) %>% 
        filter(timepoint >= 2) %>%
        mutate(timepoint = ifelse(timepoint == 2 , "baseline" , "post")) %>% 
        filter(reps == 3) %>% 
        mutate(names = paste(timepoint, angle, sep = "_")) %>% 
        select(FP, names, t_max) %>% 
        pivot_wider(names_from = names, values_from = t_max) %>% 
        mutate(change_30 = ((post_30 - baseline_30)/baseline_30)*100,
               change_60 = ((post_60 - baseline_60)/baseline_60)*100,
               change_90 = ((post_90 - baseline_90)/baseline_90)*100,)




# plot % change from baseline

df_ext_2 <- df_ext %>% 
        pivot_longer(cols = 2:10, names_to = "test", values_to = "results") %>% 
        filter(test != c("change_30", "change_60", "change_90")) %>% 
        separate(test, into = c("timepoint", "angle"))

df_ext_change <- df_ext %>% 
        pivot_longer(cols = 2:10, names_to = "test", values_to = "results") %>% 
        filter(test == c("change_30", "change_60", "change_90")) %>% 
        separate(test, into = c("timepoint", "angle"))




# % change extension plot



isom_df2 <- df_ext_2 %>% 
        pivot_wider(names_from = timepoint, values_from = results) %>% 
        mutate(change = post - baseline,
               percent_change = (change/baseline)*100) %>% 
        filter(percent_change < 80)

isom_percent <- isom_df2 %>% 
        na.omit() %>% 
        group_by(angle) %>% 
        summarise(mean_change = mean(change),
                  sd_change = sd(change),
                  mean_percent = mean(percent_change),
                  sd_percent = sd(percent_change))
        
        


ext_plot <- ggplot(data = isom_percent, aes(x = angle, y = mean_percent/100))+
        geom_point(size = 2)+
        geom_errorbar(data = isom_percent, aes(ymin = ((mean_percent - sd_percent)/100),
                                             ymax = ((mean_percent + sd_percent)/100)), width = 0.2)+
        geom_text(data = isom_percent, aes(label = paste0(format(round(mean_percent, 2), nsmall = 2),"%"), fontface = "bold", hjust = -0.3))+
        geom_point(data = isom_df2, aes(x = angle, y = percent_change/100, color = FP), size = 2)+
        theme_classic(base_size = 15)+
        scale_y_continuous(labels = percent,
                           limits = c(-0.5,1),
                           expand = c(0,0))+
        geom_hline(yintercept = 0, alpha = 0.5, linewidth = 1)+
        theme(axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.text.x = element_text(size = 15, vjust = 8),
              plot.title = element_text(size = 15),
              legend.position = "none")+
        labs(y = "% change",
             title = "ISOM extention")
                



# % flex plot

df_flex <- isom_df %>%
        select(1,4,5,6,8,9) %>% 
        filter(timepoint >= 2) %>%
        mutate(timepoint = ifelse(timepoint == 2 , "baseline" , "post")) %>% 
        filter(reps == 1) %>% 
        mutate(names = paste(timepoint, angle, sep = "_")) %>% 
        select(FP, names, t_max) %>% 
        pivot_wider(names_from = names, values_from = t_max) %>% 
        mutate(change_30 = ((post_30 - baseline_30)/baseline_30)*100,
               change_60 = ((post_60 - baseline_60)/baseline_60)*100,
               change_90 = ((post_90 - baseline_90)/baseline_90)*100,)

df_flex_2 <- df_flex %>% 
        pivot_longer(cols = 2:10, names_to = "test", values_to = "results") %>% 
        filter(test != c("change_30", "change_60", "change_90")) %>% 
        separate(test, into = c("timepoint", "angle"))



df_flex_3 <- df_flex_2 %>% 
        pivot_wider(names_from = timepoint, values_from = results) %>% 
        mutate(change = post - baseline,
               percent_change = (change/baseline)*100) %>% 
        filter(baseline > 20)           ### remove outlier

isom_percent_flex <- df_flex_3 %>% 
        na.omit() %>% 
        group_by(angle) %>% 
        summarise(mean_change = mean(change),
                  sd_change = sd(change),
                  mean_percent = mean(percent_change),
                  sd_percent = sd(percent_change))


flex_plot <- ggplot(data = isom_percent_flex, aes(x = angle, y = mean_percent/100))+
        geom_point(size = 2)+
        geom_errorbar(data = isom_percent_flex, aes(ymin = ((mean_percent - sd_percent)/100),
                                               ymax = ((mean_percent + sd_percent)/100)), width = 0.2)+
        geom_text(data = isom_percent_flex, aes(label = paste0(format(round(mean_percent, 2), nsmall = 2),"%"), fontface = "bold", hjust = -0.3))+
        geom_point(data = df_flex_3, aes(x = angle, y = percent_change/100, color = FP), size = 2)+
        theme_classic(base_size = 15)+
        scale_y_continuous(labels  = percent,
                           limits = c(-0.5,1),
                           expand = c(0,0))+
        geom_hline(yintercept = 0, alpha = 0.5, linewidth = 1)+
        theme(axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.text.x = element_text(size = 15, vjust = 8),
              plot.title = element_text(size = 15),
              legend.position = "none",
              axis.title.y = element_blank())+
        labs(title = "ISOM flection")


# plot extention and flection together

isom_plot <- plot_grid(ext_plot, flex_plot, ncol = 2)





        
df_flex_2 <- df_flex %>% 
        pivot_longer(cols = 2:10, names_to = "test", values_to = "results") %>% 
        filter(test != c("change_30", "change_60", "change_90")) %>% 
        separate(test, into = c("timepoint", "angle"))

df_flex_change <- df_flex %>% 
        pivot_longer(cols = 2:10, names_to = "test", values_to = "results") %>% 
        filter(test == c("change_30", "change_60", "change_90")) %>% 
        separate(test, into = c("timepoint", "angle"))

