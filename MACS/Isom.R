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
        

p2 <- 
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
        

p3 <- isom_df %>% 
        filter(reps == "3" & angle == "90") %>% 
        ggplot(aes(x = as.factor(timepoint), y = t_max, group = FP, color = FP))+
        geom_point()+
        geom_line()+
        scale_x_discrete(labels=timepoints)+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        labs(title = "90° peak torque (ext)",
             y = "Torque (Nm)")

p4 <- isom_df %>% 
        filter(reps == "1" & angle == "30") %>% 
        ggplot(aes(x = as.factor(timepoint), y = t_max, group = FP, color = FP))+
        geom_point(show.legend = FALSE)+
        geom_line(show.legend = FALSE)+
        scale_x_discrete(labels=timepoints)+
        theme(axis.title.x = element_blank())+
        labs(title = "30° peak torque (flex)",
             y = "Torque (Nm)")

p5 <- isom_df %>% 
        filter(reps == "1" & angle == "60") %>% 
        ggplot(aes(x = as.factor(timepoint), y = t_max, group = FP, color = FP))+
        geom_point(show.legend = FALSE)+
        geom_line(show.legend = FALSE)+
        scale_x_discrete(labels=timepoints)+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        labs(title = "60° peak torque (flex)",
             y = "Torque (Nm)")

p6 <- isom_df %>% 
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


# i like this plot stile

isom_df %>% 
        filter(reps == "3" & angle == "30") %>% 
        ggplot(aes(x = as.factor(timepoint), y = t_max))+
        geom_point(position = position_dodge(width = 0.1),
                    aes(color = FP))+
        geom_line(aes(group = FP, color = FP),
                  alpha = 0.5, position = position_dodge(width = 0.1))+
        geom_errorbar(aes(ymin = mean(t_max) - sd(t_max), ymax = mean(t_max) + sd(t_max)), width = .2)+
        geom_point(aes(y = mean(t_max)),
                   size = 3)+
        geom_line(aes(y = mean(t_max), group = mean(t_max)))+
        scale_x_discrete(labels = timepoints)+
        #theme_classic()+
        theme(axis.title.x = element_blank())+
        labs(title = "30° peak torque (ext)",
             y = "Torque (Nm)")
