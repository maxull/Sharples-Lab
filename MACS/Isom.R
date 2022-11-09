#----------------------------------------------------------------
#
#   Analysis and plotting ISOM data from MACS project
#
#
#----------------------------------------------------------------


library(ggplot2); library(tidyverse); library(readxl);library(cowplot); library(tidyr)

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
        mutate(t_max = t_max)

# recode date to familiarization, baseline, post

isom_df %>% 
        mutate(day = substr(date, 1,2),
               month = substr(date, 4,5)) %>% 
        group_by(FP) %>% 
        mutate(timepoint = if(((month == min(month)) & (day == min(day)))){print("Familiarization")})
                       


# plot torque max grouped by participant, date, angle for extension

isom_df %>% 
        filter(reps == "3" & angle == 30) %>% 
        ggplot(aes(date, t_max, group = FP, color = FP))+
        geom_point()+
        geom_line()


# to do

# recode dates to timepoints fam, baseline, post
# plot all angles both extension and flection





