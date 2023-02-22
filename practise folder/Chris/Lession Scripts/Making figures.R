

# Building Figures


library(tidyverse)
library(readxl)

cyclingstudy%>%
        select(subject, timepoint, group, sj.max)%>%
        ggplot(aes(x = timepoint, y = sj.max)) + geom_point()

# ggplot is used to visualize most data.
# (aes) deterimen what variables are assigned to the X and Y axis.
# geom_point maps data to the plot just made.


#FACTOR() - enables mutation of variable names, useful to give new descriptions to the data in the plot.
# We want pre, meso1-2-3 to be more "laymans", we first select the levels (i.e., our variables to be renamed), 
# then we set labels to determine the new names. \n is used to row break.

cyclingstudy %>%
        select(subject, timepoint, group, sj.max) %>%
        mutate(timepoint = factor(timepoint, 
                                  levels = c("pre", "meso1", "meso2", "meso3"),
                                  labels = c("Prior to\ntrainining\nperiod", "Period 1", "Period 2", "Period 3"))) %>%
        ggplot(aes(x = timepoint, y = sj.max)) + geom_point()

#This can also be done with box-plots or violin-plots as well

cyclingstudy %>%
        select(subject, timepoint, group, sj.max) %>%
        mutate(timepoint = factor(timepoint, 
                                  levels = c("pre", "meso1", "meso2", "meso3"),
                                  labels = c("Prior to\ntrainining\nperiod", "Period 1", "Period 2", "Period 3"))) %>%
        ggplot(aes(x = timepoint, y = sj.max)) + geom_boxplot() #or geom_violin


#------

#Lets say we want to plot the mean of each group, in each timepoint, in the dataset.
#Since the mean is a summary statistic, we must first make the summary before we can plot it.
#We can keep our factored mutation as this only change the names of the variables.
#thus the summarise function falls after the mutation.
# 

cyclingstudy %>%
        select(subject, group, timepoint, group, sj.max) %>%
        
        mutate(timepoint = factor(timepoint, 
                                  levels = c("pre", "meso1", "meso2", "meso3"),
                                  labels = c("Prior to\ntrainining\nperiod", "Period 1", "Period 2", "Period 3"))) %>%
        group_by(timepoint, group)%>% #i.e., make the means of each group at each timepoint.
                summarise(m.sj = mean(sj.max, na.rm = TRUE))%>% #summarise (name of summary = function( na.rm removes N/A))
                print()%>%
        
                ggplot(aes(timepoint, m.sj)) + geom_point()

#The plot does not tell us which dot belongs to what group, only the timepoint.
# aestetics can fix this 


cyclingstudy %>%
        select(subject, group, timepoint, group, sj.max) %>%
        
        mutate(timepoint = factor(timepoint, 
                                  levels = c("pre", "meso1", "meso2", "meso3"),
                                  labels = c("Prior to\ntrainining\nperiod", "Period 1", "Period 2", "Period 3"))) %>%
        group_by(timepoint, group)%>% #i.e., make the means of each group at each timepoint.
        summarise(m.sj = mean(sj.max, na.rm = TRUE))%>% #summarise (name of summary = function( na.rm removes N/A))
        print()%>%
        
        ggplot(aes(timepoint, m.sj, color = group)) + geom_point()

# now we can see each group being given their own color category!

# we can from here play around by adding lines within groups across timepoints, by geom_line()
# however, merely adding + geom_line() will give an error, this is because we must tell ggplot what grouping we want.
# so, assigning group = group, solves this, but why Im not sure...
# We can also add Standard Deviation, however, its calculation MUST come before the mean as to not calc SD from means.
# [size] is a function that set the size of both geom_point and line to be relative to the given variable.

cyclingstudy %>%
        select(subject, group, timepoint, group, sj.max) %>%
        
        mutate(timepoint = factor(timepoint, 
                                  levels = c("pre", "meso1", "meso2", "meso3"),
                                  labels = c("Prior to\ntrainining\nperiod", "Period 1", "Period 2", "Period 3"))) %>%
        group_by(timepoint, group)%>% #i.e., make the means of each group at each timepoint.
        summarise(sd.sj.max = sd(sj.max, na.rm = TRUE), 
                  m.sj = mean(sj.max, na.rm = TRUE))%>% #summarise (name of summary = function( na.rm removes N/A))
        print()%>%
        
        ggplot(aes(timepoint, m.sj, color = group, group = group, size = sd.sj.max)) + geom_point() + geom_line()


# This gives us a weird graph, so to remove the size from within the lines we set size = NULL in geomline

cyclingstudy %>%
        select(subject, group, timepoint, group, sj.max) %>%
        
        mutate(timepoint = factor(timepoint, 
                                  levels = c("pre", "meso1", "meso2", "meso3"),
                                  labels = c("Prior to\ntrainining\nperiod", "Period 1", "Period 2", "Period 3"))) %>%
        group_by(timepoint, group)%>% #i.e., make the means of each group at each timepoint.
        summarise(sd.sj.max = sd(sj.max, na.rm = TRUE), 
                  m.sj = mean(sj.max, na.rm = TRUE))%>% #summarise (name of summary = function( na.rm removes N/A))
        print()%>%
        
        ggplot(aes(timepoint, m.sj, color = group, group = group, size = sd.sj.max)) +
                geom_point() + 
                geom_line(aes(size = NULL))

# This now set each point size to be relative to the SD. 
# Importantly, to set a characteristic of either ggplot itself, or any geom function, aes() must be used.

# Typically we don't use size to represent SD, but rather error bars!¨
# geom_errorbar allows this


cyclingstudy %>%
        select(subject, group, timepoint, group, sj.max) %>%
        
        mutate(timepoint = factor(timepoint, 
                                  levels = c("pre", "meso1", "meso2", "meso3"),
                                  labels = c("Prior to\ntrainining\nperiod", "Period 1", "Period 2", "Period 3"))) %>%
        group_by(timepoint, group)%>% #i.e., make the means of each group at each timepoint.
        summarise(sd.sj.max = sd(sj.max, na.rm = TRUE), 
                  m.sj = mean(sj.max, na.rm = TRUE))%>% #summarise (name of summary = function( na.rm removes N/A))
        print()%>%
        
        ggplot(aes(timepoint, m.sj, color = group, group = group)) +
                geom_errorbar(aes(ymin = m.sj - sd.sj.max, 
                                  ymax = m.sj + sd.sj.max)) + 
                                #here we set the min/max range of the error bars, SD +/- the mean gives this range!
        geom_point() + geom_line()

#-------------
# Nice, but the plot is messy as fuck. We can clean it up by altering the positions of its elements

cyclingstudy %>%
        select(subject, group, timepoint, group, sj.max) %>%
        
        mutate(timepoint = factor(timepoint, 
                                  levels = c("pre", "meso1", "meso2", "meso3"),
                                  labels = c("Prior to\ntrainining\nperiod", "Period 1", "Period 2", "Period 3"))) %>%
        group_by(timepoint, group)%>% #i.e., make the means of each group at each timepoint.
        summarise(sd.sj.max = sd(sj.max, na.rm = TRUE), 
                  m.sj = mean(sj.max, na.rm = TRUE))%>% #summarise (name of summary = function( na.rm removes N/A))
        
        ggplot(aes(timepoint, m.sj, color = group, group = group)) +
        geom_errorbar(aes(ymin = m.sj - sd.sj.max, 
                          ymax = m.sj + sd.sj.max), 
                position = position_dodge(width = 0.2),width = 0.2) + 
                        geom_point(position = position_dodge(width = 0.2)) +
                        geom_line(position = position_dodge(width = 0.2))

#position_dodge, nudges the groups to the side by the given value.
# width determine the width of the errorbar hats

cyclingstudy %>%
        select(subject, group, timepoint, group, sj.max) %>%
        
        mutate(timepoint = factor(timepoint, 
                                  levels = c("pre", "meso1", "meso2", "meso3"),
                                  labels = c("Prior to\ntrainining\nperiod", "Period 1", "Period 2", "Period 3"))) %>%
        group_by(timepoint, group) %>%
        summarise(sd.sj.max = sd(sj.max, na.rm = TRUE),
                  m.sj = mean(sj.max, na.rm = TRUE)) %>%
        
        ggplot(aes(timepoint, m.sj, color = group, group = group))  + 
        
        geom_errorbar(aes(ymin = m.sj - sd.sj.max, 
                          ymax = m.sj + sd.sj.max), 
                      position = position_dodge(width = 0.2), 
                      width = 0.1) +
        
        geom_point(position = position_dodge(width = 0.2))  +
        geom_line(position = position_dodge(width = 0.2)) +
        
        labs(x = "Time-point", y = "Squat jump (cm)", color = "Group",
                title = "Squat jump height per group",
                subtitle = "Different patterns per group",
                caption = "Values are mean and SD") + 
        
        theme_bw()+
                theme(axis.title = element_text(size = 14) +
                theme (plot.background = element_rect(fill = "gray60")))


#Finally, we label our variables using the labs() function
# -- and we can add titles, subtitles and captions within the label function
# we can also set a theme that format our plot neatly, whereby we can change the formating extensively


#----------------------------------------------

# To save this file we must have the full script


library(tidyverse)%>%
library(readxl)%>%
        
Figure1 <- cyclingstudy %>%
        select(subject, group, timepoint, group, sj.max) %>%
        
        mutate(timepoint = factor(timepoint, 
                                  levels = c("pre", "meso1", "meso2", "meso3"),
                                  labels = c("Prior to\ntrainining\nperiod", "Period 1", "Period 2", "Period 3"))) %>%
        group_by(timepoint, group) %>%
        summarise(sd.sj.max = sd(sj.max, na.rm = TRUE),
                  m.sj = mean(sj.max, na.rm = TRUE)) %>%
        
        ggplot(aes(timepoint, m.sj, color = group, group = group))  + 
        
        geom_errorbar(aes(ymin = m.sj - sd.sj.max, 
                          ymax = m.sj + sd.sj.max), 
                      position = position_dodge(width = 0.2), 
                      width = 0.1) +
        
        geom_point(position = position_dodge(width = 0.2))  +
        geom_line(position = position_dodge(width = 0.2)) +
        
        labs(x = "Time-point", y = "Squat jump (cm)", color = "Group",
             title = "Squat jump height per group",
             subtitle = "Different patterns per group",
             caption = "Values are mean and SD") + 
        
        theme_bw()+
        theme(axis.title = element_text(size = 14) +
                      theme (plot.background = element_rect(fill = "gray60")))%>%

        ggsave("practise folder/Chris/figures/figure1.png", plot = Figure1) #didnt work for some reason..


#-----------------------------------------------------------


        















