# Basic lessions, imporitng files

#Tidy data => One OBSERVATION per ROW, one VARIABLE per COLUMN.
#Never leave a cell empty => NA is used as identifyer for missing values.

#-------------------------


#How to read an excel file:

library(readxl)
cyclingstudy <- read_excel("practise folder/Chris/Data/cyclingStudy.xlsx", na ="NA")
View(cyclingstudy)

#First we load the library (or package) containing the function we want.
# the [ <- ] points the following function to be contained within this OBJECT.
#       Our objects is called "cyclingstudy"

# Notice that the starting of the path is any folder in the main Repository,
# I.E., within the Sharples-lab folder. ALso, remember spelling. 

# Read_excel includes the following functions:
# (path, sheet = NULL, range = NULL, col_names = TRUE, col_types = NULL, na = ""
# trim_ws =TRUE, skip = 0, n_max = inf, guess_max = min(1000, n_max), 
#  progress = readxl_progress (), .name_repair = "unique")

# ==> "Path" specify the location of the file
# ==> "Sheet" specify whether the xlsx file contain multiple sheets and how many to work with
# ==> "na" specify what to be deemed missing data, i.e., what to exclude
# ==> "skip" specify the number of rows to skip before reading the data

#-----------------------------

# How to read CSV files (comma separarated values), i.e., text type files.

library(readr)
hypertrophy1 <- read_csv("practise folder/Chris/Data/hypertrophy1.csv")

#To save a CSV file from e.g., github: Open the file, then press RAW and save as
#Read_CSV also contain several function similar to those above.

#-----------------------------


# DATA WRANGLING AND SUMMARIES

# PIPE [ %>% ] enables the operation to act sequentially.
# tidyverse is a package that enables dplyr and tidyr which allow read and visualization of data.
#if this package is not loaded the pipe function will not work!

#lets load a file to play with
# Our GOAL is to calculate squat jump height per kg body mass. Look into the dataset and learn the variable names.
#sgj.bm is our new variable, sj.max is an observed squat jump height, and weight.T1 is our observation of bodymass.

library(tidyverse)%>%
library(readxl)%>%
cyclingstudy <- read_excel("practise folder/Chris/Data/cyclingStudy.xlsx", na ="NA")

libcyclingstudy %>% #this essentially takes the dataset and implies "then do=>"
        mutate(sgj.bm = sj.max/weight.T1) %>% # Here we make our new variable, then
        select(subject, group, timepoint, sgj.bm) %>% #Here we select important variables to include when we =>
        print(.) #show our results

#For some dumbass reason, you must use the TAB to manually select the variable the function is to act upon.
# E.G., to actually select the sj.max, or subject, variables, you must tab and find it. Apparently just writing it does no work.


# FILTER and DISTINCT functions
#Distinct is a function to extract information from within a column or row. E.g.,

cyclingstudy%>%
        distinct(timepoint)%>%
        print()

#This gives the output with information about what types of values are observed in the timepoint colum
# in this case pre, meso1, meso2, meso3. 
#Now we can filter by these variables

cyclingstudy%>%
        filter (timepoint=="pre")%>%
        print()

#The output is now all observations within the PRE grouping. 

#We can also invert this function by using [ ! ]

cyclingstudy%>%
        filter(timepoint != "pre")%>%
        print()

# This then filter OUT all PRE values, leaving the rest


# If we want to keep two timepoints, the function change:

cyclingstudy%>%
        filter(timepoint %in% c("pre", "meso1"))%>%
        print()

#We cannot use == as no observations are both pre and meso1. 

# we can also filter by numerical values either equal to, less or more than a given value

cyclingstudy%>%
        filter(sj.max > 30) %>%
        print()

# we could write < or == here as well

# Finally, these could be combined!

cyclingstudy%>%
        filter(timepoint=="pre", sj.max > 30)%>%
        print()

#----------------------------------------------


# ARRANGE, sort the data by a given variable

cyclingstudy%>%
        arrange(sj.max)%>%
        print()

#the [desc] #descending" function can reverse the order of this arrangement

cyclingstudy%>%
        arrange(desc(sj.max))%>%
        print()

#this can also be done by a simple [ - ] "minus" sign.

#-------------------------

#GROUP BY and SUMMARISE
# Allows summarized data within groups!

cyclingstudy%>%
        group_by(timepoint)%>%
        summarise(m = mean(sj.max, na.rm = TRUE))%>%
        print()

#If we want to exclude NA values then we add the [ na.rm = TRUE ] function =>

cyclingstudy%>%
        group_by(timepoint)%>%
        summarise((m = mean(sj.max, na.rm = TRUE)))%>%
        print()

#This grouping can contain multiple variables, and the summarize may also expand to include e.g., standard deviation

cyclingstudy%>%
        group_by(timepoint, group)%>%
        summarize(m = mean(sj.max, na.rm =TRUE), s = sd(sj.max, na.rm = TRUE))%>%
        print()

#-----------------------------

# GROUP BY and MUTATE

# Lets calculate the squat jump height as a percent of the timepoint mean

cyclingstudy%>%
        group_by(timepoint)%>%
        mutate(sqjump = (sj.max / mean(sj.max, na.rm =TRUE))*100)%>%
        select(subject, timepoint, sj.max, sqjump)%>%
        print()

#This tells us how the squat jump fares against the mean squat jump given in a percent, 
# where 100 = perfectly equal to the mean. These values are also grouped by timepoint. 

#-----------------------------------------

# PIVOT LONGER and PIVOT WIDER 

# Data is sometimes not found in tidy form, pivot enables alterations to the dataset by means of R.
# In the current dataset, several lactate values are within the same timepiont
# Lets select the variables we need to perform a lactate thershold analysis:

cyclingstudy%>%
        select(subject, timepoint, lac.125:lac.375)%>% # [ : ] is used to express a range
        print()

#This gives a non-tidy view of the data, with each lactate observation now having their own colum.
# What we want is to have each lactate observation within their own row, with a value in their own colum.

cyclingstudy%>%
        select(subject, timepoint, lac.125:lac.375)%>%
        pivot_longer(names_to = "watt", values_to = "lactate", cols = lac.125:lac.375)%>%
        print()

#Thus, we pivot longer, which starts by making a new colum for the "names" of the lactate observations. 
#In this case, each lactate observation is done at a certain wattage, so the new variable is "WATTS". 
# Then, another colum is made for the actual observed values, which we name "Lactate" indicating the amont of lactate.

# We can clean the watt column by removing the prefix =>

cyclingstudy%>%
        select(subject, timepoint, lac.125:lac.375)%>%
        pivot_longer(names_to = "watt", values_to = "lactate", cols = lac.125:lac.375,
                     names_prefix = "lac.")%>%
        print()

# Now the watt variable is numerical, yet it is identified as a character.


cyclingstudy%>%
        select(subject, timepoint, lac.125:lac.375)%>%
        pivot_longer(names_to = "watt", values_to = "lactate", cols = lac.125:lac.375,
                     names_prefix = "lac.", names_transform = list(watt = as.numeric))%>%
        print()

#Now data is tidy.

#----------------------------------

# PIVOT WIDER - acting opposite to pivot longer

cyclingstudy%>%
        select(subject, timepoint, sj.max)%>%
        print()

#We can see that the timepoint variable could actually be given its own colums:

cyclingstudy%>%
        select(subject, timepoint, sj.max)%>%
        pivot_wider(names_from = timepoint, values_from = sj.max)%>%
        print()

#Now we can calculate percentage change from the mesocycles

cyclingstudy%>%
        select(subject, timepoint, sj.max)%>%
        pivot_wider(names_from = timepoint, values_from = sj.max)%>%
        mutate(change =((meso3/pre)-1)*100)%>%
        print()

# To save a pipe operation as the one above, we simply assign the operation to an object

Percent_change <- cyclingstudy%>%
        select(subject, timepoint, sj.max)%>%
        pivot_wider(names_from = timepoint, values_from = sj.max)%>%
        mutate( change =((meso3/pre)-1)*100)%>%
        print()

# Basic lessions cumpleted!


# TO Summarize:



