#
# this will be first practise script for first basic session i R
#



library(tidyverse)
library(ggplot2)
library(readxl)


getwd()

# example: not tidy dataformat
dexa_data <- read_excel("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/DEXA/DEXA_data.xlsx")

read_excel("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/DEXA/DEXA_data.xlsx")

dexa_data %>% 
        select(FP) -> participants

dexa_data %>% 
        filter(FP != "MACS_011") -> dexa_data


x <- 1:10
y = "daniel"

