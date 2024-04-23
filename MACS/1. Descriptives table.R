############################################################################################
###                                                                                      ###
###  MACS: descriptiva table                                                             ###
###                                                                                      ###
############################################################################################

# make a table with baseline and post measurements, and with paired t.tests for difference

library(tidyverse)
library(readxl)
library(lubridate)
library(gtsummary)
library(gt)



# suggestetable

# pheno table"
# ----------------------------------------------------------------------
#         Measure      Baseline       Post     % Change     P.value
# ----------------------------------------------------------------------
#         Strength
#            30
#            60
#            90
#         dexa
#            total
#            leg
#            % fat 
#         CSA
#            RF
#            QF
# -----------------------------------------------------------------------


############################################################################################
###########       load strength data        ################################################
############################################################################################

setwd("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/ISOM/")

#list my isom files
isom_list <- list.files("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/ISOM/")


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

timepoints = c("Familiarization","Baseline","Post")

isom_table <- isom_df %>% 
        filter(timepoint != 1) %>% 
        group_by(timepoint, angle, reps) %>% 
        summarise(mean = mean(t_max),
                  sd = sd(t_max)) %>% 
        ungroup() %>% 
        mutate(timepoint = ifelse(timepoint == 2, "Baseline", "Post"),
               muscle = ifelse(reps == 1, "curl", "ext")) %>% 
        filter(muscle == "ext") %>% 
        pivot_wider(names_from = timepoint, values_from = c(mean, sd)) %>% 
        mutate(percent_change = ((mean_Post-mean_Baseline)/mean_Baseline)*100, 
               measure = paste0("Angle ", angle)) %>%
        dplyr::select(measure, 4:8) %>% 
        as.data.frame()


# paired t.tests

isom_test_df <- isom_df %>% 
        filter(reps == 3,
               timepoint != 1) %>% 
        dplyr::select(FP, angle, timepoint, t_max) %>% 
        mutate(timepoint = ifelse(timepoint == 2, "Baseline", "Post")) %>% 
        pivot_wider(names_from = timepoint, values_from = t_max) 
        

isom_30 <- t.test(isom_test_df %>% filter(angle == 30) %>% pull(Baseline), 
       isom_test_df %>% filter(angle == 30) %>% pull(Post), 
       paired = TRUE)

isom_60 <- t.test(isom_test_df %>% filter(angle == 60) %>% pull(Baseline), 
                  isom_test_df %>% filter(angle == 60) %>% pull(Post), 
                  paired = TRUE)

isom_90 <- t.test(isom_test_df %>% filter(angle == 90) %>% pull(Baseline), 
                  isom_test_df %>% filter(angle == 90) %>% pull(Post), 
                  paired = TRUE)


# merge together

isom_table <- isom_table %>% 
        mutate(p.value = c(isom_30$p.value,
                           isom_60$p.value,
                           isom_90$p.value))


############################################################################################
###########       load dexa data            ################################################
############################################################################################

dexa_data <- read_excel("C:/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/DEXA/DEXA_data.xlsx")

FP_list <- dexa_data %>% filter(Timepoint == "Post") %>% pull(FP)

total <- dexa_data %>% 
        filter(FP %in% FP_list) %>% 
        group_by(Timepoint) %>% 
        summarise(mean = (mean(`Mager(g)_Total`))/1000,
                  sd = (sd(`Mager(g)_Total`))/1000) %>% 
        pivot_wider(names_from = Timepoint, values_from = c(mean, sd)) %>% 
        mutate(measure = "Total lean mass")

legs <- dexa_data %>% 
        filter(FP %in% FP_list) %>% 
        group_by(Timepoint) %>% 
        summarise(mean = (mean(`Mager(g)_Legs`))/1000,
                  sd = (sd(`Mager(g)_Legs`))/1000) %>% 
        pivot_wider(names_from = Timepoint, values_from = c(mean, sd)) %>% 
        mutate(measure = "Leg lean mass")

fat_percent <- dexa_data %>% 
        filter(FP %in% FP_list) %>% 
        group_by(Timepoint) %>% 
        summarise(mean = (mean(`Fett(g)_Total`))/1000,
                  sd = (sd(`Fett(g)_Total`))/1000) %>% 
        pivot_wider(names_from = Timepoint, values_from = c(mean, sd)) %>% 
        mutate(measure = "Fat %")




# paired t-test 

# t test for change from baseline

df <- dexa_data %>% 
        filter(FP %in% FP_list) %>% 
        dplyr::select(1,2,6:15) %>% 
        pivot_wider(names_from = Timepoint, values_from = 3:12)

total_t.test <- t.test(df$`Mager(g)_Total_Baseline`,df$`Mager(g)_Total_Post`, paired = TRUE)

legs_t.test <- t.test(df$`Mager(g)_Legs_Baseline`, df$`Mager(g)_Legs_Post`, paired = TRUE)

fat_t.test <- t.test(df$`Fett(g)_Total_Baseline`, df$`Fett(g)_Total_Post`, paired = TRUE)

# merge together
dexa_table <- rbind(total, legs, fat_percent) %>% 
        dplyr::select(5,1:4) %>% 
        mutate(percent_change = ((mean_Post-mean_Baseline)/mean_Baseline)*100,
               p.value = c(round(total_t.test$p.value, 4),
                           round(legs_t.test$p.value, 5),
                           round(fat_t.test$p.value, 2)))




############################################################################################
###########       load CSA data             ################################################
############################################################################################

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


csa_table <- df %>% 
        group_by(FP, muscle, timepoint) %>% 
        summarise("cm2" = mean(cm2)) %>% 
        ungroup() %>%
        as.data.frame() %>% 
        filter(timepoint != "T1") %>% 
        group_by(muscle, timepoint) %>% 
        summarise(mean = mean(cm2),
                  sd = sd(cm2)) %>% 
        ungroup() %>% 
        mutate(timepoint = ifelse(timepoint == "T2", "Baseline", "Post")) %>% 
        pivot_wider(names_from = timepoint, values_from = c(mean, sd)) %>% 
        mutate(percent_change = ((mean_Post-mean_Baseline)/mean_Baseline)*100) %>% 
        dplyr::select(measure = muscle, 2:6)


# t.tests

csa_test_df <- df %>% 
        group_by(FP, muscle, timepoint) %>% 
        summarise("cm2" = mean(cm2)) %>% 
        ungroup() %>%
        as.data.frame() %>% 
        filter(timepoint != "T1") %>% 
        mutate(timepoint = ifelse(timepoint == "T2", "Baseline", "Post")) %>% 
        pivot_wider(names_from = timepoint, values_from = cm2)

csa_RF <- t.test(csa_test_df %>% filter(muscle == "RF") %>% pull(Baseline),
                 csa_test_df %>% filter(muscle == "RF") %>% pull(Post),
                 paired = TRUE)

csa_VL <- t.test(csa_test_df %>% filter(muscle == "VL") %>% pull(Baseline),
                 csa_test_df %>% filter(muscle == "VL") %>% pull(Post),
                 paired = TRUE)


# merge

csa_table <- csa_table %>% 
        mutate(p.value = c(round(csa_RF$p.value, 4),
                           round(csa_VL$p.value, 5)))




############################################################################################
###########       merge all tables          ################################################
############################################################################################


pheno_table <- rbind(isom_table, dexa_table, csa_table)

pheno_table <- pheno_table %>%
        mutate(
                Baseline = paste0(format(round(mean_Baseline, 1), nsmall = 0, big.mark = ","), " ± ", format(round(sd_Baseline,1), nsmall = 0, big.mark = ",")),
                Post = paste0(format(round(mean_Post, 1), nsmall = 1, big.mark = ","), " ± ", format(round(sd_Post, 1), nsmall = 1, big.mark = ",")),
                percent_change = percent_change/100, 
                Group = as.factor(c("Isometric strength","Isometric strength","Isometric strength",
                           "DEXA","DEXA","DEXA",
                           "CSA", "CSA"))
        )

### Step 2: Create the Table with gt

gt_table <- pheno_table %>%
        select(measure, Baseline, Post, percent_change, p.value, Group) %>%
        gt(groupname_col = "Group") %>%
        cols_label(
                measure = "Measure",
                Baseline = "Baseline (Mean ± SD)",
                Post = "Post (Mean ± SD)",
                percent_change = "Percent Change",
                p.value = "P Value") %>%
        fmt_percent(
                columns = c(percent_change),
                decimals = 2) %>%
        fmt_number(
                columns = c(p.value),
                decimals = 3) %>%
        tab_header(
                title = "Table 1. Summary of Phenotypical Measures") %>%
        tab_source_note(
                source_note = "Note: Data are presented as mean ± SD. Percent changes are shown as percentages, Isom strength as Nm and P values are from paired t-tests.") %>%
        tab_options(
                column_labels.font.size = px(14), column_labels.font.weight = "bold",
                row_group.font.size = px(14), 
                column_labels.border.bottom.width = 3,column_labels.border.bottom.color = "black",
                column_labels.border.top.width = 3, column_labels.border.top.color = "black",
                table.border.top.style = "none", table.border.bottom.style = "none",
                table_body.hlines.style = "none", 
                table_body.border.bottom.width = 3, table_body.border.bottom.color = "black", 
                heading.align = "left", data_row.padding = "none", 
                row_group.border.top.style = "none",row_group.border.bottom.style = "none", row_group.font.weight = "bold", row_group.padding = "none")

# Print the table
print(gt_table)



gtsave(gt_table, filename = "Pheno_table.pdf", path = "/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Figures/")


               