#####################################################################################################
###
###     epiage for cancer data (thormods cancer project)
###
#####################################################################################################


BiocManager::install("SummarizedExperiment")
BiocManager::install("MEAT")

library(SummarizedExperiment); library(MEAT); library(wateRmelon); library(ggplot2); library(tidyverse)



###########################
####
#### cancer data


bVals <- read.csv("C:/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/Cancer data/Cancer data normalised B_values for clock analysis_transposed.csv")


### col 1 changed to row names


samp2 <- bVals[,-1]
rownames(samp2) <- bVals[,1]



### estimate horwath clock

epiage_horvath <- agep(samp2, coeff=NULL, method="horvath")



### read participant data

library(readxl);library(tibble)
df <- tibble::rownames_to_column(epiage_horvath, "ID")

participants <- read_xlsx("C:/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/Cancer data/array no & conditions.xlsx")


### combine horvath clock and participant information


epiage <- cbind(df, participants)

epiage <- epiage %>% 
        select(c(1,2,3,5,6,7))



############################3
###
### running MEAT 2.0 on the same data



pheno <- read_xlsx("C:/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/Cancer data/Phenotypes.xlsx")

cancer <- SummarizedExperiment(assays = list(beta = samp2),
                               colData = pheno)                 ### add chronological age her if you want AAdiff and AAresid, and then add the name in "age_col_name = ..."

cancer2 <- SummarizedExperiment(assays = list(beta = samp2))


### clean beta matwix so it only contains the 18747 elements (CpGs) needed

cancer_clean <- clean_beta(SE = cancer,
                           version = "MEAT2.0")

cancer2_clean <- clean_beta(SE = cancer2,
                            version = "MEAT2.0")



### BMIQ calibration 


calibrated <- BMIQcalibration(cancer2_clean)



### estimate epigenetic age

epiage_meat <- epiage_estimation(SE = cancer2_clean,
                                 age_col_name = NULL,
                                 version = "MEAT2.0")

epiage_meat <- epiage_estimation(SE = calibrated,
                                 age_col_name = NULL,
                                 version = "MEAT2.0")

### add the epi age estimation to the epiage dataframe

DNAmage <- as.data.frame(epiage_meat$DNAmage)

DNAmage$`epiage_meat$DNAmage`= DNAmage$`epiage_meat$DNAmage`+20

epiage <- cbind(epiage, DNAmage)


### fix dataset for visualization


epiage<- epiage %>% 
        mutate(Timepoint = factor(TIMEPOINT, levels = c("Pre", "Post")),
               Condition = factor(CONDITION, levels = c("Healthy Aged-matched Trained", "Cancer Untrained", "Cancer Trained")),
               DNAmage = `epiage_meat$DNAmage`)


### boxplot of epiage estimations

epi_plot <- epiage %>% 
        ggplot(aes(x = Condition, y = DNAmage, fill = Timepoint, color = Timepoint))+
        geom_boxplot(fatten = 1,lwd = 1.5)+
        theme_classic()+
        scale_fill_manual(values=c("#97aad2", "#daceb8"))+
        scale_color_manual(values=c("#041334", "#6d5b37"))+
        theme(legend.key.size = unit(1, 'cm'),
              legend.key.height = unit(1, 'cm'),
              legend.key.width = unit(1, 'cm'))+
        theme(axis.title.x  = element_blank(),
              axis.text.x = element_text(size = 14, face = "bold"),
              axis.title.y = element_text(size = 18,face = "bold"))+
        ylab("DNAm Age")+
        ylim(47,70)

### get boxplot data

?ggplot_build(epi_plot)$data %>% 
        as.data.frame() %>% 
        mutate(group = c("Healthy Age-matched PRE",
                         "Healthy Age-matched POST",
                         "Cancer Untrained PRE",
                         "Cancer Untrained POST",
                         "Cancer Trained PRE",
                         "Cancer Trained POST")) %>% 
        dplyr::select(group,
                      3:10) %>% 
        mutate(outliers = ifelse(outliers == 0, "NA", outliers)) %>% 
        as.data.frame() -> fig_9_data

write.csv(fig_9_data,"C:/Users/maxul/Documents/Skole/Master 21-22/Master/Cancer data/fig_9_data.csv")

fig_9_data$outliers <- as.character(fig_9_data$outliers)

# create legends for manual copying
epiage %>% 
        ggplot(aes(x = DNAmage, fill = Timepoint, color = Timepoint))+
        geom_histogram(size = 1.15)+
        theme_classic()+
        scale_fill_manual(values=c("#97aad2", "#daceb8"))+
        scale_color_manual(values=c("#041334", "#6d5b37"))+
        theme(legend.key.size = unit(5, 'cm'),
              legend.key.height = unit(0.7, 'cm'),
              legend.key.width = unit(0.7, 'cm'))
fig_9_data$


# colors

# fill pre = #97aad2
# border pre = #041334
# fill post = #daceb8
# border post = #6d5b37



########################################
##########################################3
#######################################


### clean dataset and save

epiage_data <- epiage %>% 
        select(1,2,3,6,7,8,9,10)

write.csv(epiage_data, "C:/Users/maxul/Documents/Skole/Master 21-22/Master/Cancer data/Cancer data/cancer_epiage_corrected_mistake.csv",
          row.names = FALSE)
write.csv(epiage_data, "C:/Users/maxul/Documents/Skole/Master 21-22/Master/Cancer data/epiage_data.csv")

###
### AA diff


epiage_data <- epiage_data %>% 
        mutate(AAdiff = ifelse(epiage_data$Condition == "Cancer Trained", 60.4 - epiage_data$DNAmage, 
                               ifelse(epiage_data$Condition == "Cancer Untrained", 62 - epiage_data$DNAmage,
                                      58.3 - epiage_data$DNAmage)))

epiage_data %>% 
        ggplot(aes(x = Condition, y = AAdiff, fill = Timepoint))+
        geom_boxplot()+
        theme_classic()



#################
###
### comparing means of uncalibrated and calibrated signals

library(reshape2)

calibrated_df = as.data.frame(calibrated@assays@data@listData[["beta"]])

precalibrated_df = as.data.frame(cancer2_clean@assays@data@listData[["beta"]])

epiage$all <- paste(epiage$FP,epiage$CONDITION, epiage$TIMEPOINT, "precal", sep = "_")

names(calibrated_df) <- as.list(epiage$all)

names(precalibrated_df) <- as.list(epiage$all)


barplot(colMeans(calibrated_df))

barplot(colMeans(precalibrated_df))

### changed col names to include all fp info, continue with making boxplot from here



calibrated_df$CpG <- row.names(calibrated_df)

df_c <- melt(calibrated_df)



df_c %>% 
        ggplot(aes(x = value, color = variable))+
        geom_density()+
        labs(title = "all 18747 CpGs")

df_c %>% 
        ggplot(aes(x = variable, y = value))+
        geom_boxplot()


df_uc <- melt(precalibrated_df)

df_uc %>% 
        ggplot(aes(x = value, color = variable))+
        geom_density()+
        labs(title = "all 18747 CpGs")

df_uc %>% 
        ggplot(aes(x = variable, y = value))+
        geom_boxplot()+
        theme_classic()

df <- cbind(x = df_c, y = df_uc)

df %>% 
        ggplot(aes(x = x.value, y = y.value))+
        geom_point()+
        theme_classic()+
        xlab("calibrated")+
        ylab("precalibrated")+
        labs(title = "all 18747 CpGs")

melt_df <- melt(df)
melt_df <- melt_df %>% 
        mutate(x.variable = as.character(x.variable))

melt_df$x.variable = substr(melt_df$x.variable, 1, nchar(melt_df$x.variable)-4)

melt_df %>% 
        ggplot(aes(x = x.variable, y = value, color = variable))+
        geom_boxplot()+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.title.x=element_blank())+
        scale_color_discrete(name = "BMIQcalibration",  labels = c("Calibrated", "Pre-Calibrated"))+
        ylab("Beta Value")+
        labs(title = "all 18747 CpGs")

### works

### chech only the CpG sites used in MEAT 2.0

### load list of 200 cpgs used in MEAT2.0

cpg <- read_xlsx("C:/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/Cancer data/200cpg.xlsx")


### keep only 200 cpgs in df og calibrated and precalibrated data

cpg <- cpg %>% 
        select(1)

calibrated_df1 <- calibrated_df %>% 
        mutate(CpG = row.names(calibrated_df))

merge_cal <- as.data.frame(cpg %>% left_join(calibrated_df1, by = "CpG"))

merge_cal <- merge_cal %>% remove_rownames %>% column_to_rownames(var="CpG")

melt_cal_200 <- melt(merge_cal)

melt_cal_200 %>% 
        ggplot(aes(x = value, color = variable))+
        geom_density()+
        labs(title = "200 CpGs in MEAT2.0")


precalibrated_df1 <- precalibrated_df %>% 
        mutate(CpG = row.names(precalibrated_df))

merge_precal <- as.data.frame(cpg %>% left_join(precalibrated_df1, by = "CpG"))

merge_precal <- merge_precal %>% remove_rownames %>% column_to_rownames(var="CpG")

melt_precal_200 <- melt(merge_precal)

melt_precal_200 %>% 
        ggplot(aes(x = value, color = variable))+
        geom_density()+
        labs(title = "200 CpGs in MEAT2.0")


df_200 <- cbind(x = melt_cal_200, y = melt_precal_200)

df_200 %>% 
        ggplot(aes(x = x.value, y = y.value))+
        geom_point()+
        theme_classic()+
        xlab("calibrated")+
        ylab("precalibrated")+
        labs(title = "200 CpGs in MEAT2.0")

melt_df_200 <- melt(df_200)
melt_df_200 <- melt_df_200 %>% 
        mutate(x.variable = as.character(x.variable))

melt_df_200$x.variable = substr(melt_df_200$x.variable, 1, nchar(melt_df_200$x.variable)-4)

melt_df_200 %>% 
        ggplot(aes(x = x.variable, y = value, color = variable))+
        geom_boxplot()+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.title.x=element_blank())+
        scale_color_discrete(name = "BMIQcalibration",  labels = c("Calibrated", "Pre-Calibrated"))+
        ylab("Beta Value")+
        labs(title = "200 CpGs in MEAT2.0")



###################################33
###
### grouped boxplots of calibrated and precalibrated
library(cowplot); library(grid); library(gridExtra)

all_df <- melt_df %>% 
        mutate(x.variable = factor(gsub("^[^_]*_", "", melt_df$x.variable),
                                   levels = c("Healthy Aged-matched Trained_Pre",
                                              "Healthy Aged-matched Trained_Post",
                                              "Cancer Untrained_Pre",
                                              "Cancer Untrained_Post",
                                              "Cancer Trained_Pre",
                                              "Cancer Trained_Post")))

cpg_200_df <- melt_df_200 %>% 
        mutate(x.variable = factor(gsub("^[^_]*_", "", melt_df$x.variable),
                                   levels = c("Healthy Aged-matched Trained_Pre",
                                              "Healthy Aged-matched Trained_Post",
                                              "Cancer Untrained_Pre",
                                              "Cancer Untrained_Post",
                                              "Cancer Trained_Pre",
                                              "Cancer Trained_Post")))


all_df %>% 
        ggplot(aes(x = x.variable, y = value, color = variable))+
        geom_boxplot()+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.title.x=element_blank())+
        scale_color_discrete(name = "BMIQcalibration",  labels = c("Calibrated", "Pre-Calibrated"))+
        ylab("Beta Value")+
        labs(title = "all 18747 CpGs")

cpg_200_df %>% 
        ggplot(aes(x = x.variable, y = value, color = variable))+
        geom_boxplot()+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.title.x=element_blank())+
        scale_color_discrete(name = "BMIQcalibration",  labels = c("Calibrated", "Pre-Calibrated"))+
        ylab("Beta Value")+
        labs(title = "200 CpGs in MEAT2.0")



filtered_calibrated_df <- calibrated_df %>% 
        select(2,3,7,9,10)

f_c_df <- melt(filtered_calibrated_df)

f_c_df %>% 
        ggplot(aes(x = value, color = variable))+
        geom_density()


filtered_precalibrated_df <- precalibrated_df %>% 
        select(2,3,7,9,10)

f_pc_df <- melt(filtered_precalibrated_df)

f_pc_df %>% 
        ggplot(aes(x = value, color = variable))+
        geom_density()


#############################################################################
###########################################################################################
###########################################################################
###
###     stat tests between groups



library(readr)
library(ggplot2)
library(multcompView)
library(dplyr)


epiage_data$Timepoint <- as.factor(epiage_data$Timepoint)
epiage_data$Condition <- as.factor(epiage_data$Condition)

epiage_data$DNAmage <- epiage_data$DNAmage + 20

anova <- aov(DNAmage ~ Condition * Timepoint, data = epiage_data)
summary(anova)

data_summary <- group_by(epiage_data, Condition, Timepoint) %>%
        summarise(mean=mean(DNAmage), sd=sd(DNAmage)) %>%
        arrange(desc(mean))
print(data_summary)

tukey <- TukeyHSD(anova)
print(tukey)

epiage_data %>%
        select(4, 6,7,8) %>% 
        pivot_wider(names_from = Timepoint, values_from = DNAmage) %>% 
        na.omit() %>% 
        mutate(diff = Post - Pre) -> diff

diff_summary <- group_by(epiage_data, Condition) %>%
        summarise(mean=mean(diff), sd=sd(diff)) %>%
        arrange(desc(mean))
print(data_summary)





