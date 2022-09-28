##################################################################
###
### epigenetic age (horvath pan tissue)
###
###################################################################

library(wateRmelon); library(methylumi); library(minfi);; library(tidyverse);library(ggplot2)


### MsetExProbes dataset from "Oshlack workflow - filtering and QC"


bVals <- getBeta(MsetExProbes)


epiage_horvath <- agep(bVals, coeff=NULL, method="horvath")



### add sample names to epiage dataframe

df <- as.data.frame(MsetExProbes$Sample_Name)

ea <- cbind(df, epiage_horvath)

df2 <- as.data.frame(MsetExProbes$Sample_Group)

epiage <- cbind(df2, ea)

epiage <- epiage %>% 
        mutate(Group = `MsetExProbes$Sample_Group`,
               Sample = `MsetExProbes$Sample_Name`) %>% 
        select(Sample, Group, horvath.age)


### boxplot of ages for the groups

epiage %>% 
        ggplot(aes(x = Group, y = horvath.age, fill = Group))+
        geom_boxplot()+
        theme_classic()






################################################################################
################################################################################
################################################################################
###
###     MEAT 2.0 muscle tissue clock
###
################################################################################
################################################################################
################################################################################


library(SummarizedExperiment); library(MEAT)



### get bValues as data frame

bVals <- read.csv("C:/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/Cancer data/Cancer data normalised B_values for clock analysis_transposed.csv")


### if col 1 is the CpGs use code below

samp <- bVals[,-1]
rownames(samp2) <- bVals[,1]

### create Summarized experiment element for the epiage_estimation 
### as seen below it can be run with and without a pheno dataframe

cancer <- SummarizedExperiment(assays = list(beta = samp),
                               colData = pheno)                 

cancer2 <- SummarizedExperiment(assays = list(beta = samp))


### then the CpGs need to be "cleaned" so the dataset only contains the relevant 18747 CpG sites

cancer_clean <- clean_beta(SE = cancer,
                           version = "MEAT2.0")

### calibrate beta values 



cancer_clean_calibrated <- BMIQcalibration(cancer_clean)


### then you estimate the epigenetic age with this function, where the "age_col_name" is optional, 
###     but reqired if you wish to get difference and residuals between chronological age and predicted age


epiage_meat <- epiage_estimation(SE = cancer_clean_calibrated,
                                 age_col_name = NULL,
                                 version = "MEAT2.0")

### get only the age estimations


DNAmage <- as.data.frame(epiage_meat$DNAmage)

### add to earlyer created participant data with horvath clock estimations



epiage <- cbind(epiage, DNAmage)


epiage_data <- epiage %>% 
        select(1,2,3,6,7,8,9,10)

write.csv(epiage_data, "C:/Users/maxul/OneDrive/Dokumenter/Skole/Master 21-22/Master/Cancer data/cancer_epiage.csv",
          row.names = FALSE)







