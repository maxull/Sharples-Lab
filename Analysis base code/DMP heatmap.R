################################################################
###                                                          ###
### DMP based heatmap                                        ###
###                                                          ###
################################################################



################################################################################

### DMP heatmap

library(heatmaply)
library(pheatmap)
library(tidyverse)

setwd("/Users/maxul/Documents/Coding/Sharples-Lab/")

mat <- as.matrix(head(f_bVals, n = 100000))

heatmaply(mat,
          dendrogram = "none",
          showticklabels = FALSE)


heatmap(mat, Colv = NA, Rowv = NA)    # heatmap without clustering works, but give no info...


# create heatmap of dmps

dmp_list <- rownames(myDMP$Baseline_to_7wk_Loading)


dmp_bvals <- f_bVals[dmp_list,]   # filter f_bVals for only significant DMPs

df_bVals <- head(dmp_bvals, n = 1000)


pheatmap(df_bVals, scale = "row") #scale across columns or row, scaling across rows highlights the differences

# add annotation of condition

dfh <- data.frame(sample = as.character(colnames(dmp_bvals)), condition = "Timepoint") %>% 
        column_to_rownames("sample")
dfh$condition <- ifelse(rownames(dfh) == pd_gt$colnames, pd_gt$Sample_Group, "")

dfh$condition <- factor(dfh$condition, levels = c("Baseline", "7wk_Loading"))


pheatmap(t(dmp_bvals), annotation_row = dfh, cutree_rows = 2, 
         show_colnames = FALSE, annotation_names_row = FALSE, 
         color = viridis(n = 100), scale = "row")

