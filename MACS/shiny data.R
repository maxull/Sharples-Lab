
# shiny data

# prepare data for shiny app

# Sample data for illustration purposes. Replace with your actual datasets.
myo_data <- DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1))


myo_int_data <- DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1))



type1 <- c("ANKRD2",
           "ATP2A2",
           "CA3",
           "CASQ2",
           "CD36",
           "CYB5R1",
           "FABP3",
           "LDHB",
           "MYH7",
           "MYL12A",
           "MYL2",
           "MYL3",
           "MYL6B",
           "MYOZ2",
           "PDLIM1",
           "PLN",
           "TNNC1",
           "TNNI1",
           "TNNT1",
           "TPM3"
)

type2a <- c("ALDOA",
            "ATP2A1",
            "DDIT4L",
            "ENO3",
            "G0S2",
            "GAPDH",
            "LDHA",
            "MYBPC2",
            "MYH1",
            "MYH2",
            "MYL1",
            "MYLPF",
            "PFKM",
            "PGM1",
            "PKM",
            "SLN",
            "TNNC2",
            "TNNI2",
            "TNNT3",
            "TPM1"
)

Myonuclei = Myonuclei %>% pull(symbol)

gene_lists <- list(
        type1 = type1,
        type2a = type2a,
        Myonuclei = Myonuclei
)

GO_MYOINT <- read.csv("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/GO_all_pathways_homogenate.csv")
GO_MYO <- read.csv("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/GO_all_pathways_myonuclei.csv")
KEGG_MYOINT <- read.csv("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/kegg_all_pathways_homogenate.csv")
KEGG_MYO <- read.csv("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/kegg_all_pathways_myonuclei.csv")

GO_MYOINT <- GO_MYOINT %>% filter(P.DE < 0.05) %>% dplyr::select(2:9)
GO_MYO <- GO_MYO %>% filter(P.DE < 0.05) %>% dplyr::select(2:9)
KEGG_MYO <- KEGG_MYO %>% filter(P.DE < 0.05) %>% dplyr::select(2:8)
KEGG_MYOINT <- KEGG_MYOINT %>% filter(P.DE < 0.05) %>% dplyr::select(2:8)

# add subset to KEGG pathways

names(KEGG_new$kg.sets) %>% as.data.frame() %>% 
        dplyr::select("pathway" = 1) -> subset

subset[KEGG_new$sig.idx,"subset"] <- "Signalling"
subset[KEGG_new$met.idx,"subset"] <- "Metabolic"
subset[KEGG_new$dise.idx,"subset"] <- "Disease"

KEGG_MYO <- KEGG_MYO %>% merge(., subset, by = "pathway")
KEGG_MYOINT <- KEGG_MYOINT %>% merge(., subset, by = "pathway")

# save to directory

setwd("/Users/maxul/Documents/Coding/Sharples-Lab/MACS/Shiny/")

write.csv(myo_data, "./myo_data.csv")
write.csv(myo_int_data, "./myo_int_data.csv")
write.csv(GO_MYO, "./GO_MYO.csv")
write.csv(GO_MYOINT, "./GO_MYOINT.csv")
write.csv(KEGG_MYO, "./KEGG_MYO.csv")
write.csv(KEGG_MYOINT, "./KEGG_MYOINT.csv")
saveRDS(gene_lists, "./gene_lists.RData")