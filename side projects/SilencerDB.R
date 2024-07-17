¢###################################################################################################################
##################       Gene silencer methylation                   ##############################################
###################################################################################################################  

# load full proba info and location
Illumina_anno <- read.csv("Annotation file Illumina/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7, header = TRUE)

Illumina_anno <- Illumina_anno %>% 
        mutate(CHR = paste0("chr", CHR))

Illumina_anno <- Illumina_anno %>% 
        mutate(MAPINFO = as.numeric(MAPINFO))

# load silencer database
Silencer_DB <- read_tsv("Silencer_DB_Homo_sapiens_Muscle_UCSC.bed", skip = 1)

#fix colnames of silencer_db

Silencer_DB <- rbind(Silencer_DB,colnames(Silencer_DB))

Silencer_DB <- Silencer_DB[,c(1:7,10)]

colnames(Silencer_DB) <- c("Chr", "Start", "End", "Silencer_ID", "Length", "Strand", "Muscle_origin", "Detection_method")

str(Silencer_DB)

Silencer_DB <- Silencer_DB %>% 
        mutate(Start = as.numeric(Start),
               End = as.numeric(End))


####
#       Identify which illumina probes are located in silencer regions
####

# Initialize the 'Within_Silencer' column
Illumina_anno$Within_Silencer <- NA

# Loop through each row in the Illumina_anno data frame
for (i in 1:nrow(Illumina_anno)) {
        # Filter silencers on the same chromosome
        relevant_silencers <- Silencer_DB %>%
                filter(Chr == as.character(Illumina_anno$CHR[i]))
        
        # Check each silencer for overlap
        for (j in 1:nrow(relevant_silencers)) {
                if (Illumina_anno$MAPINFO[i] >= relevant_silencers$Start[j] &&
                    Illumina_anno$MAPINFO[i] <= relevant_silencers$End[j]) {
                        Illumina_anno$Within_Silencer[i] <- relevant_silencers$Silencer_ID[j]
                        
                }
        }
        print(i)
}

# second iteration

# Setup a parallel cluster
cl <- makeCluster(detectCores() - 1)  # Leave one core free for system processes
clusterExport(cl, list("Illumina_anno", "Silencer_DB"))  # Export necessary objects to each worker

# Run pblapply in parallel
results <- pblapply(seq_len(nrow(Illumina_anno)), function(i) {
        relevant_silencers <- Silencer_DB[Silencer_DB$Chr == as.character(Illumina_anno$CHR[i]), ]
        for (j in seq_len(nrow(relevant_silencers))) {
                if (Illumina_anno$MAPINFO[i] >= relevant_silencers$Start[j] &&
                    Illumina_anno$MAPINFO[i] <= relevant_silencers$End[j]) {
                        return(relevant_silencers$Silencer_ID[j])
                }
        }
        return(NA)
}, cl = cl)

# Stop the cluster
stopCluster(cl)

# unslit results

results <- unlist(results)

results <- as.data.frame(results)


Illumina_anno <- cbind(Illumina_anno, results) %>% 
        dplyr::select(1:51, "Silencer" = results)

# save annotation file

write_csv(Illumina_anno, "Annotation file Illumina/infinium-methylationepic-v-1-0-b5-manifest-file_with_silencers.csv")

# check file

silencer_anno <- Illumina_anno %>% 
        filter(Silencer != "NA") %>% 
        dplyr::select("IlmnID","CHR", "MAPINFO", "Strand", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Group", "Silencer_ID" = Silencer) %>% 
        merge(., Silencer_DB, by = "Silencer_ID") %>% 
        dplyr::select("IlmnID","CHR", "MAPINFO", "Strand.x", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Group", "Silencer_ID", 10:16)

write_csv(silencer_anno, "Annotation file Illumina/infinium-methylationepic-v-1-0-b5_silencer_anno.csv")