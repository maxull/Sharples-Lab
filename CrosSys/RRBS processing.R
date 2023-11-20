#
#
#   Align and preprocess CrosSys RRBS data
#
#


# pipeline: https://bioinformatics-core-shared-training.github.io/RNAseq-R/align-and-count.nb.html

BiocManager::install("Rsubread")

library(Rsubread)


# set working direcory to the cirectory with the fastq.gz files

setwd("/Users/maxul/Documents/Skole/Lab/CrosSys/test/")

fastq.files <- list.files(path = "./", pattern = ".fq.gz$", full.names = TRUE)
fastq.files


# build index file

buildindex()


align(index = "./gencode.v38.chr_patch_hapl_scaff.basic.annotation.gtf.gz", 
      readfile1 = fastq.files, 
      output_format = "BAM", 
      nthreads = 8, 
      useAnnotation = TRUE,
      annot.ext = "./gencode.v38.chr_patch_hapl_scaff.basic.annotation.gtf.gz",
      isGTF = TRUE)



































