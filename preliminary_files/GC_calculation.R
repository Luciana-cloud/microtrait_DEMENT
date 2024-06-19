# Calculate GC content

# Instaling packages ----
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Biostrings")
library(Biostrings)
#install.packages("metaCluster", 
#                 repos = c("https://diprosinha.r-universe.dev", 
#                           "https://cran.r-project.org"))
library(metaCluster)
library(seqinr)

# Calling data ----

filenames   = list.files("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/seqtk_MEGAHITcoa_bins", 
                        pattern="*.fa", full.names=TRUE)
sample_data = read.fasta(file = system.file(filenames[[1]], package = "metaCluster"),
                          seqtype = "DNA")
gc          = GC.content(genomes_files[1])

legacyfile <- system.file(filenames[[1]], package = "metaCluster")
legacyseq <- read.fasta(file = filenames[[1]], as.string = TRUE)
gc          = GC.content(legacyseq[["k141_2440_flag_1_multi_6.0000_len_47"]])

