# Calculate GC content

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
library(Biostrings)
install.packages("metaCluster", 
                 repos = c("https://diprosinha.r-universe.dev", 
                           "https://cran.r-project.org"))
library(metaCluster)
