# CRAN package dependencies:
library("R.utils")
library("RColorBrewer")
#install.packages("ape")
library("ape")
library("assertthat")
library("checkmate")
#install.packages("coRdon")
#library("coRdon")
library("corrplot")
library("doParallel")
library("dplyr")
library("futile.logger")
library("grid")
library("gtools")
#install.packages("kmed")
library("kmed")
library("lazyeval")
library("magrittr")
library("parallel")
#install.packages("pheatmap")
library("pheatmap")
library("readr")
library("stringr")
library("tibble")
#install.packages("tictoc")
library("tictoc")
library("tidyr")
library("devtools")

# Bioconductor package dependencies:

library("BiocManager")
#BiocManager::install("Biostrings")
#BiocManager::install("coRdon")
#BiocManager::install("ComplexHeatmap")
library("Biostrings")
library("coRdon")
library("ComplexHeatmap")

# Microtrait

library("gRodon")
library("microtrait")

genomes_files = list.files('/pub/lucianac/isolates_2')
message("Number of cores:", parallel::detectCores(), "\n")

tictoc::tic.clearlog()
tictoc::tic(paste0("Running microtrait for ", length(genomes_files)))
microtrait_results = extract.traits.parallel(genomes_files, dirname(genomes_files), ncores = floor(parallel::detectCores()*0.7))
tictoc::toc(log = "TRUE")

