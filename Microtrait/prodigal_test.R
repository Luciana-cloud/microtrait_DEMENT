library("BiocManager")
#install.packages("devtools")
library("Biostrings")
#BiocManager::install("coRdon")
library("coRdon")
library("matrixStats")
#install.packages("matrixStats")
#remotes::install_github("jlw-ecoevo/gRodon2")
library("gRodon")
#install.packages('micropan')
library(micropan)

# Run prodigal----

genomes_files = list.files('C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/isolates/batch_1/batch_1a', pattern='fna')
gene_paths    = paste0('C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/isolates/batch_1/batch_1a/',"2636416112.fna")

#genomes_files = list.files('C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/Burn severity MAGs', pattern='fa')
#gene_paths    = paste0('C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/Burn severity MAGs/',genomes_files[1])
genes         = readDNAStringSet(gene_paths)
highly_expressed = grepl("ribosomal protein",names(genes),ignore.case = T)
predictGrowth(genes, highly_expressed)

# Subset your sequences to those that code for proteins
CDS_IDs <- readLines("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/isolates/2636416112.txt")
gene_IDs <- gsub(" .*","",names(genes)) #Just look at first part of name before the space
genes <- genes[gene_IDs %in% CDS_IDs]

#Search for genes annotated as ribosomal proteins
highly_expressed <- grepl("^(?!.*(methyl|hydroxy)).*0S ribosomal protein",names(genes),ignore.case = T, perl = TRUE)



