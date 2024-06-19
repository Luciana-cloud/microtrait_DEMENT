# Packages ####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(stats)
library(corrr)
library(olsrr)
library(vegan)
library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
# devtools::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)
# install.packages("remotes")
# remotes::install_github("gavinsimpson/ggvegan")
library(ggvegan)
library(fmsb)
library(readxl)

#

mat_ori    = read.csv("/pub/lucianac/r_work_3/hmm_MAG_IMG.csv",dec=".")
df_2       = read_excel("/pub/lucianac/r_work_3/IMG_bindata_withmeta.xlsx")
mat_ori.1  = mat_ori[mat_ori$id %in% df_2$Bin.ID, ]
df_2.1     = df_2[df_2$Bin.ID %in% mat_ori.1$id, ]

set.seed(1)
#distance.total  = vegdist(mat_ori.1[3:1723], method = "jaccard", binary = TRUE)
#cluster.total   = hclust(distance.total, method="ward.D2")
nmds            = metaMDS(mat_ori.1[3:1723], distance = "jaccard",k = 2,trymax=100)
fort.1          = fortify(nmds)
write.csv(fort.1, file = "/pub/lucianac/r_work_3/total.fort.1.csv")

