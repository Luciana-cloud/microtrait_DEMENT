library(dplyr)
library(tidyverse)
library(ggplot2)
library(stats)
#install.packages("corrr")
library(corrr)
#install.packages("olsrr")
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

total.mags                   = read.csv("/pub/lucianac/r_work_3/total.mags.csv",dec=".")

set.seed(1)
distance.total  = vegdist(total.mags[3:507], method = "jaccard", binary = TRUE)
cluster.total   = hclust(distance.total, method="ward.D2")

v                      = cutree(cluster.total,k=110)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_1               = pairwiseAdonis::pairwise.adonis(distance.total,genome2guild$guild,
                                                         perm = 999,p.adjust.m='BH')
adonis_1

v                      = cutree(cluster.total,k=120)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_1               = pairwiseAdonis::pairwise.adonis(distance.total,genome2guild$guild,
                                                         perm = 999,p.adjust.m='BH')
adonis_1

v                      = cutree(cluster.total,k=130)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_1               = pairwiseAdonis::pairwise.adonis(distance.total,genome2guild$guild,
                                                         perm = 999,p.adjust.m='BH')
adonis_1

v                      = cutree(cluster.total,k=140)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_1               = pairwiseAdonis::pairwise.adonis(distance.total,genome2guild$guild,
                                                         perm = 999,p.adjust.m='BH')
adonis_1

