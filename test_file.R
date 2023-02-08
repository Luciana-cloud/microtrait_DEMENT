# Packages #

library(dplyr)
library(tidyverse)
library(ggplot2)

# Calling data #

mat_trait  = read.csv("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/litter_mags_trait_matrixatgranularity.csv",dec=".")
trait_sd   = mat_trait %>% select(1,43:61) # substrate degradation
trait_tp   = mat_trait %>% select(1:42) # transporters or substrate uptake
trait_st   = mat_trait %>% select(1,162:190) # stress related traits
  
# Exploratory plots #

# complex compounds
boxplot(trait_sd%>% select(2:7),ylab="gene counts",
        names=c("cellulose","chitin", "heteromannan", "linkage-glucan", "xylan", "xyloglucan"))
# protein + complex compounds
boxplot(trait_sd%>% select(2,6,8),ylab="gene counts",
        names=c("cellulose","xylan","protein"))
# DEMENT-related transporters
boxplot(trait_tp%>% select(8,4,7,18,24:27,31:35,38),ylab="gene counts",
        names=c("monosacc","amino sugar","carb_phos",
                "f_aminoac","amine","ammonium","nitrate","nitrite",
                "nucleobase","nucleoside","nucleotide","ribonucle",
                "organP","peptide"))
# DEMENT-related stress-related traits
boxplot(trait_st%>% select(13:16,27,28),ylab="gene counts",
        names=c("sol_transp","sol_synt","EPS biosyn (S)","osmo_sensing",
                "O sensing","misfolded proteins"))
# Heat-related traits
boxplot(trait_st%>% select(5:7),ylab="gene counts",
        names=c("heat shock proteins","ATP proteases","transcription factor"))

# Targeted traits

trait_tar = mat_trait %>% select(4,6:8,24:27,38,43:49,163,174:177,93:99,124,125) # 

# Forming functional guilds (https://uc-r.github.io/kmeans_clustering#fn:scale)

library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

distance <- get_dist(trait_tar)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

trait_tar1 <- (trait_tar %>% select(2:31))
trait_tar1 <- na.omit(trait_tar1)
trait_tar1 <- scale(trait_tar1)
trait_tar1 <- (trait_tar1 %>% select(1:11,13))

k2 <- kmeans(trait_tar1, centers = 100, nstart = 25)
str(k2)

