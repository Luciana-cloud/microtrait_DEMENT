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
