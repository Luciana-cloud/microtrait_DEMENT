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

trait_tar = mat_trait %>% select(4,6:8,18,23:25,31:35,38,43:49,163,174:177,93:99,125) # 
genomes   = as.list(mat_trait[1])
test      = trait_tar %>% summarize_if(is.numeric, sum, na.rm=TRUE)

# Change to binary matrix

trait_tar_bin_1     = as.data.frame(ifelse(trait_tar>0,1,0))
trait_tar_bin_2     = as.data.frame(cbind(genomes,trait_tar_bin_1))
trait_tar_T         = as.numeric(as.matrix((trait_tar_bin_1)))

# Forming functional guilds (https://www.youtube.com/watch?v=GPOUGpF-Sno)
# https://github.com/ukaraoz/microtrait

library(vegan)
set.seed(1)
distance  = vegdist(trait_tar_bin_1, method = "chisq", binary = TRUE)
distance1 = vegdist(trait_tar_bin_1, method = "jaccard", binary = TRUE)

# chisq = avgdist(distance, dmethod = "chisq",sample = 10000)

nguilds = seq(2, nrow(trait_tar_bin_1), 2)
cluster = hclust(distance, method="ward")
plot(cluster)

# Two clusters (0.08827973)
v                      = cutree(cluster,k=2)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_2               = vegan::adonis2(distance ~ guild, data = genome2guild, perm = 1)
adonis_2$R2

# Four clusters (0.1836519)
v                      = cutree(cluster,k=4)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_4               = vegan::adonis2(distance ~ guild, data = genome2guild, perm = 1)
adonis_4$R2

# 50 clusters (0.5726445)
v                      = cutree(cluster,k=50)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_50              = vegan::adonis2(distance ~ guild, data = genome2guild, perm = 1)
adonis_50$R2

# 100 clusters (0.7078275)
v                      = cutree(cluster,k=100)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_100             = vegan::adonis2(distance ~ guild, data = genome2guild, perm = 1)
adonis_100$R2

# 120 clusters (0.7462432)
v                      = cutree(cluster,k=120)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_120             = vegan::adonis2(distance ~ guild, data = genome2guild, perm = 1)
adonis_120$R2

# 150 clusters (0.7945925)
v                      = cutree(cluster,k=150)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_150             = vegan::adonis2(distance ~ guild, data = genome2guild, perm = 1)
adonis_150$R2

# Raw plot

guild_plot <- data.frame(guild_N = c (2,4,50,100,120,150), 
              R2 = c(0.08827973,0.1836519,0.5726445,0.7078275,0.7462432,0.7945925))

ggplot(data=guild_plot,aes(x=guild_N,y=R2)) + geom_line()

# Working with 100 guilds

mat_trait_f     = as.data.frame(cbind(genome2guild,trait_tar))  
trait_tar_bin_f = as.data.frame(cbind(genome2guild,trait_tar_bin_2))

temp = mat_trait_f %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)

mat = (as.matrix(temp[2:35]))

heatmap(mat, Colv = NA, Rowv = NA, scale="column")

my_colnames2 <- names(temp)




