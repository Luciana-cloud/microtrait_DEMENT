# Packages #

library(dplyr)
library(tidyverse)
library(ggplot2)

# Calling data #

mat_ori    = read.csv("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/litter_mags_trait_matrixatgranularity.csv",dec=".")
gen_size   = read.delim("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/litter_mags_metadata.txt",dec=".") 

# Normalized data
mat_trait    = mat_ori %>% select(2:191)/gen_size$length
mat_trait    = as.data.frame(cbind(mat_ori$id,mat_trait))
trait_sd     = mat_trait %>% select(1,43:61) # substrate degradation
trait_tp     = mat_trait %>% select(1:42) # transporters or substrate uptake
trait_st     = mat_trait %>% select(1,162:190) # stress related traits

# Raw data
mat_trait_1  = mat_ori 
trait_sd_1   = mat_trait %>% select(1,43:61) # substrate degradation
trait_tp_1   = mat_trait %>% select(1:42) # transporters or substrate uptake
trait_st_1   = mat_trait %>% select(1,162:190) # stress related traits

# Exploratory plots - gene costs (normalized genes) #

# complex compounds
boxplot(trait_sd%>% select(2:7),ylab="gene cost",
        names=c("cellulose","chitin", "heteromannan", "linkage-glucan", "xylan", "xyloglucan"))
# protein + complex compounds
boxplot(trait_sd%>% select(2,6,8),ylab="gene cost",
        names=c("cellulose","xylan","protein"))
# DEMENT-related transporters
boxplot(trait_tp%>% select(8,4,7,18,24:27,31:35,38),ylab="gene cost",
        names=c("monosacc","amino sugar","carb_phos",
                "f_aminoac","amine","ammonium","nitrate","nitrite",
                "nucleobase","nucleoside","nucleotide","ribonucle",
                "organP","peptide"))
# DEMENT-related stress-related traits
boxplot(trait_st%>% select(13:16,28),ylab="gene cost",
        names=c("sol_transp","sol_synt","EPS biosyn (S)","osmo_sensing",
                "misfolded proteins"))
# Heat-related traits
boxplot(trait_st%>% select(5:7),ylab="gene cost",
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

# Presence/absence

distance  = vegdist(trait_tar_bin_1, method = "chisq", binary = TRUE)
distance1 = vegdist(trait_tar_bin_1, method = "jaccard", binary = TRUE)
# chisq = avgdist(distance, dmethod = "chisq",sample = 10000)
nguilds = seq(2, nrow(trait_tar_bin_1), 2)
cluster = hclust(distance, method="ward")
plot(cluster)

my_vec = c()

for(i in nguilds) {
  v                      = cutree(cluster,k=i)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_2               = vegan::adonis2(distance ~ guild, data = genome2guild, perm = 1)
  temp                   = adonis_2$R2
  temp_1                 = temp[1]
  my_out                 = temp_1        
  my_vec <- c(my_vec, my_out)    
}

mat_r2  = as.data.frame(cbind(nguilds,my_vec))
ggplot(data=mat_r2,aes(x=nguilds,y=my_vec)) + geom_line() +
       xlab("# of guilds") + ylab("Similarity within guilds") +
       theme(text = element_text(size = 20)) + 
       geom_hline(yintercept=0.70782750, linetype="dashed", color = "red") + 
       geom_vline(xintercept = 100, linetype="dashed", color = "red")

# Working with 100 guilds
v                      = cutree(cluster,k=100)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
mat_trait_f     = as.data.frame(cbind(genome2guild,trait_tar))  
trait_tar_bin_f = as.data.frame(cbind(genome2guild,trait_tar_bin_2))
temp = mat_trait_f %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat = (as.matrix(temp[2:35]))
heatmap(mat, Colv = NA, Rowv = NA, scale="column")
my_colnames2 <- names(temp)

# Gene cost

distanceg  = vegdist(trait_tar, method = "bray", binary = FALSE)
clusterg   = hclust(distanceg, method="ward")
plot(clusterg)

my_vec = c()

for(i in nguilds) {
  v                      = cutree(clusterg,k=i)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_2               = vegan::adonis2(distance ~ guild, data = genome2guild, perm = 1)
  temp                   = adonis_2$R2
  temp_1                 = temp[1]
  my_out                 = temp_1        
  my_vec <- c(my_vec, my_out)    
}

mat_r2_c  = as.data.frame(cbind(nguilds,my_vec))
ggplot(data=mat_r2_c,aes(x=nguilds,y=my_vec)) + geom_line() +
  xlab("# of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 20)) + 
  geom_hline(yintercept=0.7023060, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 200, linetype="dashed", color = "red")

# Working with 100 guilds
v                      = cutree(clusterg,k=100)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
mat_trait_f     = as.data.frame(cbind(genome2guild,trait_tar))  
trait_tar_bin_f = as.data.frame(cbind(genome2guild,trait_tar_bin_2))
temp = mat_trait_f %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat = (as.matrix(temp[2:35]))
heatmap(mat, Colv = NA, Rowv = NA, scale="column")
my_colnames2 <- names(temp)

# Working with 200 guilds
v                      = cutree(clusterg,k=200)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
mat_trait_f     = as.data.frame(cbind(genome2guild,trait_tar))  
trait_tar_bin_f = as.data.frame(cbind(genome2guild,trait_tar_bin_2))
temp = mat_trait_f %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat = (as.matrix(temp[2:35]))
heatmap(mat, Colv = NA, Rowv = NA, scale="column")
my_colnames2 <- names(temp)

##################################################################v




