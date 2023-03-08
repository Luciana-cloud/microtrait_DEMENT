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

# Working with 200 guilds
v                      = cutree(cluster,k=200)
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

# Working with 200 guilds
v                      = cutree(clusterg,k=200)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
mat_trait_f     = as.data.frame(cbind(genome2guild,mat_ori$id,trait_tar))  
temp = mat_trait_f %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat = (as.matrix(temp[2:35]))
heatmap(mat)
# heatmap(mat, Colv = NA, Rowv = NA, scale="column")
# my_colnames2 <- names(temp)

# Preliminary arrangements
a = as.data.frame(colnames(temp))
r_acqui = rowSums(temp %>% select(2:22), na.rm=FALSE)
s_tol   = rowSums(temp %>% select(23:27), na.rm=FALSE)
r_use   = rowSums(temp %>% select(28:35), na.rm=FALSE)

temp1  = as.data.frame(cbind(temp$guild,r_acqui,
                             s_tol,
                             r_use))
mat1 = (as.matrix(temp1[2:4]))
heatmap(mat1)

library(fmsb)

temp2  = as.data.frame(rbind(rep(max(r_acqui),3),
                             rep(min(r_use),3),temp1[2:4]))
temp3  = as.data.frame(temp2[c(1,2,3,4),])
radarchart(temp2)

temp4  = as.data.frame(temp2[c(1,2,8),])
radarchart(temp4)

par(mfrow=c(1,3))
plot(temp1$r_acqui,temp1$s_tol)
plot(temp1$s_tol,temp1$r_use)
plot(temp1$r_acqui,temp1$r_use)
par(mfrow=c(1,1))

resource_acquisition1 = temp %>% select(2:22)
colnames(resource_acquisition1) <- c("1","2","3","4","5","6","7","8","9","10","11","12",
                         "13","14","15","16","16","18","19","20","21")
RA_total   = as.data.frame(rbind(rep(max(resource_acquisition1),3),
                                rep(min(resource_acquisition1),3),resource_acquisition1))
radarchart(RA_total)


stress_tolerance1 = temp %>% select(24:27)
colnames(stress_tolerance1) <- c("solute.transport",
                         "solute.synthesis","EPS.biosynthesis(S)",
                         "osmotic.sensors")
st_total   = as.data.frame(rbind(rep(max(stress_tolerance1),3),
                                rep(min(stress_tolerance1),3),stress_tolerance1))
radarchart(st_total)

###############################################################################

# Working with 3 guilds
v                      = cutree(clusterg,k=3)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
mat_trait_f     = as.data.frame(cbind(genome2guild,mat_ori$id,trait_tar))  
temp = mat_trait_f %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat = (as.matrix(temp[2:35]))
heatmap(mat)
# heatmap(mat, Colv = NA, Rowv = NA, scale="column")
# my_colnames2 <- names(temp)

# Preliminary arrangements
a = as.data.frame(colnames(temp))
r_acqui = rowSums(temp %>% select(2:22), na.rm=FALSE)
s_tol   = rowSums(temp %>% select(23:27), na.rm=FALSE)
r_use   = rowSums(temp %>% select(28:35), na.rm=FALSE)

temp1  = as.data.frame(cbind(temp$guild,r_acqui,
                             s_tol,
                             r_use))
mat1 = (as.matrix(temp1[2:4]))
heatmap(mat1)

library(fmsb)

temp2  = as.data.frame(rbind(rep(max(r_acqui),3),
                             rep(min(r_use),3),temp1[2:4]))
radarchart(temp2)


par(mfrow=c(1,3))
plot(temp1$r_acqui,temp1$s_tol)
plot(temp1$s_tol,temp1$r_use)
plot(temp1$r_acqui,temp1$r_use)
par(mfrow=c(1,1))

resource_acquisition1 = temp %>% select(2:22)
colnames(resource_acquisition1) <- c("1","2","3","4","5","6","7","8","9","10","11","12",
                                     "13","14","15","16","16","18","19","20","21")
RA_total   = as.data.frame(rbind(rep(max(resource_acquisition1),3),
                                 rep(min(resource_acquisition1),3),resource_acquisition1))
radarchart(RA_total)


stress_tolerance1 = temp %>% select(24:27)
colnames(stress_tolerance1) <- c("solute.transport",
                                 "solute.synthesis","EPS.biosynthesis(S)",
                                 "osmotic.sensors")
st_total   = as.data.frame(rbind(rep(max(stress_tolerance1),3),
                                 rep(min(stress_tolerance1),3),stress_tolerance1))
radarchart(st_total)


###############################################################################

# Guild 6 (invest the most in resource acquisition)

guild_6 = mat_trait_f %>% filter(guild == 6)

a = as.data.frame(colnames(guild_6))
resource_acquisition1 = rowSums(guild_6 %>% select(3:23), na.rm=FALSE)
stress_tolerance1     = rowSums(guild_6 %>% select(24:28), na.rm=FALSE)
resource_use1         = rowSums(guild_6 %>% select(29:36), na.rm=FALSE)

temp1a  = as.data.frame(cbind(guild_6$guild,resource_acquisition1,
                             resource_use1,
                             stress_tolerance1))
temp2a  = as.data.frame(rbind(rep(max(resource_acquisition),3),
                             rep(min(resource_use),3),temp1a[2:4]))
radarchart(temp2a)

# Stress tolerance
guild_6st = guild_6 %>% select(24:28)
colnames(guild_6st) <- c("EPS.biosynthesis(G)","solute.transport",
                         "solute.synthesis","EPS.biosynthesis(S)",
                         "osmotic.sensors")
temp6st   = as.data.frame(rbind(rep(max(guild_6st),3),
                              rep(min(guild_6st),3),guild_6st))
radarchart(temp6st)

# Resource Acquisition
guild_6ra = guild_6 %>% select(3:23)
colnames(guild_6ra) <- c("1","2","3","4","5","6","7","8","9","10","11","12",
                         "13","14","15","16","16","18","19","20","21")
temp6ra   = as.data.frame(rbind(rep(max(guild_6ra),3),
                                rep(min(guild_6ra),3),guild_6ra))
radarchart(temp6ra)

###############################################################################

# library(GPareto) # Pending
# set.seed(25468)


# Guilds:

guild_1 = mat_trait_f %>% filter(guild == 1)
guild_2 = mat_trait_f %>% filter(guild == 2)
guild_3 = mat_trait_f %>% filter(guild == 3)
guild_4 = mat_trait_f %>% filter(guild == 4)
guild_5 = mat_trait_f %>% filter(guild == 5)
guild_6 = mat_trait_f %>% filter(guild == 6)
guild_7 = mat_trait_f %>% filter(guild == 7)
guild_8 = mat_trait_f %>% filter(guild == 8)
guild_9 = mat_trait_f %>% filter(guild == 9)
guild_10 = mat_trait_f %>% filter(guild == 10)
guild_11 = mat_trait_f %>% filter(guild == 11)
guild_12 = mat_trait_f %>% filter(guild == 12)
guild_13 = mat_trait_f %>% filter(guild == 13)
guild_14 = mat_trait_f %>% filter(guild == 14)
guild_15 = mat_trait_f %>% filter(guild == 15)
guild_16 = mat_trait_f %>% filter(guild == 16)
guild_17 = mat_trait_f %>% filter(guild == 17)
guild_18 = mat_trait_f %>% filter(guild == 18)
guild_19 = mat_trait_f %>% filter(guild == 19)
guild_20 = mat_trait_f %>% filter(guild == 20)
guild_21 = mat_trait_f %>% filter(guild == 21)
guild_22 = mat_trait_f %>% filter(guild == 22)
guild_23 = mat_trait_f %>% filter(guild == 23)
guild_24 = mat_trait_f %>% filter(guild == 24)
guild_25 = mat_trait_f %>% filter(guild == 25)
guild_26 = mat_trait_f %>% filter(guild == 26)
guild_27 = mat_trait_f %>% filter(guild == 27)
guild_28 = mat_trait_f %>% filter(guild == 28)
guild_29 = mat_trait_f %>% filter(guild == 29)
guild_30 = mat_trait_f %>% filter(guild == 30)
guild_31 = mat_trait_f %>% filter(guild == 31)
guild_32 = mat_trait_f %>% filter(guild == 32)
guild_33 = mat_trait_f %>% filter(guild == 33)
guild_34 = mat_trait_f %>% filter(guild == 34)
guild_35 = mat_trait_f %>% filter(guild == 35)
guild_36 = mat_trait_f %>% filter(guild == 36)
guild_37 = mat_trait_f %>% filter(guild == 37)
guild_38 = mat_trait_f %>% filter(guild == 38)
guild_39 = mat_trait_f %>% filter(guild == 39)
guild_40 = mat_trait_f %>% filter(guild == 40)
guild_41 = mat_trait_f %>% filter(guild == 41)
guild_42 = mat_trait_f %>% filter(guild == 42)
guild_43 = mat_trait_f %>% filter(guild == 43)
guild_44 = mat_trait_f %>% filter(guild == 44)
guild_45 = mat_trait_f %>% filter(guild == 45)
guild_46 = mat_trait_f %>% filter(guild == 46)
guild_47 = mat_trait_f %>% filter(guild == 47)
guild_48 = mat_trait_f %>% filter(guild == 48)
guild_49 = mat_trait_f %>% filter(guild == 49)
guild_50 = mat_trait_f %>% filter(guild == 50)

# Working with 200 guilds
v                      = cutree(clusterg,k=200)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
mat_trait_f     = as.data.frame(cbind(genome2guild,trait_tar))  
temp = mat_trait_f %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat = (as.matrix(temp[2:35]))
heatmap(mat, Colv = NA, Rowv = NA, scale="column")
my_colnames2 <- names(temp)
heatmap(mat)
##################################################################v




