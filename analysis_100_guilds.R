# Packages ####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(stats)

# Calling data ####

mat_ori    = read.csv("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/litter_mags_trait_matrixatgranularity.csv",dec=".")
gen_size   = read.delim("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/litter_mags_metadata.txt",dec=".") 

mag_stat   = read.delim("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/mag_stats.txt") 
mag_abun   = read.delim("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/mag_adundance.txt") 

mag_stat   = mag_stat %>% full_join(mag_abun)

# Erase weird MAGs, which were found in preliminary analysis: 116, 175

mat_ori    = mat_ori[-c(116,175),] 
gen_size   = gen_size[-c(116,175),] 
mag_stat   = mag_stat[-c(116,175),]

mat_trait    = mat_ori %>% select(2:191)/gen_size$length
mat_trait    = as.data.frame(cbind(mat_ori$id,mat_trait))
mat_trait    = mat_trait %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))

# Objective 1: Traits responsible for changes in abundance ####

# MAG abundance Analysis  ####

# Grass-ambient vs grass-drought  ####

# Changes in MAGs abundance

mag_stat = mag_stat %>% mutate(ratio1 = log10(mag_stat$Average.1/mag_stat$Average))
mag_stat = mag_stat %>% mutate(pos = ratio1 >= 0)
mag_stat = mag_stat %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))
p<-ggplot(data=mag_stat, aes(x=id_2, y=ratio1, fill = pos)) +
   geom_bar(stat="identity") + theme(legend.position="none") + ylab("Log-transformed Ratio") + 
   xlab("MAG ID")
p

# Correlation plots

mat_trait   = mat_trait %>% mutate(rat_am_dro_gra = mag_stat$ratio1)
mat_trait   = mat_trait %>% mutate(rat_am_dro_gra = ifelse(is.na(rat_am_dro_gra),0,rat_am_dro_gra),
                                   rat_am_dro_gra = ifelse(rat_am_dro_gra==Inf,0,rat_am_dro_gra),
                                   rat_am_dro_gra = ifelse(rat_am_dro_gra ==-Inf, 0, rat_am_dro_gra))
mat_trait_1 = as.data.frame(mat_trait[,-192])
mat = (as.matrix(mat_trait_1[2:192]))

corre = as.data.frame(cor(mat[,1:190], mat[,191]))

# Shrub-ambient vs shrub-drought  ####

# Changes in MAGs abundance

mag_stat.1 = mag_stat %>% mutate(ratio1 = log10(mag_stat$Average.3/mag_stat$Average.2))
mag_stat.1 = mag_stat.1 %>% mutate(pos = ratio1 >= 0)
mag_stat.1 = mag_stat.1 %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))
p<-ggplot(data=mag_stat.1, aes(x=id_2, y=ratio1, fill = pos)) +
  geom_bar(stat="identity") + theme(legend.position="none") + ylab("Log-transformed Ratio") + 
  xlab("MAG ID")
p
# Correlation plots

mat_traits   = mat_trait %>% mutate(rat_am_dro_shr = mag_stat.1$ratio1)
mat_traits   = mat_traits %>% mutate(rat_am_dro_shr = ifelse(is.na(rat_am_dro_shr),0,rat_am_dro_shr),
                                   rat_am_dro_shr = ifelse(rat_am_dro_shr==Inf,0,rat_am_dro_shr),
                                   rat_am_dro_shr = ifelse(rat_am_dro_shr ==-Inf, 0, rat_am_dro_shr))
mat_trait_1s = as.data.frame(mat_traits[,-192])
mat_trait_1s = as.data.frame(mat_trait_1s[,-192])
mats = (as.matrix(mat_trait_1s[2:192]))

corre1 = as.data.frame(cor(mats[,1:190], mats[,191]))

# Grass-ambient vs Shrubland-ambient  ####

# Changes in MAGs abundance

mag_stat.2 = mag_stat %>% mutate(ratio1 = log10(mag_stat$Average.2/mag_stat$Average))
mag_stat.2 = mag_stat.2 %>% mutate(pos = ratio1 >= 0)
mag_stat.2 = mag_stat.2 %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))
p<-ggplot(data=mag_stat.2, aes(x=id_2, y=ratio1, fill = pos)) +
  geom_bar(stat="identity") + theme(legend.position="none") + ylab("Log-transformed Ratio") + 
  xlab("MAG ID")
p
# Correlation plots

mat_traitsg   = mat_trait %>% mutate(rat_shr_gras = mag_stat.2$ratio1)
mat_traitsg   = mat_traitsg %>% mutate(rat_shr_gras = ifelse(is.na(rat_shr_gras),0,rat_shr_gras),
                                     rat_shr_gras = ifelse(rat_shr_gras==Inf,0,rat_shr_gras),
                                     rat_shr_gras = ifelse(rat_shr_gras ==-Inf, 0, rat_shr_gras))
mat_trait_1gs = as.data.frame(mat_traitsg[,-192])
mat_trait_1gs = as.data.frame(mat_trait_1gs[,-192])
matgs = (as.matrix(mat_trait_1gs[2:192]))

corre2 = as.data.frame(cor(matgs[,1:190], matgs[,191]))

# Objective 2: Grouping based on principal traits  ####






# Objective 3: Traits-tradeoff analysis  ####





###############################################################################
###############################################################################
###############################################################################




###############################################################################
# Normalized data
###############################################################################

mat_trait    = mat_ori %>% select(2:191)/gen_size$length
mat_trait    = as.data.frame(cbind(mat_ori$id,mat_trait))
trait_sd     = mat_trait %>% select(1,43:61) # substrate degradation
trait_tp     = mat_trait %>% select(1:42) # transporters or substrate uptake
trait_st     = mat_trait %>% select(1,162:190) # stress related traits

###############################################################################
# Exploratory plots - gene costs (normalized genes) #
###############################################################################

# complex compounds
boxplot(trait_sd%>% select(2:7,8),ylab="gene cost",
        names=c("cellulose","chitin", "heteromannan", "linkage-glucan", "xylan",
                "xyloglucan","protein"))
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

# Resource - Use
boxplot(mat_trait%>% select(93:99,125),ylab="gene cost",
        names=c("oxi.pentose","nonoxi.pentose","glyoxylate.cycle",
                "TCA.cycle","ETC.complex.I","ETC.complex.II","ETC.complex.III",
                "nitrite.oxidation"))

###############################################################################
# Targeted traits for DEMENT modeling
###############################################################################

trait_tar = mat_trait %>% select(4,6:8,18,23:25,31:35,38,43:49,174:177,93:99,125) # 
genomes   = as.list(mat_trait[1])
test      = trait_tar %>% summarize_if(is.numeric, sum, na.rm=TRUE)

###############################################################################
# Forming functional guilds 
# (https://www.youtube.com/watch?v=GPOUGpF-Sno)
# https://github.com/ukaraoz/microtrait
###############################################################################

library(vegan)
set.seed(1)

# Gene cost
###############################################################################
distanceg  = vegdist(trait_tar, method = "bray", binary = FALSE)
clusterg   = hclust(distanceg, method="ward")
nguilds    = seq(2, nrow(trait_tar), 2)
plot(clusterg)

# Similiarity within guilds
###############################################################################
my_vec = c()

for(i in nguilds) {
  v                      = cutree(clusterg,k=i)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_2               = vegan::adonis2(distanceg ~ guild, data = genome2guild, perm = 1)
  temp                   = adonis_2$R2
  temp_1                 = temp[1]
  my_out                 = temp_1        
  my_vec <- c(my_vec, my_out)    
}

# Plotting
mat_r2_c  = as.data.frame(cbind(nguilds,my_vec))
ggplot(data=mat_r2_c,aes(x=nguilds,y=my_vec)) + geom_line() +
  xlab("# of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 20)) + 
  geom_hline(yintercept=0.8349532, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 100, linetype="dashed", color = "red")

###############################################################################
# Finding the optimal number of guilds
# https://www.r-bloggers.com/2019/01/10-tips-for-choosing-the-optimal-number-of-clusters/
# https://r-tastic.co.uk/post/optimal-number-of-clusters/
###############################################################################

# Manual Test 
# Silhouette width : average dissimilarity between and object and all other 
# objects within a cluster to which it belongs
# https://www.youtube.com/watch?v=ezn6bFz-Phk&t=2s (19:00)
###############################################################################
clusterg_2 = clusterg
Si = numeric(nrow(trait_tar))
for (k in 2:(nrow(trait_tar) - 1))
{
  sil = silhouette(cutree(clusterg_2,k=k),distanceg)
  Si[k] = summary(sil)$avg.width
}

k.best = which.max(Si)
plot(1:nrow(trait_tar),Si,type="h",main="Silhouette-optimal number of clusters",
     xlab = "Number of clusters",ylab = "Average silhouette width")
axis(1,k.best,paste("optimum",k.best,sep="\n"),col="red",font=2,col.axis="red")
points(k.best,max(Si),pch=16,col="red",cex=1.5)

dend = as.dendrogram(clusterg_2)
heatmap(as.matrix(distanceg),Rowv = dend, symm = T, margin = c(3,3))

or = tabasco(trait_tar, scale="log",dend)
or = tabasco(trait_tar,dend)

# CALINSKY CRITERION
###############################################################################
cal_fit2 <- cascadeKM(trait_tar, 1, 10, iter = 1000)
plot(cal_fit2, sortg = TRUE, grpmts.plot = TRUE)
calinski.best2 <- as.numeric(which.max(cal_fit2$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best2, "\n")

# Elbow method
###############################################################################
library(factoextra)
# install.packages("NbClust")
library(NbClust)

fviz_nbclust(trait_tar, kmeans, method = "wss", k.max = 10) + 
  theme_minimal() + ggtitle("the Elbow Method")

# Gap Statistic
###############################################################################
library(cluster)
gap_stat <- clusGap(trait_tar, kmeans, 100, B = 100, verbose = interactive())
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")

# silhouette
###############################################################################
fviz_nbclust(trait_tar, kmeans, method = "silhouette", k.max = 100) + 
  theme_minimal() + ggtitle("The Silhouette Plot")

###############################################################################
# Working with 3 guilds
###############################################################################

v                      = cutree(clusterg,k=4) # Clusters
genome2guild_5         = data.frame(guild = factor(v))
rownames(genome2guild_5) = names(v)
mat_trait_5     = as.data.frame(cbind(genome2guild_5,mat_ori$id,trait_tar))  
temp = mat_trait_5 %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat = (as.matrix(temp[2:34]))
a   = as.data.frame(colnames(mat))
colnames(mat) = c("amino.sug.trans","glycosi.trans",
                  "FOS.fayt.trans","monosac.trans","aminoac.trans",
                  "lipid.trans","amide.trans","NH4.trans","nucleob.trans",
                  "nucleos.trans","nucleot.trans","ribonucle.trans",
                  "organophos.trans","peptide.trans","cellulose","chitin",
                  "heteroman","glucan.mix","xylan","xyloglucan","protein",
                  "sol.trans","sol.synt","EPS.synt","osmotic.sensor",
                  "oxi.pentose","nonoxi.pentose","glyoxylate","TCA","ETC.I",
                  "ETC.II","ETC.III","nitrite")

# Ordination
###############################################################################
set.seed(2)
colnames(trait_tar) = c("amino.sug.trans","glycosi.trans",
                        "FOS.fayt.trans","monosac.trans","aminoac.trans",
                        "lipid.trans","amide.trans","NH4.trans","nucleob.trans",
                        "nucleos.trans","nucleot.trans","ribonucle.trans",
                        "organophos.trans","peptide.trans","cellulose","chitin",
                        "heteroman","glucan.mix","xylan","xyloglucan","protein",
                        "sol.trans","sol.synt","EPS.synt","osmotic.sensor",
                        "oxi.pentose","nonoxi.pentose","glyoxylate","TCA","ETC.I",
                        "ETC.II","ETC.III","nitrite")
temp_1 = metaMDS(trait_tar,autotransform = F,trymax = 2000)
ordiplot(temp_1)
ordiplot(temp_1,type = "t")
library(ggvegan)

# One Panel
###############################################################################
fort = fortify(temp_1)

ggplot() + geom_point(data=subset(fort,Score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =mat_trait_5$guild,size = 5),
                      alpha=0.5) + 
  geom_segment(data=subset(fort,Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort,Score=="species"),
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(temp_1$stress,2)))

# Substrate
###############################################################################
substrate = as.data.frame(substring(mat_ori$id,1,1)) 
colnames(substrate) = c("subst")

ggplot() + geom_point(data=subset(fort,Score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =mat_trait_5$guild,
                                    shape = substrate$subst,size = 5),alpha=0.5) + 
  geom_segment(data=subset(fort,Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort,Score=="species"),
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(nmds$stress,2)))

# Other codes
###############################################################################
other = as.data.frame(substring(mat_ori$id,20)) 
colnames(other) = c("other")
other = other %>% mutate(other = replace(other,other == ".orig","orig"))
other = other %>% mutate(other = replace(other,other == ".permissive","permissive"))
other = other %>% mutate(other = replace(other,other == ".strict","strict"))

ggplot() + geom_point(data=subset(fort,Score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =mat_trait_5$guild,
                                    shape = other$other,size = 5),alpha=0.5) + 
  geom_segment(data=subset(fort,Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort,Score=="species"),
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(nmds$stress,2)))

# Spider plots for YAS framework traits
###############################################################################
b   = as.data.frame(colnames(temp)) # extract columns to extract the data for each category
r_acqui = rowSums(temp %>% select(2:22), na.rm=FALSE)
s_tol   = rowSums(temp %>% select(23:26), na.rm=FALSE)
r_use   = rowSums(temp %>% select(27:34), na.rm=FALSE)

library(fmsb)

temp1  = as.data.frame(cbind(temp$guild,r_acqui,s_tol,r_use))
total  = rowSums(temp1 %>% select(2:4), na.rm=FALSE)
temp2  = as.data.frame(rbind(rep(max(total),3),
                             rep(min(r_use),3),temp1[2:4]))
colnames(temp2) = c("Resource Acquisition","Stress Tolerance",
                    "Resource Use")
radarchart(temp2) # Spider plots using 100 guilds and total costs

temp1a  = as.data.frame(cbind(temp$guild,r_acqui/21,s_tol/4,r_use/8))
totala  = rowSums(temp1 %>% select(2:4), na.rm=FALSE)
temp2a  = as.data.frame(rbind(rep(max(temp1a$V2),3),
                              rep(min(temp1a$V4),3),temp1a[2:4]))
colnames(temp2a) = c("Resource Acquisition","Stress Tolerance",
                     "Resource Use")
radarchart(temp2a) # Spider plots using 100 guilds and costs per trait

# Binary plots scaled by trait number
###############################################################################
r_acqui1 = rowSums(temp %>% select(2:22), na.rm=FALSE)/21
s_tol1   = rowSums(temp %>% select(23:26), na.rm=FALSE)/4
r_use1   = rowSums(temp %>% select(27:34), na.rm=FALSE)/8

par(mfrow=c(1,3))
plot(r_acqui1,s_tol1,xlab = "Resource Acquisition", ylab = "Stress Tolerance",
     col = "red",pch = 15,cex.lab = 1.5,cex = 3)
plot(s_tol1,r_use1,xlab = "Stress Tolerance", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5,cex = 3)
plot(r_acqui1,r_use1,xlab = "Resource Acquisition", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5,cex = 3)
par(mfrow=c(1,1))


############

###############################################################################
# Working with 100 guilds
###############################################################################

v                      = cutree(clusterg,k=100) # Clusters
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
mat_trait_f     = as.data.frame(cbind(genome2guild,mat_ori$id,trait_tar))  
temp = mat_trait_f %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat = (as.matrix(temp[2:34]))
a   = as.data.frame(colnames(mat))
colnames(mat) = c("amino.sug.trans","glycosi.trans",
                  "FOS.fayt.trans","monosac.trans","aminoac.trans",
                  "lipid.trans","amide.trans","NH4.trans","nucleob.trans",
                  "nucleos.trans","nucleot.trans","ribonucle.trans",
                  "organophos.trans","peptide.trans","cellulose","chitin",
                  "heteroman","glucan.mix","xylan","xyloglucan","protein",
                  "sol.trans","sol.synt","EPS.synt","osmotic.sensor",
                  "oxi.pentose","nonoxi.pentose","glyoxylate","TCA","ETC.I",
                  "ETC.II","ETC.III","nitrite")

# Heatmap of trait distribution
###############################################################################
heatmap(mat,Colv = NA) # including main cluster
heatmap(mat) # including both clusters

###############################################################################
# Ordination
###############################################################################
nmds = metaMDS(distanceg,trymax = 2000)
nmds
stressplot(nmds)
nmds$points

NMDS1 <- nmds$points[,1] 
NMDS2 <- nmds$points[,2]
df_nor2.plot<-cbind(mat_trait_f, NMDS1, NMDS2)

ggplot(data=df_nor2.plot,aes(x=NMDS1,y=NMDS2,color =guild))+geom_point(size = 4) + 
  annotate("text", x=-0.4, y=-0.6, label=paste('Stress =',round(nmds$stress,2)))

###############################################################################
# Binary plots with total values
###############################################################################
b   = as.data.frame(colnames(temp)) # extract columns to extract the data for each cathegory
r_acqui = rowSums(temp %>% select(2:22), na.rm=FALSE)
s_tol   = rowSums(temp %>% select(23:26), na.rm=FALSE)
r_use   = rowSums(temp %>% select(27:34), na.rm=FALSE)

par(mfrow=c(1,3))
plot(r_acqui,s_tol,xlab = "Resource Acquisition", ylab = "Stress Tolerance",
     col = "red",pch = 15,cex.lab = 1.5)
plot(s_tol,r_use,xlab = "Stress Tolerance", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5)
plot(r_acqui,r_use,xlab = "Resource Acquisition", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5)
par(mfrow=c(1,1))

###############################################################################
# Binary plots scaled by trait number
###############################################################################
r_acqui1 = rowSums(temp %>% select(2:22), na.rm=FALSE)/21
s_tol1   = rowSums(temp %>% select(23:26), na.rm=FALSE)/4
r_use1   = rowSums(temp %>% select(27:34), na.rm=FALSE)/8

par(mfrow=c(1,3))
plot(r_acqui1,s_tol1,xlab = "Resource Acquisition", ylab = "Stress Tolerance",
     col = "red",pch = 15,cex.lab = 1.5)
plot(s_tol1,r_use1,xlab = "Stress Tolerance", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5)
plot(r_acqui1,r_use1,xlab = "Resource Acquisition", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5)
par(mfrow=c(1,1))

###############################################################################
# Spider plots for YAS framework traits
###############################################################################

library(fmsb)

temp1  = as.data.frame(cbind(temp$guild,r_acqui,s_tol,r_use))
total  = rowSums(temp1 %>% select(2:4), na.rm=FALSE)
temp2  = as.data.frame(rbind(rep(max(total),3),
                             rep(min(r_use),3),temp1[2:4]))
colnames(temp2) = c("Resource Acquisition","Stress Tolerance",
                    "Resource Use")
radarchart(temp2) # Spider plots using 100 guilds and total costs

temp1a  = as.data.frame(cbind(temp$guild,r_acqui/21,s_tol/4,r_use/8))
totala  = rowSums(temp1 %>% select(2:4), na.rm=FALSE)
temp2a  = as.data.frame(rbind(rep(max(r_acqui/21),3),
                             rep(min(r_use),3),temp1a[2:4]))
colnames(temp2a) = c("Resource Acquisition","Stress Tolerance",
                    "Resource Use")
radarchart(temp2a) # Spider plots using 100 guilds and costs per trait

###############################################################################
# Spider plots for within YAS framework traits
###############################################################################

# Resource acquisition
###############################################################################
resource_acquisition1 = temp %>% select(2:22)
totala  = rowSums(resource_acquisition1, na.rm=FALSE)
colnames(resource_acquisition1) <- c("1","2","3","4","5","6","7","8","9","10","11","12",
                                     "13","14","15","16","17","18","19","20","21")
RA_total   = as.data.frame(rbind(rep(max(totala),3),
                                 rep(min(resource_acquisition1),3),resource_acquisition1))
radarchart(RA_total)

par(mfrow=c(1,2))
plot(temp$Resource.Acquisition.Substrate.degradation.simple.compound.degradation.protein.degradation,
     temp$Resource.Acquisition.Substrate.uptake.free.amino.acids.transport,xlab = "protein degradation", 
     ylab = "free aminoacids transporters",
     col = "red",pch = 15,cex.lab = 1.5)
plot(temp$Resource.Acquisition.Substrate.degradation.simple.compound.degradation.protein.degradation,
     temp$Resource.Acquisition.Substrate.degradation.complex.carbohydrate.depolymerization.xylan.and.heteroxylan.breakdown,xlab = "protein degradation", 
     ylab = "xylan-heteroxylan degradation",
     col = "red",pch = 15,cex.lab = 1.5)
par(mfrow=c(1,1))

# Stress Tolerance
###############################################################################
stress_tolerance1 = temp %>% select(23:26)
totala  = rowSums(stress_tolerance1, na.rm=FALSE)
colnames(stress_tolerance1) <- c("solute.transport",
                                 "solute.synthesis","EPS.biosynthesis(S)",
                                 "osmotic.sensors")
st_total   = as.data.frame(rbind(rep(max(totala),3),
                                 rep(min(stress_tolerance1),3),stress_tolerance1))
radarchart(st_total)

par(mfrow=c(1,2))
plot(temp$Stress.Tolerance.Specific.desiccation.osmotic.salt.stress.accumulation.of.compatible.solutes.synthesis.of.compatible.solutes,
     temp$Stress.Tolerance.Specific.desiccation.osmotic.salt.stress.accumulation.of.compatible.solutes.transport.of.compatible.solutes,xlab = "solute synthesis", 
     ylab = "solute transport",
     col = "red",pch = 15,cex.lab = 1.5)
plot(temp$Stress.Tolerance.Specific.desiccation.osmotic.salt.stress.accumulation.of.compatible.solutes.synthesis.of.compatible.solutes,
     temp$Stress.Tolerance.Specific.desiccation.osmotic.salt.stress.enhancement.of.microenvironment.EPS.biosynthesis.export,xlab = "solute synthesis", 
     ylab = "EPS biosynthesis",
     col = "red",pch = 15,cex.lab = 1.5)
par(mfrow=c(1,1))

# Resource use
###############################################################################
resource_use1 = temp %>% select(27:34)
totala  = rowSums(resource_use1, na.rm=FALSE)
colnames(resource_use1) <- c("1","2","3","4","5","6","7","8")
st_total   = as.data.frame(rbind(rep(max(totala),3),
                                 rep(min(resource_use1),3),resource_use1))
radarchart(st_total)

par(mfrow=c(1,2))
plot(temp$Resource.Use.Chemotrophy.chemoorganoheterotrophy.aerobic.respiration.electron.transport.chain..ETC.complex.III,
     temp$Resource.Use.Chemotrophy.chemoorganoheterotrophy.aerobic.respiration.electron.transport.chain..ETC.complex.II,xlab = "ETC complex III", 
     ylab = "ETC complex II",
     col = "red",pch = 15,cex.lab = 1.5)
plot(temp$Resource.Use.Chemotrophy.chemoorganoheterotrophy.aerobic.respiration.electron.transport.chain..ETC.complex.III,
     temp$Resource.Use.Chemotrophy.chemoorganoheterotrophy.aerobic.respiration.oxidative.pentose.phosphate.pathway,xlab = "ETC complex III", 
     ylab = "oxidative pentose phosphate pathway ",
     col = "red",pch = 15,cex.lab = 1.5)
par(mfrow=c(1,1))

###############################################################################
# Guild 65 (invest the most in resource acquisition)
###############################################################################

guild_65 = mat_trait_f %>% filter(guild == 65)

a = as.data.frame(colnames(guild_65))
resource_acquisition1 = rowSums(guild_65 %>% select(3:23), na.rm=FALSE)
stress_tolerance1     = rowSums(guild_65 %>% select(24:27), na.rm=FALSE)
resource_use1         = rowSums(guild_65 %>% select(28:35), na.rm=FALSE)

temp1a  = as.data.frame(cbind(guild_65$guild,resource_acquisition1,
                              resource_use1,
                              stress_tolerance1))
temp2a  = as.data.frame(rbind(rep(max(resource_acquisition1),3),
                              rep(min(resource_use1),3),temp1a[2:4]))
radarchart(temp1a)

# Stress tolerance
###############################################################################
guild_6st = guild_65 %>% select(24:27)
colnames(guild_6st) <- c("solute.transport",
                         "solute.synthesis","EPS.biosynthesis(S)",
                         "osmotic.sensors")
temp6st   = as.data.frame(rbind(rep(max(guild_6st),3),
                                rep(min(guild_6st),3),guild_6st))
radarchart(temp6st)

# Resource Acquisition
###############################################################################
guild_6ra = guild_65 %>% select(3:23)
colnames(guild_6ra) <- c("1","2","3","4","5","6","7","8","9","10","11","12",
                         "13","14","15","16","16","18","19","20","21")
temp6ra   = as.data.frame(rbind(rep(max(guild_6ra),3),
                                rep(min(guild_6ra),3),guild_6ra))
radarchart(temp6ra)

# Resource use
###############################################################################
guild_6ru = guild_65 %>% select(28:35)
colnames(guild_6ru) <- c("1","2","3","4","5","6","7","8")
temp6ru   = as.data.frame(rbind(rep(max(guild_6ru),3),
                                 rep(min(guild_6ru),3),guild_6ru))
radarchart(temp6ru)

###############################################################################
# Correlation plots
###############################################################################

library(ggcorrplot)

# Correlation plots - all samples
###############################################################################

colnames(trait_tar) = c("amino.sug.trans","glycosi.trans",
                        "FOS.fayt.trans","monosac.trans","aminoac.trans",
                        "lipid.trans","amide.trans","NH4.trans","nucleob.trans",
                        "nucleos.trans","nucleot.trans","ribonucle.trans",
                        "organophos.trans","peptide.trans","cellulose","chitin",
                        "heteroman","glucan.mix","xylan","xyloglucan","protein",
                        "sol.trans","sol.synt","EPS.synt","osmotic.sensor",
                        "oxi.pentose","nonoxi.pentose","glyoxylate","TCA","ETC.I",
                        "ETC.II","ETC.III","nitrite")
corr <- round(cor(trait_tar), 1)
p.mat <- cor_pmat(trait_tar)
ggcorrplot(corr)
# using hierarchical clustering
ggcorrplot(corr, hc.order = TRUE, outline.col = "white")
# Leave blank on no significant coefficient
ggcorrplot(corr, p.mat = p.mat, hc.order = TRUE,
           type = "lower", insig = "blank")

# Correlation plots - 100 guilds
###############################################################################

temp_100 = temp %>% select(2:34)
colnames(temp_100) = c("amino.sug.trans","glycosi.trans",
                        "FOS.fayt.trans","monosac.trans","aminoac.trans",
                        "lipid.trans","amide.trans","NH4.trans","nucleob.trans",
                        "nucleos.trans","nucleot.trans","ribonucle.trans",
                        "organophos.trans","peptide.trans","cellulose","chitin",
                        "heteroman","glucan.mix","xylan","xyloglucan","protein",
                        "sol.trans","sol.synt","EPS.synt","osmotic.sensor",
                        "oxi.pentose","nonoxi.pentose","glyoxylate","TCA","ETC.I",
                        "ETC.II","ETC.III","nitrite")
corr100 <- round(cor(temp_100), 1)
p.mat100 <- cor_pmat(temp_100)
ggcorrplot(corr100)
ggcorrplot(corr100, hc.order = TRUE, outline.col = "white")
# Leave blank on no significant coefficient
ggcorrplot(corr100, p.mat = p.mat100, hc.order = TRUE,
           type = "lower", insig = "blank")

###############################################################################
# CAZymes-related traits
###############################################################################

trait_tarCAZ = mat_trait %>% select(6,8,43:48) # 
genomes      = as.list(mat_trait[1])
test         = trait_tar %>% summarize_if(is.numeric, sum, na.rm=TRUE)
names(trait_tarCAZ)
colnames(trait_tarCAZ) = c("glycoside.trans","monosac.trans","cellulose","chitin",
                           "heteromannan","mixed.glucan","xylan","xyloglucan")
boxplot(trait_tarCAZ,ylab="gene cost")

###############################################################################
# Forming functional guilds 
###############################################################################

distanceg  = vegdist(trait_tarCAZ, method = "bray", binary = FALSE)
clusterg   = hclust(distanceg, method="ward")
nguilds    = seq(2, nrow(trait_tar), 2)
plot(clusterg)

# Similiarity within guilds
###############################################################################
my_vec = c()

for(i in nguilds) {
  v                      = cutree(clusterg,k=i)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_2               = vegan::adonis2(distanceg ~ guild, data = genome2guild, perm = 1)
  temp                   = adonis_2$R2
  temp_1                 = temp[1]
  my_out                 = temp_1        
  my_vec <- c(my_vec, my_out)    
}

# Plotting
mat_r2_c  = as.data.frame(cbind(nguilds,my_vec))
ggplot(data=mat_r2_c,aes(x=nguilds,y=my_vec)) + geom_line() +
  xlab("# of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 20)) + 
  geom_hline(yintercept=0.9200740, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 100, linetype="dashed", color = "red")

barplot(c(29,67,4),ylab="Number of guilds",names.arg=c("Grass","Both","Shrub"))