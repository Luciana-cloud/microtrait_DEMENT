###############################################################################
# Packages #
###############################################################################

library(dplyr)
library(tidyverse)
library(ggplot2)

###############################################################################
# Calling data #
###############################################################################

mat_ori    = read.csv("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/litter_mags_trait_matrixatgranularity.csv",dec=".")
gen_size   = read.delim("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/litter_mags_metadata.txt",dec=".") 

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

###############################################################################
# Targeted traits
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
  geom_hline(yintercept=0.8424206, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 100, linetype="dashed", color = "red")

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
# Test
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

#install.packages("remotes")
#remotes::install_github("gavinsimpson/ggvegan")
#autoplot(temp_1)

# full control with fortified ordination
###############################################################################
fort = fortify(temp_1)

# One Panel
###############################################################################

ggplot() + geom_point(data=subset(fort,Score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =mat_trait_f$guild),
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
        axis.line=element_line(colour="black"))

# Two Panels
###############################################################################

p1 = ggplot() + geom_point(data=subset(fort,Score=="sites"),
                           mapping = aes(x=NMDS1,y=NMDS2),
                           colour="black",
                           alpha=0.5) + 
  geom_segment(data=subset(fort,Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0,
               alpha=0) + 
  geom_text(data=subset(fort,Score=="species"),
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1),alpha=0) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"))

p2 = ggplot() + geom_point(data=subset(fort,Score=="sites"),
                           mapping = aes(x=NMDS1,y=NMDS2,color =mat_trait_f$guild),
                           colour="black",
                           alpha=0) + 
  geom_segment(data=subset(fort,Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0) + 
  geom_text(data=subset(fort,Score=="species"),
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + scale_x_continuous(limit = c(-0.25,0.25))

p3 = ggplot() + geom_point(data=subset(fort,Score=="sites"),
                           mapping = aes(x=NMDS1,y=NMDS2,color =mat_trait_f$guild),
                           alpha=0.5) + 
  geom_segment(data=subset(fort,Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0,
               alpha=0) + 
  geom_text(data=subset(fort,Score=="species"),
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1),alpha=0) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"))

# Goodness of fit
###############################################################################

gof <- goodness(temp_1)
plot(temp_1, type="t", main = "goodness of fit")
points(temp_1, display="sites", cex=gof*100)


en = envfit(temp_1, trait_tar, permutations = 999, na.rm = TRUE)
en

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
